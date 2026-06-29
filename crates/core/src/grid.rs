//! Gridded PV potential over a DEM (feature `terrain`).
//!
//! Builds on the validated point chain ([`crate::solpos`], [`crate::irradiance`],
//! [`crate::pv`]) and the SurtGIS terrain layer: each DEM cell is treated as a
//! ground-following surface whose tilt and azimuth come from the terrain
//! `slope`/`aspect`, with beam shading from the surrounding topography via
//! `horizon_angles`.
//!
//! For a self-contained "clear-sky PV potential" map, global horizontal
//! irradiance is generated with the Haurwitz clear-sky model and decomposed with
//! Erbs; supplying measured / TMY series is left to v0.2.
//!
//! The scene is assumed small enough that solar geometry at the scene centre
//! applies across the grid (the usual local-DEM approximation); per-cell
//! latitude is a v0.2 refinement.

use rayon::prelude::*;
use surtgis_algorithms::terrain::{
    aspect, horizon_angles, slope, AspectOutput, HorizonParams, SlopeParams, SlopeUnits,
};
use surtgis_core::Raster;

use crate::error::{Error, Result};
use crate::irradiance::{
    erbs, extra_radiation, haurwitz_clearsky_ghi, poa_irradiance, relative_airmass, Decomposition,
    SkyModel,
};
use crate::pv::{ac_power, PvSystem};
use crate::solpos::{solar_position, DateTimeUtc, Location};

const DEG: f64 = std::f64::consts::PI / 180.0;

/// Configuration for a gridded PV-potential run over a single day.
#[derive(Debug, Clone)]
pub struct GridConfig {
    /// Scene-centre location used for solar geometry.
    pub center: Location,
    /// Calendar day to integrate (the time-of-day fields are ignored).
    pub date: DateTimeUtc,
    /// Integration time step in minutes (e.g. 15 or 30).
    pub time_step_minutes: u32,
    /// Sky-diffuse transposition model.
    pub sky_model: SkyModel,
    /// Ground albedo (e.g. 0.25 desert, 0.2 vegetation, 0.7 snow).
    pub albedo: f64,
    /// PV system whose yield is mapped per cell.
    pub system: PvSystem,
    /// Ambient air temperature (°C), uniform over the scene for v0.1.
    pub temp_air: f64,
    /// Wind speed (m/s), uniform over the scene for v0.1.
    pub wind: f64,
    /// Horizon-angle parameters (search radius in cells, number of directions).
    pub horizon: HorizonParams,
}

impl GridConfig {
    /// Sensible defaults for a small DEM: 15-minute steps, Perez sky, desert
    /// albedo, the 1 kW reference system, 18 °C / 2 m·s⁻¹.
    pub fn new(center: Location, date: DateTimeUtc) -> Self {
        Self {
            center,
            date,
            time_step_minutes: 15,
            sky_model: SkyModel::Perez,
            albedo: 0.25,
            system: PvSystem::reference_1kw(),
            temp_air: 18.0,
            wind: 2.0,
            horizon: HorizonParams::default(),
        }
    }
}

/// Precomputed solar state for one daylight time step (scene-centre geometry
/// plus the clear-sky decomposition that is identical for every cell).
struct SunStep {
    zenith: f64,
    azimuth: f64,
    elev_rad: f64,
    az_rad: f64,
    airmass: f64,
    base: Decomposition,
}

/// Per-cell results of a gridded PV-potential run.
pub struct GridResult {
    /// Plane-of-array insolation, Wh/m² over the day.
    pub poa_wh: Raster<f64>,
    /// AC energy for the configured system, Wh over the day.
    pub ac_wh: Raster<f64>,
    /// Specific yield, kWh/kWp over the day (= AC energy ÷ nameplate),
    /// the orientation-independent PVGIS metric.
    pub specific_yield: Raster<f64>,
    /// Number of daylight time steps integrated.
    pub sun_steps: usize,
}

/// Compute the daily clear-sky PV potential over a DEM.
///
/// Returns per-cell POA insolation, AC energy and specific yield. Reuses SurtGIS
/// `slope`/`aspect`/`horizon_angles`; the per-cell physics is the same validated
/// point chain used elsewhere in the crate.
pub fn pv_potential(dem: &Raster<f64>, cfg: &GridConfig) -> Result<GridResult> {
    let (rows, cols) = dem.shape();

    // Terrain orientation and skyline from SurtGIS.
    let slope_rad = slope(dem, SlopeParams { units: SlopeUnits::Radians, ..Default::default() })
        .map_err(|e| Error::Terrain(format!("slope: {e}")))?;
    let asp_deg = aspect(dem, AspectOutput::Degrees)
        .map_err(|e| Error::Terrain(format!("aspect: {e}")))?;
    let horizon = horizon_angles(dem, cfg.horizon.clone())
        .map_err(|e| Error::Terrain(format!("horizon_angles: {e}")))?;

    let dt_hours = cfg.time_step_minutes as f64 / 60.0;
    let steps = (24 * 60) / cfg.time_step_minutes.max(1);
    let doy = cfg.date.day_of_year();
    let dni_extra = extra_radiation(doy);

    // Precompute the solar state at the scene centre for every daylight step,
    // so the per-cell loop never recomputes geometry or decomposition.
    let mut sun_steps_vec: Vec<SunStep> = Vec::new();
    for i in 0..steps {
        let total_min = i * cfg.time_step_minutes;
        let when = DateTimeUtc::new(
            cfg.date.year,
            cfg.date.month,
            cfg.date.day,
            total_min / 60,
            total_min % 60,
            0,
        )?;
        let sun = solar_position(when, cfg.center);
        if sun.apparent_elevation <= 0.0 {
            continue; // night
        }
        let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
        sun_steps_vec.push(SunStep {
            zenith: sun.apparent_zenith,
            azimuth: sun.azimuth,
            elev_rad: sun.apparent_elevation * DEG,
            az_rad: sun.azimuth * DEG,
            airmass: relative_airmass(sun.apparent_zenith),
            base: erbs(ghi, sun.apparent_zenith, doy),
        });
    }
    let sun_steps = sun_steps_vec.len();

    // Per-cell accumulation, parallelised over the flattened grid. Each cell is
    // independent; slope/aspect/horizon are read-only and Sync.
    let energies: Vec<(f64, f64)> = (0..rows * cols)
        .into_par_iter()
        .map(|idx| {
            let (r, c) = (idx / cols, idx % cols);

            // Border cells (Horn's method undefined) come back as NaN → flat.
            let raw_tilt = slope_rad.get(r, c).unwrap_or(0.0);
            let tilt = if raw_tilt.is_finite() { raw_tilt.to_degrees() } else { 0.0 };
            let a = asp_deg.get(r, c).unwrap_or(-1.0);
            let surface_azimuth = if a.is_finite() && a >= 0.0 { a } else { 0.0 };

            let mut poa_acc = 0.0;
            let mut ac_acc = 0.0;
            for s in &sun_steps_vec {
                // Beam shading: cell is in cast shadow when the sun sits below
                // the terrain horizon at its azimuth.
                let h = horizon.interpolate(r, c, s.az_rad);
                let decomp = if s.elev_rad < h {
                    Decomposition { ghi: s.base.dhi, dni: 0.0, dhi: s.base.dhi }
                } else {
                    s.base
                };
                let poa = poa_irradiance(
                    cfg.sky_model,
                    tilt,
                    surface_azimuth,
                    decomp,
                    s.zenith,
                    s.azimuth,
                    cfg.albedo,
                    dni_extra,
                    s.airmass,
                );
                poa_acc += poa.global * dt_hours;
                ac_acc += ac_power(&cfg.system, poa.global, cfg.temp_air, cfg.wind) * dt_hours;
            }
            (poa_acc, ac_acc)
        })
        .collect();

    let poa_vec: Vec<f64> = energies.iter().map(|e| e.0).collect();
    let ac_vec: Vec<f64> = energies.iter().map(|e| e.1).collect();
    // Specific yield = AC energy (Wh) / nameplate (W) = kWh/kWp.
    let sy_vec: Vec<f64> = ac_vec.iter().map(|&ac| ac / cfg.system.pdc0).collect();

    let to_raster = |v: Vec<f64>| -> Raster<f64> {
        let mut out = Raster::from_vec(v, rows, cols).expect("dimensions match dem");
        out.set_transform(*dem.transform());
        out
    };

    Ok(GridResult {
        poa_wh: to_raster(poa_vec),
        ac_wh: to_raster(ac_vec),
        specific_yield: to_raster(sy_vec),
        sun_steps,
    })
}
