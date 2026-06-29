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

    let mut poa_wh = dem.like(0.0);
    let mut ac_wh = dem.like(0.0);

    let dt_hours = cfg.time_step_minutes as f64 / 60.0;
    let steps = (24 * 60) / cfg.time_step_minutes.max(1);
    let doy = cfg.date.day_of_year();
    let dni_extra = extra_radiation(doy);
    let mut sun_steps = 0usize;

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
        sun_steps += 1;

        let airmass = relative_airmass(sun.apparent_zenith);
        let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
        let base = erbs(ghi, sun.apparent_zenith, doy);
        let sun_elev_rad = sun.apparent_elevation * DEG;
        let sun_az_rad = sun.azimuth * DEG;

        for r in 0..rows {
            for c in 0..cols {
                // Border cells (Horn's method undefined) come back as NaN;
                // treat them as flat.
                let raw_tilt = slope_rad.get(r, c).unwrap_or(0.0);
                let tilt = if raw_tilt.is_finite() { raw_tilt.to_degrees() } else { 0.0 };
                // Flat cells are tagged -1 by aspect(); their azimuth is moot.
                let a = asp_deg.get(r, c).unwrap_or(-1.0);
                let surface_azimuth = if a.is_finite() && a >= 0.0 { a } else { 0.0 };

                // Beam shading: the cell is in cast shadow when the sun sits
                // below the terrain horizon at its azimuth.
                let h = horizon.interpolate(r, c, sun_az_rad);
                let decomp = if sun_elev_rad < h {
                    // No beam; only the diffuse sky and ground reach the cell.
                    Decomposition { ghi: base.dhi, dni: 0.0, dhi: base.dhi }
                } else {
                    base
                };

                let poa = poa_irradiance(
                    cfg.sky_model,
                    tilt,
                    surface_azimuth,
                    decomp,
                    sun.apparent_zenith,
                    sun.azimuth,
                    cfg.albedo,
                    dni_extra,
                    airmass,
                );

                let e_poa = poa.global * dt_hours;
                let e_ac = ac_power(&cfg.system, poa.global, cfg.temp_air, cfg.wind) * dt_hours;
                // Accumulate (errors only on out-of-bounds, impossible here).
                let _ = poa_wh.set(r, c, poa_wh.get(r, c).unwrap_or(0.0) + e_poa);
                let _ = ac_wh.set(r, c, ac_wh.get(r, c).unwrap_or(0.0) + e_ac);
            }
        }
    }

    // Specific yield = AC energy (Wh) / nameplate (W) = kWh/kWp.
    let mut specific_yield = dem.like(0.0);
    for r in 0..rows {
        for c in 0..cols {
            let sy = ac_wh.get(r, c).unwrap_or(0.0) / cfg.system.pdc0;
            let _ = specific_yield.set(r, c, sy);
        }
    }

    Ok(GridResult { poa_wh, ac_wh, specific_yield, sun_steps })
}
