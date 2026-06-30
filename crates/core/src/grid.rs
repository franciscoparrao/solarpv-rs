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
//! Erbs ([`pv_potential`], [`pv_potential_annual`]). To drive the engine from a
//! measured or TMY irradiance series instead, use [`pv_potential_series`].
//!
//! Scene geometry is taken at the scene centre by default; [`LatitudeMode`]
//! switches to per-cell latitude/longitude for scenes large enough that it
//! matters. Modules can be ground-following, fixed-tilt or tracking ([`Mount`]).

use rayon::prelude::*;
use surtgis_algorithms::terrain::{
    aspect, horizon_angles, slope, AspectOutput, HorizonAngles, HorizonParams, SlopeParams,
    SlopeUnits,
};
use surtgis_core::Raster;

use crate::error::{Error, Result};
use crate::irradiance::{
    erbs, extra_radiation, haurwitz_clearsky_ghi, poa_irradiance, relative_airmass, Decomposition,
    SkyModel,
};
use crate::pv::{ac_power, PvSystem};
use crate::solpos::{solar_ephemeris, solar_position, solar_position_at, DateTimeUtc, Location};
use crate::tracking::{single_axis, SingleAxisTracker};

const DEG: f64 = std::f64::consts::PI / 180.0;

/// How the modules are mounted on each cell, which sets their orientation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mount {
    /// Ground-following: panel tilt/azimuth equal the terrain slope/aspect.
    FixedTerrain,
    /// Fixed racks at a chosen tilt and azimuth (degrees, azimuth clockwise
    /// from North), the same on every cell.
    FixedTilt { tilt: f64, surface_azimuth: f64 },
    /// Single-axis tracker; the orientation follows the sun each time step and
    /// the terrain only contributes horizon shading.
    SingleAxis(SingleAxisTracker),
}

/// How the latitude/longitude used for solar geometry is chosen across the grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatitudeMode {
    /// One scene-centre location for the whole grid (the small-DEM
    /// approximation). Fastest; correct for scenes spanning a fraction of a
    /// degree.
    Center,
    /// Per-cell latitude/longitude read from the DEM's geographic transform
    /// (requires a lon/lat DEM, e.g. EPSG:4326). Use for scenes large enough
    /// that latitude varies meaningfully across the grid.
    PerCellGeographic,
}

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
    /// How latitude/longitude is chosen per cell for solar geometry.
    pub latitude_mode: LatitudeMode,
    /// How the modules are mounted (fixed to terrain, fixed tilt, or tracking).
    pub mount: Mount,
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
            latitude_mode: LatitudeMode::Center,
            mount: Mount::FixedTerrain,
        }
    }
}

/// Precomputed solar state plus the irradiance/weather for one time step
/// (scene geometry and the decomposition shared by every cell).
struct SunStep {
    zenith: f64,
    azimuth: f64,
    elev_rad: f64,
    az_rad: f64,
    airmass: f64,
    dni_extra: f64,
    base: Decomposition,
    temp_air: f64,
    wind: f64,
    /// Hours this step represents in the energy integration.
    dt_hours: f64,
}

/// Shading + transposition + PV for one cell at one [`SunStep`]. Returns the
/// `(POA Wh, AC Wh)` energy contributed by this step.
fn step_energy(
    cfg: &GridConfig,
    horizon: &HorizonAngles,
    r: usize,
    c: usize,
    terrain_tilt: f64,
    terrain_azimuth: f64,
    s: &SunStep,
) -> (f64, f64) {
    let (tilt, surface_azimuth) = match cfg.mount {
        Mount::FixedTerrain => (terrain_tilt, terrain_azimuth),
        Mount::FixedTilt { tilt, surface_azimuth } => (tilt, surface_azimuth),
        Mount::SingleAxis(tracker) => match single_axis(&tracker, s.zenith, s.azimuth) {
            Some(o) => (o.surface_tilt, o.surface_azimuth),
            None => return (0.0, 0.0),
        },
    };
    // Beam shading: cell is in cast shadow when the sun sits below the terrain
    // horizon at its azimuth.
    let h = horizon.interpolate(r, c, s.az_rad);
    let decomp = if s.elev_rad < h {
        Decomposition { ghi: s.base.dhi, dni: 0.0, dhi: s.base.dhi }
    } else {
        s.base
    };
    let poa = poa_irradiance(
        cfg.sky_model, tilt, surface_azimuth, decomp, s.zenith, s.azimuth, cfg.albedo,
        s.dni_extra, s.airmass,
    );
    let e_poa = poa.global * s.dt_hours;
    let e_ac = ac_power(&cfg.system, poa.global, s.temp_air, s.wind) * s.dt_hours;
    (e_poa, e_ac)
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

/// Precomputed terrain layers (slope, aspect, skyline), independent of the day —
/// computed once and reused across every day of an annual run.
struct Terrain {
    slope_rad: Raster<f64>,
    asp_deg: Raster<f64>,
    horizon: HorizonAngles,
}

fn build_terrain(dem: &Raster<f64>, horizon_params: &HorizonParams) -> Result<Terrain> {
    let slope_rad = slope(dem, SlopeParams { units: SlopeUnits::Radians, ..Default::default() })
        .map_err(|e| Error::Terrain(format!("slope: {e}")))?;
    let asp_deg = aspect(dem, AspectOutput::Degrees)
        .map_err(|e| Error::Terrain(format!("aspect: {e}")))?;
    let horizon = horizon_angles(dem, horizon_params.clone())
        .map_err(|e| Error::Terrain(format!("horizon_angles: {e}")))?;
    Ok(Terrain { slope_rad, asp_deg, horizon })
}

/// Build the output rasters (POA, AC and specific yield) from flattened energy
/// vectors, copying the DEM's georeferencing.
fn build_result(
    dem: &Raster<f64>,
    cfg: &GridConfig,
    poa_vec: Vec<f64>,
    ac_vec: Vec<f64>,
    sun_steps: usize,
) -> GridResult {
    let (rows, cols) = dem.shape();
    let sy_vec: Vec<f64> = ac_vec.iter().map(|&ac| ac / cfg.system.pdc0).collect();
    let to_raster = |v: Vec<f64>| -> Raster<f64> {
        let mut out = Raster::from_vec(v, rows, cols).expect("dimensions match dem");
        out.set_transform(*dem.transform());
        out
    };
    GridResult {
        poa_wh: to_raster(poa_vec),
        ac_wh: to_raster(ac_vec),
        specific_yield: to_raster(sy_vec),
        sun_steps,
    }
}

/// Per-cell `(POA Wh, AC Wh)` energy for a single `date`, plus the number of
/// daylight steps at the scene centre. Terrain is supplied precomputed.
fn day_energies(
    dem: &Raster<f64>,
    cfg: &GridConfig,
    terrain: &Terrain,
    date: DateTimeUtc,
) -> Result<(Vec<(f64, f64)>, usize)> {
    let (rows, cols) = dem.shape();
    let slope_rad = &terrain.slope_rad;
    let asp_deg = &terrain.asp_deg;
    let horizon = &terrain.horizon;

    let dt_hours = cfg.time_step_minutes as f64 / 60.0;
    let steps = (24 * 60) / cfg.time_step_minutes.max(1);
    let doy = date.day_of_year();
    let dni_extra = extra_radiation(doy);

    // Precompute, for every step, the time-only ephemeris (shared by all cells)
    // and the scene-centre solar state for the fast Center path.
    let mut ephemerides: Vec<crate::solpos::SolarEphemeris> = Vec::with_capacity(steps as usize);
    let mut sun_steps_vec: Vec<SunStep> = Vec::new();
    for i in 0..steps {
        let total_min = i * cfg.time_step_minutes;
        let when = DateTimeUtc::new(
            date.year,
            date.month,
            date.day,
            total_min / 60,
            total_min % 60,
            0,
        )?;
        let eph = solar_ephemeris(when);
        let sun = solar_position(when, cfg.center);
        if sun.apparent_elevation > 0.0 {
            let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
            sun_steps_vec.push(SunStep {
                zenith: sun.apparent_zenith,
                azimuth: sun.azimuth,
                elev_rad: sun.apparent_elevation * DEG,
                az_rad: sun.azimuth * DEG,
                airmass: relative_airmass(sun.apparent_zenith),
                dni_extra,
                base: erbs(ghi, sun.apparent_zenith, doy),
                temp_air: cfg.temp_air,
                wind: cfg.wind,
                dt_hours,
            });
        }
        ephemerides.push(eph);
    }
    let sun_steps = sun_steps_vec.len();

    let transform = *dem.transform();

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
            match cfg.latitude_mode {
                LatitudeMode::Center => {
                    for s in &sun_steps_vec {
                        let (ep, ea) = step_energy(cfg, horizon, r, c, tilt, surface_azimuth, s);
                        poa_acc += ep;
                        ac_acc += ea;
                    }
                }
                LatitudeMode::PerCellGeographic => {
                    // (x, y) = (lon, lat) from the geographic transform.
                    let (lon, lat) = transform.pixel_to_geo(c, r);
                    let loc = Location { latitude: lat, longitude: lon };
                    for eph in &ephemerides {
                        let sun = solar_position_at(eph, loc);
                        if sun.apparent_elevation <= 0.0 {
                            continue;
                        }
                        let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
                        let s = SunStep {
                            zenith: sun.apparent_zenith,
                            azimuth: sun.azimuth,
                            elev_rad: sun.apparent_elevation * DEG,
                            az_rad: sun.azimuth * DEG,
                            airmass: relative_airmass(sun.apparent_zenith),
                            dni_extra,
                            base: erbs(ghi, sun.apparent_zenith, doy),
                            temp_air: cfg.temp_air,
                            wind: cfg.wind,
                            dt_hours,
                        };
                        let (ep, ea) = step_energy(cfg, horizon, r, c, tilt, surface_azimuth, &s);
                        poa_acc += ep;
                        ac_acc += ea;
                    }
                }
            }
            (poa_acc, ac_acc)
        })
        .collect();

    Ok((energies, sun_steps))
}

/// Compute the daily clear-sky PV potential over a DEM.
///
/// Returns per-cell POA insolation (Wh/m²/day), AC energy (Wh/day) and specific
/// yield (kWh/kWp/day). Reuses SurtGIS `slope`/`aspect`/`horizon_angles`; the
/// per-cell physics is the same validated point chain used elsewhere.
pub fn pv_potential(dem: &Raster<f64>, cfg: &GridConfig) -> Result<GridResult> {
    let terrain = build_terrain(dem, &cfg.horizon)?;
    let (energies, sun_steps) = day_energies(dem, cfg, &terrain, cfg.date)?;
    let poa_vec = energies.iter().map(|e| e.0).collect();
    let ac_vec = energies.iter().map(|e| e.1).collect();
    Ok(build_result(dem, cfg, poa_vec, ac_vec, sun_steps))
}

/// How the year is sampled when integrating annual PV potential.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DaySampling {
    /// Evaluate every `n`-th day of the year; each sampled day represents `n`
    /// days (the tail is clamped so the weights sum to the year length).
    EveryNDays(u32),
    /// One representative day per month (Duffie & Beckman recommended average
    /// days), each weighted by the number of days in its month. Only 12 daily
    /// evaluations — fast and standard for annual estimates.
    MonthlyRepresentative,
}

/// The (date, weight-in-days) pairs to integrate for a year, per [`DaySampling`].
fn sampled_days(year: i32, sampling: DaySampling) -> Result<Vec<(DateTimeUtc, f64)>> {
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let year_len = if leap { 366 } else { 365 };
    match sampling {
        DaySampling::EveryNDays(n) => {
            let n = n.max(1);
            let mut out = Vec::new();
            let mut doy = 1;
            while doy <= year_len {
                let weight = n.min(year_len - doy + 1) as f64;
                out.push((DateTimeUtc::from_ordinal(year, doy)?, weight));
                doy += n;
            }
            Ok(out)
        }
        DaySampling::MonthlyRepresentative => {
            // Duffie & Beckman average day-of-month, (month, day).
            const REP: [(u32, u32); 12] = [
                (1, 17), (2, 16), (3, 16), (4, 15), (5, 15), (6, 11),
                (7, 17), (8, 16), (9, 15), (10, 15), (11, 14), (12, 10),
            ];
            let mdays = [31u32, if leap { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            REP.iter()
                .enumerate()
                .map(|(i, &(m, d))| Ok((DateTimeUtc::new(year, m, d, 0, 0, 0)?, mdays[i] as f64)))
                .collect()
        }
    }
}

/// Compute the annual clear-sky PV potential over a DEM.
///
/// Terrain (slope/aspect/horizon) is computed once and reused across the sampled
/// days. Output rasters carry **annual** totals: POA insolation (Wh/m²/year), AC
/// energy (Wh/year) and specific yield (kWh/kWp/year). `sun_steps` reports the
/// number of representative days integrated.
pub fn pv_potential_annual(
    dem: &Raster<f64>,
    cfg: &GridConfig,
    sampling: DaySampling,
) -> Result<GridResult> {
    let (rows, cols) = dem.shape();
    let terrain = build_terrain(dem, &cfg.horizon)?;
    let days = sampled_days(cfg.date.year, sampling)?;

    let mut poa_tot = vec![0.0f64; rows * cols];
    let mut ac_tot = vec![0.0f64; rows * cols];
    for (date, weight) in &days {
        let (energies, _) = day_energies(dem, cfg, &terrain, *date)?;
        for (i, (p, a)) in energies.iter().enumerate() {
            poa_tot[i] += p * weight;
            ac_tot[i] += a * weight;
        }
    }
    Ok(build_result(dem, cfg, poa_tot, ac_tot, days.len()))
}

/// One time-stamped weather observation driving a measured / TMY run.
///
/// Only `ghi` is required. If `dni`/`dhi` are `None` they are estimated from
/// `ghi` with the Erbs model; if `temp_air`/`wind` are `None` the
/// [`GridConfig`] scalars are used.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct WeatherRecord {
    /// Observation time (UTC).
    pub when: DateTimeUtc,
    /// Global horizontal irradiance, W/m².
    pub ghi: f64,
    /// Direct normal irradiance, W/m² (optional).
    pub dni: Option<f64>,
    /// Diffuse horizontal irradiance, W/m² (optional).
    pub dhi: Option<f64>,
    /// Air temperature, °C (optional).
    pub temp_air: Option<f64>,
    /// Wind speed, m/s (optional).
    pub wind: Option<f64>,
}

/// Compute PV potential over a DEM from a measured / TMY irradiance series.
///
/// Unlike [`pv_potential`] / [`pv_potential_annual`], which synthesise clear-sky
/// irradiance, this consumes an explicit `records` series for the scene and
/// integrates it (each record represents `dt_hours`). Scene-centre geometry is
/// used ([`LatitudeMode::Center`]); the terrain still modulates each cell via
/// orientation and horizon shading. Output totals span whatever period the
/// series covers (e.g. a full year for an 8760-hour TMY).
pub fn pv_potential_series(
    dem: &Raster<f64>,
    cfg: &GridConfig,
    records: &[WeatherRecord],
    dt_hours: f64,
) -> Result<GridResult> {
    let (rows, cols) = dem.shape();
    let terrain = build_terrain(dem, &cfg.horizon)?;
    let slope_rad = &terrain.slope_rad;
    let asp_deg = &terrain.asp_deg;
    let horizon = &terrain.horizon;

    // Build the per-step solar + weather state at the scene centre.
    let mut steps: Vec<SunStep> = Vec::with_capacity(records.len());
    for rec in records {
        let sun = solar_position(rec.when, cfg.center);
        if sun.apparent_elevation <= 0.0 || rec.ghi <= 0.0 {
            continue; // night or no irradiance contributes nothing
        }
        let z = sun.apparent_zenith;
        let doy = rec.when.day_of_year();
        let base = match (rec.dni, rec.dhi) {
            (Some(dni), Some(dhi)) => Decomposition { ghi: rec.ghi, dni, dhi },
            _ => erbs(rec.ghi, z, doy),
        };
        steps.push(SunStep {
            zenith: z,
            azimuth: sun.azimuth,
            elev_rad: sun.apparent_elevation * DEG,
            az_rad: sun.azimuth * DEG,
            airmass: relative_airmass(z),
            dni_extra: extra_radiation(doy),
            base,
            temp_air: rec.temp_air.unwrap_or(cfg.temp_air),
            wind: rec.wind.unwrap_or(cfg.wind),
            dt_hours,
        });
    }
    let n_steps = steps.len();

    let energies: Vec<(f64, f64)> = (0..rows * cols)
        .into_par_iter()
        .map(|idx| {
            let (r, c) = (idx / cols, idx % cols);
            let raw_tilt = slope_rad.get(r, c).unwrap_or(0.0);
            let tilt = if raw_tilt.is_finite() { raw_tilt.to_degrees() } else { 0.0 };
            let a = asp_deg.get(r, c).unwrap_or(-1.0);
            let surface_azimuth = if a.is_finite() && a >= 0.0 { a } else { 0.0 };

            let mut poa_acc = 0.0;
            let mut ac_acc = 0.0;
            for s in &steps {
                let (ep, ea) = step_energy(cfg, horizon, r, c, tilt, surface_azimuth, s);
                poa_acc += ep;
                ac_acc += ea;
            }
            (poa_acc, ac_acc)
        })
        .collect();

    let poa_vec = energies.iter().map(|e| e.0).collect();
    let ac_vec = energies.iter().map(|e| e.1).collect();
    Ok(build_result(dem, cfg, poa_vec, ac_vec, n_steps))
}
