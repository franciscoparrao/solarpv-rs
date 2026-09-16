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

use crate::crs::{inverse_utm, utm_zone_from_epsg};
use crate::error::{Error, Result};
use crate::irradiance::{
    erbs, extra_radiation, haurwitz_clearsky_ghi, poa_irradiance, relative_airmass, Decomposition,
    SkyModel,
};
use crate::losses::{IamModel, SpectralLoss};
use crate::pv::{ac_power, PvSystem};
use crate::solpos::{
    angle_of_incidence, solar_ephemeris, solar_position, solar_position_at, DateTimeUtc, Location,
    SolarPosition,
};
use crate::spa::{solar_position_spa, SpaParams};
use crate::tracking::{dual_axis, single_axis, DualAxisTracker, SingleAxisTracker};

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
    /// Dual-axis tracker; the module normal tracks the sun at every time step.
    /// Terrain slope/aspect are ignored for orientation but still shade beam.
    DualAxis(DualAxisTracker),
}

/// How the latitude/longitude used for solar geometry is chosen across the grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LatitudeMode {
    /// One scene-centre location for the whole grid (the small-DEM
    /// approximation). Fastest; correct for scenes spanning a fraction of a
    /// degree.
    Center,
    /// Per-cell latitude/longitude derived from the DEM's georeferencing. Use
    /// for scenes large enough that latitude varies meaningfully across the
    /// grid.
    ///
    /// The DEM may be geographic (lon/lat degrees, e.g. EPSG:4326) or a WGS84
    /// UTM grid (e.g. EPSG:32719), in which case easting/northing are inverted
    /// to latitude/longitude analytically. A projected DEM whose CRS is unknown
    /// or non-UTM is rejected rather than misread as degrees — see
    /// [`crate::crs`].
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
    /// Local atmospheric pressure (Pa), used to convert relative airmass to
    /// absolute airmass for the spectral-loss model. Default: sea level.
    pub pressure_pa: f64,
    /// Horizon-angle parameters (search radius in cells, number of directions).
    pub horizon: HorizonParams,
    /// How latitude/longitude is chosen per cell for solar geometry.
    pub latitude_mode: LatitudeMode,
    /// How the modules are mounted (fixed to terrain, fixed tilt, or tracking).
    pub mount: Mount,
    /// Use the NREL SPA for scene-centre solar position when `Some` (highest
    /// accuracy); `None` uses the faster Michalsky model. Ignored by
    /// [`LatitudeMode::PerCellGeographic`], which relies on shared ephemerides.
    pub spa: Option<SpaParams>,
    /// Incidence-angle modifier applied to the beam component when `Some`
    /// (angular reflection loss). `None` leaves the beam unmodified.
    pub iam: Option<IamModel>,
    /// SAPM spectral mismatch modifier applied to the effective irradiance when
    /// `Some`. `None` leaves the irradiance unmodified.
    pub spectral: Option<SpectralLoss>,
    /// Apply a topographic sky-view factor to the diffuse-sky component.
    ///
    /// The SVF is computed from the terrain horizon angles and reduces the
    /// sky-diffuse irradiance reaching each cell (e.g. in valleys or near
    /// ridges). Default `false` preserves the unobstructed-sky assumption used
    /// for pvlib parity.
    pub apply_sky_view_factor: bool,
    /// Process the DEM in square tiles of this many interior cells per side,
    /// bounding peak memory (the horizon array is `8 × directions × rows × cols`
    /// bytes — the dominant cost on large scenes). `None` computes the whole DEM
    /// at once. Each tile is grown by a halo of `horizon.radius` cells so that
    /// interior cells still see the surrounding topography; results are
    /// identical to the untiled run.
    pub tile: Option<usize>,
    /// Rescale the clear-sky model to an observed mean-daily GHI raster (e.g.
    /// the Explorador Solar), capturing clouds and the coastal *camanchaca* that
    /// the clear-sky model cannot. `None` keeps the pure clear-sky potential.
    ///
    /// Per cell and day the horizontal clear-sky GHI is scaled by a clearness
    /// index `k = observed / clear-sky` before Erbs decomposition and
    /// transposition, so the diurnal shape and terrain shading are preserved
    /// while the integrated horizontal energy matches the observation.
    pub observed_ghi: Option<ObservedGhi>,
}

/// Observed mean-daily GHI (Wh/m²/day), aligned cell-for-cell to the DEM, used
/// to rescale the clear-sky model (see [`GridConfig::observed_ghi`]).
#[derive(Debug, Clone)]
pub enum ObservedGhi {
    /// One raster of annual mean-daily GHI; the same clearness index scales
    /// every day.
    Annual(Raster<f64>),
    /// Twelve rasters of monthly mean-daily GHI (index 0 = January); each
    /// sampled day is scaled by its own month, capturing the seasonal cloud
    /// cycle. Preferred where seasonal variation matters.
    Monthly(Box<[Raster<f64>; 12]>),
}

impl ObservedGhi {
    /// Observed mean-daily GHI (Wh/m²/day) for a cell in a given month, or
    /// `None` where the raster has no data there.
    fn daily(&self, month: u32, r: usize, c: usize) -> Option<f64> {
        let raster = match self {
            ObservedGhi::Annual(x) => x,
            ObservedGhi::Monthly(m) => &m[((month.max(1) - 1) % 12) as usize],
        };
        let v = raster.get(r, c).ok()?;
        (v.is_finite() && v > 0.0).then_some(v)
    }
}

impl GridConfig {
    /// Sensible defaults for a small DEM: 15-minute steps, Perez sky, desert
    /// albedo, the 1 kW reference system, 18 °C / 2 m·s⁻¹, sea-level pressure.
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
            pressure_pa: 101_325.0,
            horizon: HorizonParams::default(),
            latitude_mode: LatitudeMode::Center,
            mount: Mount::FixedTerrain,
            spa: None,
            iam: None,
            spectral: None,
            apply_sky_view_factor: false,
            tile: None,
            observed_ghi: None,
        }
    }
}

/// Scene-centre solar position using the configured algorithm (SPA or Michalsky).
fn scene_sun(cfg: &GridConfig, when: DateTimeUtc) -> SolarPosition {
    match cfg.spa {
        Some(p) => solar_position_spa(when, cfg.center, &p),
        None => solar_position(when, cfg.center),
    }
}

/// Precomputed solar state plus the irradiance/weather for one time step
/// (scene geometry and the decomposition shared by every cell).
#[derive(Clone, Copy)]
struct SunStep {
    zenith: f64,
    azimuth: f64,
    elev_rad: f64,
    az_rad: f64,
    airmass: f64,
    dni_extra: f64,
    /// Horizontal clear-sky GHI (W/m²) driving `base`, kept so the step can be
    /// rescaled by an observed clearness index without recomputing geometry.
    clearsky_ghi: f64,
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
    svf: Option<&Raster<f64>>,
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
        Mount::DualAxis(tracker) => match dual_axis(&tracker, s.zenith, s.azimuth) {
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
    // Topographic sky-view factor: reduces the visible fraction of the sky dome
    // and therefore the diffuse-sky component reaching the plane of array.
    let svf = svf.map_or(1.0, |raster| raster.get(r, c).unwrap_or(1.0).clamp(0.0, 1.0));
    let sky_diffuse = poa.sky_diffuse * svf;
    let global = poa.direct + sky_diffuse + poa.ground_diffuse;
    // Reported POA is geometric; IAM reduces only the effective irradiance the
    // cells convert (beam reflection loss), keeping the two quantities distinct.
    let effective = match cfg.iam {
        Some(model) => {
            let aoi = angle_of_incidence(tilt, surface_azimuth, s.zenith, s.azimuth);
            poa.direct * model.iam(aoi) + sky_diffuse + poa.ground_diffuse
        }
        None => global,
    };
    // SAPM spectral mismatch: pressure-corrected airmass times relative airmass.
    let effective = cfg.spectral.map_or(effective, |sp| {
        let am_abs = s.airmass * (cfg.pressure_pa / 101_325.0);
        effective * sp.factor(am_abs)
    });
    let e_poa = global * s.dt_hours;
    let e_ac = ac_power(&cfg.system, effective, s.temp_air, s.wind) * s.dt_hours;
    (e_poa, e_ac)
}

/// [`step_energy`] with an optional observed clearness index `k`: when `Some`,
/// the clear-sky GHI is scaled by `k` and re-decomposed (Erbs) before
/// transposition, so an observed GHI raster drives the energy while the clear-
/// sky diurnal shape and terrain shading are preserved.
#[allow(clippy::too_many_arguments)]
fn step_energy_scaled(
    cfg: &GridConfig,
    horizon: &HorizonAngles,
    svf: Option<&Raster<f64>>,
    r: usize,
    c: usize,
    terrain_tilt: f64,
    terrain_azimuth: f64,
    s: &SunStep,
    doy: u32,
    k: Option<f64>,
) -> (f64, f64) {
    match k {
        None => step_energy(cfg, horizon, svf, r, c, terrain_tilt, terrain_azimuth, s),
        Some(k) => {
            let mut scaled = *s;
            scaled.base = erbs(k * s.clearsky_ghi, s.zenith, doy);
            step_energy(cfg, horizon, svf, r, c, terrain_tilt, terrain_azimuth, &scaled)
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

/// Precomputed terrain layers (slope, aspect, skyline), independent of the day —
/// computed once and reused across every day of an annual run.
struct Terrain {
    slope_rad: Raster<f64>,
    asp_deg: Raster<f64>,
    horizon: HorizonAngles,
    /// Optional topographic sky-view factor raster (0 = fully obstructed,
    /// 1 = full hemisphere visible).
    svf: Option<Raster<f64>>,
}

/// Sky-view factor from horizon angles: `SVF = 1 − mean(sin²(horizon))`.
///
/// A flat horizon yields `1.0`; a fully obstructed sky yields `0.0`. The raster
/// copies the DEM's georeferencing.
fn compute_sky_view_factor(horizon: &HorizonAngles, dem: &Raster<f64>) -> Raster<f64> {
    let (rows, cols) = horizon.shape();
    let n = horizon.directions();
    let mut data = vec![0.0; rows * cols];
    for r in 0..rows {
        for c in 0..cols {
            let mut sum_sin2 = 0.0;
            let mut count = 0usize;
            for d in 0..n {
                let h = horizon.get(d, r, c);
                if h.is_finite() {
                    sum_sin2 += h.sin().powi(2);
                    count += 1;
                }
            }
            data[r * cols + c] = if count > 0 { 1.0 - sum_sin2 / count as f64 } else { 1.0 };
        }
    }
    let mut raster = Raster::from_vec(data, rows, cols).expect("dimensions match horizon");
    raster.set_transform(*dem.transform());
    raster
}

fn build_terrain(dem: &Raster<f64>, cfg: &GridConfig) -> Result<Terrain> {
    let mut slope_params = SlopeParams::default();
    slope_params.units = SlopeUnits::Radians;
    let slope_rad = slope(dem, slope_params)
        .map_err(|e| Error::Terrain(format!("slope: {e}")))?;
    let asp_deg = aspect(dem, AspectOutput::Degrees)
        .map_err(|e| Error::Terrain(format!("aspect: {e}")))?;
    let horizon = horizon_angles(dem, cfg.horizon.clone())
        .map_err(|e| Error::Terrain(format!("horizon_angles: {e}")))?;
    let svf = if cfg.apply_sky_view_factor {
        Some(compute_sky_view_factor(&horizon, dem))
    } else {
        None
    };
    Ok(Terrain { slope_rad, asp_deg, horizon, svf })
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

/// How a DEM cell's `(col, row)` becomes a `(latitude, longitude)` in degrees
/// for [`LatitudeMode::PerCellGeographic`].
///
/// Resolved once per run from the DEM's CRS (falling back to the coordinate
/// magnitudes when no CRS is tagged), so the hot per-cell loop only does the
/// arithmetic.
#[derive(Debug, Clone, Copy)]
enum LatLonSource {
    /// The transform is already geographic: `(x, y) = (lon, lat)` in degrees.
    GeographicDegrees,
    /// The transform is a WGS84 UTM grid; invert easting/northing.
    Utm { zone: u8, north: bool },
}

impl LatLonSource {
    #[inline]
    fn lat_lon(self, transform: &surtgis_core::GeoTransform, col: usize, row: usize) -> Location {
        let (x, y) = transform.pixel_to_geo(col, row);
        match self {
            LatLonSource::GeographicDegrees => Location { latitude: y, longitude: x },
            LatLonSource::Utm { zone, north } => {
                let (lat, lon) = inverse_utm(x, y, zone, north);
                Location { latitude: lat, longitude: lon }
            }
        }
    }
}

/// Decide how to turn DEM coordinates into latitude/longitude for per-cell solar
/// geometry, or fail with a clear message when the DEM cannot support it.
///
/// A projected UTM DEM (the usual case for terrain work) is inverted to
/// lat/lon; a geographic DEM is used directly. A projected DEM whose CRS is
/// missing or non-UTM is **rejected** rather than silently read as degrees —
/// which is the bug that turned UTM northings into millions-of-degrees
/// "latitudes" and produced a per-row-striped yield map.
fn resolve_latlon_source(dem: &Raster<f64>) -> Result<LatLonSource> {
    if let Some(crs) = dem.crs() {
        if crs.is_geographic() {
            return Ok(LatLonSource::GeographicDegrees);
        }
        if let Some(epsg) = crs.epsg() {
            if let Some((zone, north)) = utm_zone_from_epsg(epsg) {
                return Ok(LatLonSource::Utm { zone, north });
            }
            return Err(Error::Terrain(format!(
                "per-cell latitude needs a geographic or WGS84 UTM DEM; EPSG:{epsg} is projected \
                 but not UTM. Reproject the DEM to its UTM zone, or use scene-centre latitude \
                 (LatitudeMode::Center)."
            )));
        }
        // CRS present but neither geographic nor EPSG-coded: fall through to the
        // coordinate-magnitude check below.
    }

    // No usable CRS. Infer from the corner coordinates: real lon/lat stay within
    // ±360 / ±90, so anything larger is projected metres we cannot interpret
    // without knowing the projection — refuse instead of guessing.
    let t = dem.transform();
    let (rows, cols) = dem.shape();
    let looks_geographic = [(0, 0), (cols - 1, 0), (0, rows - 1), (cols - 1, rows - 1)]
        .iter()
        .all(|&(c, r)| {
            let (x, y) = t.pixel_to_geo(c, r);
            x.abs() <= 360.0 && y.abs() <= 90.0
        });
    if looks_geographic {
        Ok(LatLonSource::GeographicDegrees)
    } else {
        Err(Error::Terrain(
            "per-cell latitude requested but the DEM has no CRS and its coordinates are not \
             longitude/latitude (they look like projected metres, e.g. UTM). Tag the DEM with its \
             CRS (e.g. EPSG:32719) so latitude can be derived, or use scene-centre latitude \
             (LatitudeMode::Center)."
                .into(),
        ))
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
            let svf = terrain.svf.as_ref();

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
        let sun = scene_sun(cfg, when);
        if sun.apparent_elevation > 0.0 {
            let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
            sun_steps_vec.push(SunStep {
                zenith: sun.apparent_zenith,
                azimuth: sun.azimuth,
                elev_rad: sun.apparent_elevation * DEG,
                az_rad: sun.azimuth * DEG,
                airmass: relative_airmass(sun.apparent_zenith),
                dni_extra,
                clearsky_ghi: ghi,
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

    // For per-cell latitude, decide once how to map DEM coordinates to lat/lon
    // (and fail fast on an un-interpretable DEM). Center mode needs no CRS.
    let latlon_source = match cfg.latitude_mode {
        LatitudeMode::PerCellGeographic => Some(resolve_latlon_source(dem)?),
        LatitudeMode::Center => None,
    };

    // Observed-GHI rescaling context: the month this day belongs to and, for
    // Center mode, the scene-centre clear-sky daily GHI (Wh/m²/day) that the
    // observed value is divided by to get each cell's clearness index.
    let observed = cfg.observed_ghi.as_ref();
    let month = date.month;
    let cs_daily_centre: f64 = sun_steps_vec.iter().map(|s| s.clearsky_ghi * s.dt_hours).sum();

    // Per-cell accumulation, parallelised over the flattened grid. Each cell is
    // independent; slope/aspect/horizon are read-only and Sync.
    let energies: Vec<(f64, f64)> = (0..rows * cols)
        .into_par_iter()
        .map(|idx| {
            let (r, c) = (idx / cols, idx % cols);

            // Observed GHI (if configured): a cell with no coverage yields nothing.
            let observed_daily = match observed {
                Some(o) => match o.daily(month, r, c) {
                    Some(v) => Some(v),
                    None => return (0.0, 0.0),
                },
                None => None,
            };

            // Border cells (Horn's method undefined) come back as NaN → flat.
            let raw_tilt = slope_rad.get(r, c).unwrap_or(0.0);
            let tilt = if raw_tilt.is_finite() { raw_tilt.to_degrees() } else { 0.0 };
            let a = asp_deg.get(r, c).unwrap_or(-1.0);
            let surface_azimuth = if a.is_finite() && a >= 0.0 { a } else { 0.0 };

            // Clearness index k = observed / clear-sky daily horizontal GHI.
            let clearness = |cs_daily: f64| -> Option<f64> {
                observed_daily.map(|obs| if cs_daily > 0.0 { obs / cs_daily } else { 0.0 })
            };

            let mut poa_acc = 0.0;
            let mut ac_acc = 0.0;
            match cfg.latitude_mode {
                LatitudeMode::Center => {
                    let k = clearness(cs_daily_centre);
                    for s in &sun_steps_vec {
                        let (ep, ea) = step_energy_scaled(
                            cfg, horizon, svf, r, c, tilt, surface_azimuth, s, doy, k,
                        );
                        poa_acc += ep;
                        ac_acc += ea;
                    }
                }
                LatitudeMode::PerCellGeographic => {
                    // Latitude/longitude derived from the DEM georeferencing
                    // (geographic direct, or UTM inverted); resolved once above.
                    let loc = latlon_source
                        .expect("per-cell source resolved for PerCellGeographic")
                        .lat_lon(&transform, c, r);
                    let build_step = |sun: &SolarPosition| -> SunStep {
                        let ghi = haurwitz_clearsky_ghi(sun.apparent_zenith);
                        SunStep {
                            zenith: sun.apparent_zenith,
                            azimuth: sun.azimuth,
                            elev_rad: sun.apparent_elevation * DEG,
                            az_rad: sun.azimuth * DEG,
                            airmass: relative_airmass(sun.apparent_zenith),
                            dni_extra,
                            clearsky_ghi: ghi,
                            base: erbs(ghi, sun.apparent_zenith, doy),
                            temp_air: cfg.temp_air,
                            wind: cfg.wind,
                            dt_hours,
                        }
                    };
                    if observed_daily.is_none() {
                        // Fast path: no rescaling, stream steps without buffering.
                        for eph in &ephemerides {
                            let sun = solar_position_at(eph, loc);
                            if sun.apparent_elevation <= 0.0 {
                                continue;
                            }
                            let s = build_step(&sun);
                            let (ep, ea) =
                                step_energy(cfg, horizon, svf, r, c, tilt, surface_azimuth, &s);
                            poa_acc += ep;
                            ac_acc += ea;
                        }
                    } else {
                        // Rescaling needs this cell's clear-sky daily GHI, so
                        // buffer the day's steps, then apply k.
                        let mut steps: Vec<SunStep> = Vec::with_capacity(ephemerides.len());
                        for eph in &ephemerides {
                            let sun = solar_position_at(eph, loc);
                            if sun.apparent_elevation <= 0.0 {
                                continue;
                            }
                            steps.push(build_step(&sun));
                        }
                        let cs_daily: f64 =
                            steps.iter().map(|s| s.clearsky_ghi * s.dt_hours).sum();
                        let k = clearness(cs_daily);
                        for s in &steps {
                            let (ep, ea) = step_energy_scaled(
                                cfg, horizon, svf, r, c, tilt, surface_azimuth, s, doy, k,
                            );
                            poa_acc += ep;
                            ac_acc += ea;
                        }
                    }
                }
            }
            (poa_acc, ac_acc)
        })
        .collect();

    Ok((energies, sun_steps))
}

/// Copy a `[r0..r0+h) × [c0..c0+w)` window of the DEM into a standalone raster,
/// carrying a correctly shifted transform plus the CRS and nodata, so terrain
/// and per-cell latitude resolve exactly as on the parent.
fn extract_subdem(dem: &Raster<f64>, r0: usize, c0: usize, h: usize, w: usize) -> Raster<f64> {
    let src = dem.data();
    let mut data = Vec::with_capacity(h * w);
    for rr in 0..h {
        for cc in 0..w {
            data.push(src[[r0 + rr, c0 + cc]]);
        }
    }
    let mut sub = Raster::from_vec(data, h, w).expect("sub-DEM dimensions are valid");
    // Shift the origin to the window's top-left corner; keep pixel size and any
    // rotation so georeferencing (and thus per-cell UTM inversion) is exact.
    let mut t = *dem.transform();
    let (ox, oy) = t.pixel_to_geo_corner(c0, r0);
    t.origin_x = ox;
    t.origin_y = oy;
    sub.set_transform(t);
    if let Some(crs) = dem.crs() {
        sub.set_crs(Some(crs.clone()));
    }
    if let Some(nd) = dem.nodata() {
        sub.set_nodata(Some(nd));
    }
    sub
}

/// Drive a whole-DEM computation tile by tile, bounding peak memory.
///
/// The DEM is walked in `tile × tile` interior blocks, each grown by a
/// `horizon.radius` halo (clamped to the DEM) so that interior cells see the
/// surrounding relief. `per_tile` computes the per-cell `(POA, AC)` energy for a
/// sub-DEM and its terrain; only interior cells are scattered into the full
/// output, so the result is identical to the untiled computation.
fn tiled_energies<F>(
    dem: &Raster<f64>,
    cfg: &GridConfig,
    tile: usize,
    per_tile: F,
) -> Result<(Vec<f64>, Vec<f64>, usize)>
where
    F: Fn(&Raster<f64>, &Terrain) -> Result<(Vec<(f64, f64)>, usize)>,
{
    let (rows, cols) = dem.shape();
    // Fail fast on an un-interpretable DEM before doing any tile work.
    if cfg.latitude_mode == LatitudeMode::PerCellGeographic {
        resolve_latlon_source(dem)?;
    }

    let halo = cfg.horizon.radius;
    let tile = tile.max(1);
    let mut poa = vec![0.0f64; rows * cols];
    let mut ac = vec![0.0f64; rows * cols];
    let mut sun_steps = 0usize;

    let mut r0 = 0;
    while r0 < rows {
        let ih = tile.min(rows - r0);
        let er0 = r0.saturating_sub(halo);
        let er1 = (r0 + ih + halo).min(rows);
        let mut c0 = 0;
        while c0 < cols {
            let iw = tile.min(cols - c0);
            let ec0 = c0.saturating_sub(halo);
            let ec1 = (c0 + iw + halo).min(cols);

            let sub = extract_subdem(dem, er0, ec0, er1 - er0, ec1 - ec0);
            let terrain = build_terrain(&sub, cfg)?;
            let (energies, steps) = per_tile(&sub, &terrain)?;
            sun_steps = steps;

            // Scatter interior cells (offset by the halo) into the full grid.
            let sub_cols = ec1 - ec0;
            let ir = r0 - er0;
            let ic = c0 - ec0;
            for rr in 0..ih {
                for cc in 0..iw {
                    let sidx = (ir + rr) * sub_cols + (ic + cc);
                    let gidx = (r0 + rr) * cols + (c0 + cc);
                    poa[gidx] = energies[sidx].0;
                    ac[gidx] = energies[sidx].1;
                }
            }
            c0 += tile;
        }
        r0 += tile;
    }
    Ok((poa, ac, sun_steps))
}

/// Compute the daily clear-sky PV potential over a DEM.
///
/// Returns per-cell POA insolation (Wh/m²/day), AC energy (Wh/day) and specific
/// yield (kWh/kWp/day). Reuses SurtGIS `slope`/`aspect`/`horizon_angles`; the
/// per-cell physics is the same validated point chain used elsewhere. Set
/// [`GridConfig::tile`] to bound memory on large scenes.
pub fn pv_potential(dem: &Raster<f64>, cfg: &GridConfig) -> Result<GridResult> {
    let (poa_vec, ac_vec, sun_steps) = match cfg.tile {
        Some(t) => tiled_energies(dem, cfg, t, |sub, terrain| {
            day_energies(sub, cfg, terrain, cfg.date)
        })?,
        None => {
            let terrain = build_terrain(dem, cfg)?;
            let (energies, sun_steps) = day_energies(dem, cfg, &terrain, cfg.date)?;
            (
                energies.iter().map(|e| e.0).collect(),
                energies.iter().map(|e| e.1).collect(),
                sun_steps,
            )
        }
    };
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
    let days = sampled_days(cfg.date.year, sampling)?;

    // Sum the sampled days (weighted) for one sub-DEM + terrain.
    let annual_tile = |sub: &Raster<f64>, terrain: &Terrain| -> Result<(Vec<(f64, f64)>, usize)> {
        let (sr, sc) = sub.shape();
        let mut acc = vec![(0.0f64, 0.0f64); sr * sc];
        for (date, weight) in &days {
            let (energies, _) = day_energies(sub, cfg, terrain, *date)?;
            for (i, (p, a)) in energies.iter().enumerate() {
                acc[i].0 += p * weight;
                acc[i].1 += a * weight;
            }
        }
        Ok((acc, days.len()))
    };

    let (poa_tot, ac_tot, _) = match cfg.tile {
        Some(t) => tiled_energies(dem, cfg, t, annual_tile)?,
        None => {
            let (rows, cols) = dem.shape();
            let terrain = build_terrain(dem, cfg)?;
            let (energies, _) = annual_tile(dem, &terrain)?;
            let mut poa = vec![0.0f64; rows * cols];
            let mut ac = vec![0.0f64; rows * cols];
            for (i, (p, a)) in energies.iter().enumerate() {
                poa[i] = *p;
                ac[i] = *a;
            }
            (poa, ac, days.len())
        }
    };
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
    let terrain = build_terrain(dem, cfg)?;
    let slope_rad = &terrain.slope_rad;
    let asp_deg = &terrain.asp_deg;
    let horizon = &terrain.horizon;
    let svf = terrain.svf.as_ref();

    // Build the per-step solar + weather state at the scene centre.
    let mut steps: Vec<SunStep> = Vec::with_capacity(records.len());
    for rec in records {
        let sun = scene_sun(cfg, rec.when);
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
            // Measured GHI already drives this path; observed-GHI rescaling does
            // not apply, but keep the field consistent with the decomposition.
            clearsky_ghi: rec.ghi,
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
                let (ep, ea) = step_energy(cfg, horizon, svf, r, c, tilt, surface_azimuth, s);
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
