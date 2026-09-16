//! `solarpv` — command-line interface for the terrain PV-potential engine.
//!
//! Reads a DEM GeoTIFF, computes the daily clear-sky PV potential over it and
//! writes per-cell POA insolation, AC energy and specific-yield rasters.

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};

use solarpv_core::grid::{
    pv_potential, pv_potential_annual, pv_potential_series, DaySampling, GridConfig, LatitudeMode,
    Mount, WeatherRecord,
};
use solarpv_core::irradiance::SkyModel;
use solarpv_core::losses::SpectralLoss;
use solarpv_core::pv::PvSystem;
use solarpv_core::solpos::{DateTimeUtc, Location};
use solarpv_core::tracking::SingleAxisTracker;
use surtgis_core::io::{read_geotiff, write_geotiff};

/// Sky-diffuse transposition model (CLI flag).
#[derive(Copy, Clone, Debug, ValueEnum)]
enum Sky {
    Isotropic,
    Haydavies,
    Perez,
}

impl From<Sky> for SkyModel {
    fn from(s: Sky) -> Self {
        match s {
            Sky::Isotropic => SkyModel::Isotropic,
            Sky::Haydavies => SkyModel::HayDavies,
            Sky::Perez => SkyModel::Perez,
        }
    }
}

/// Module mounting (CLI flag).
#[derive(Copy, Clone, Debug, ValueEnum)]
enum MountKind {
    /// Panels follow the terrain slope/aspect.
    Terrain,
    /// Fixed tilt (see `--tilt` / `--surface-azimuth`).
    Tilt,
    /// Horizontal single-axis tracker with backtracking.
    Tracker,
    /// Dual-axis tracker (follows sun in azimuth and elevation).
    DualAxis,
}

/// Incidence-angle-modifier model (CLI flag).
#[derive(Copy, Clone, Debug, ValueEnum)]
enum IamKind {
    /// No angular reflection loss.
    None,
    /// ASHRAE model (b = 0.05).
    Ashrae,
    /// Martín–Ruiz model (a_r = 0.16).
    MartinRuiz,
    /// Physical Fresnel model (n = 1.526).
    Physical,
}

/// Terrain photovoltaic potential from a DEM ("PVGIS lite").
#[derive(Parser, Debug)]
#[command(name = "solarpv", version, about, allow_negative_numbers = true)]
struct Cli {
    /// Input DEM GeoTIFF (single band, elevation in metres).
    #[arg(long)]
    dem: String,

    /// Scene-centre latitude in degrees (north positive).
    #[arg(long)]
    lat: f64,

    /// Scene-centre longitude in degrees (east positive).
    #[arg(long)]
    lon: f64,

    /// Day to evaluate, `YYYY-MM-DD` (UTC).
    #[arg(long)]
    date: String,

    /// Output path prefix; writes `<prefix>_poa.tif`, `_ac.tif`,
    /// `_specific_yield.tif`.
    #[arg(long)]
    out_prefix: String,

    /// Integration time step in minutes.
    #[arg(long, default_value_t = 15)]
    step: u32,

    /// Sky-diffuse transposition model.
    #[arg(long, value_enum, default_value_t = Sky::Perez)]
    sky: Sky,

    /// Ground albedo.
    #[arg(long, default_value_t = 0.25)]
    albedo: f64,

    /// System nameplate DC power at STC (W).
    #[arg(long, default_value_t = 1000.0)]
    pdc0: f64,

    /// Power temperature coefficient (fraction per °C).
    #[arg(long, default_value_t = -0.004)]
    gamma: f64,

    /// Ambient air temperature (°C).
    #[arg(long, default_value_t = 18.0)]
    temp_air: f64,

    /// Wind speed (m/s).
    #[arg(long, default_value_t = 2.0)]
    wind: f64,

    /// Horizon search radius in cells.
    #[arg(long, default_value_t = 100)]
    horizon_radius: usize,

    /// Number of horizon azimuth directions.
    #[arg(long, default_value_t = 36)]
    horizon_dirs: usize,

    /// Use per-cell latitude/longitude from the DEM georeferencing (a
    /// geographic lon/lat DEM, or a WGS84 UTM DEM whose easting/northing are
    /// inverted analytically, e.g. EPSG:32719). A projected DEM with no CRS is
    /// rejected. Default: scene-centre only.
    #[arg(long, default_value_t = false)]
    per_cell_lat: bool,

    /// Integrate over the whole year (the date's year is used; outputs become
    /// annual totals). Without this, a single day is evaluated.
    #[arg(long, default_value_t = false)]
    annual: bool,

    /// Annual day sampling stride: 0 = one representative day per month;
    /// N > 0 = every N-th day of the year. Only used with `--annual`.
    #[arg(long, default_value_t = 0)]
    day_stride: u32,

    /// Module mounting.
    #[arg(long, value_enum, default_value_t = MountKind::Terrain)]
    mount: MountKind,

    /// Fixed-tilt angle in degrees (only for `--mount tilt`).
    #[arg(long, default_value_t = 0.0)]
    tilt: f64,

    /// Fixed-tilt surface azimuth, degrees clockwise from North
    /// (only for `--mount tilt`; 0 = North-facing).
    #[arg(long, default_value_t = 0.0)]
    surface_azimuth: f64,

    /// Tracker ground coverage ratio (only for `--mount tracker`).
    #[arg(long, default_value_t = 2.0 / 7.0)]
    gcr: f64,

    /// Dual-axis tracker maximum tilt from horizontal, degrees (only for
    /// `--mount dual-axis`). 90 = full hemispherical tracking.
    #[arg(long, default_value_t = 90.0)]
    max_tilt: f64,

    /// Drive the run from a measured / TMY irradiance CSV instead of clear-sky.
    /// Header columns: year,month,day,hour,ghi[,dni,dhi,temp_air,wind] (UTC).
    #[arg(long)]
    weather: Option<String>,

    /// Cadence of the weather series in hours (e.g. 1.0 for hourly TMY).
    #[arg(long, default_value_t = 1.0)]
    weather_dt: f64,

    /// Use the high-accuracy NREL SPA for solar position (default: Michalsky).
    #[arg(long, default_value_t = false)]
    spa: bool,

    /// Incidence-angle-modifier model for angular reflection loss.
    #[arg(long, value_enum, default_value_t = IamKind::None)]
    iam: IamKind,

    /// Apply the SAPM crystalline-silicon spectral mismatch factor to the
    /// effective irradiance (uses pressure-corrected absolute airmass).
    #[arg(long, default_value_t = false)]
    spectral: bool,

    /// Atmospheric pressure in pascals (used for absolute airmass and spectral
    /// loss; default: sea level 101325 Pa).
    #[arg(long, default_value_t = 101325.0)]
    pressure: f64,

    /// Apply the topographic sky-view factor to the diffuse-sky component.
    /// Off by default to preserve the unobstructed-sky assumption used for
    /// pvlib parity.
    #[arg(long, default_value_t = false)]
    svf: bool,

    /// System DC loss fraction (soiling, wiring, mismatch, …). Default 0.14.
    #[arg(long, default_value_t = 0.14)]
    loss: f64,
}

/// Minimal CSV reader for the weather series. Maps columns by header name;
/// requires `year,month,day,hour,ghi`, optional `dni,dhi,temp_air,wind`.
fn read_weather_csv(path: &str) -> Result<Vec<WeatherRecord>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = lines.next().context("empty weather CSV")?;
    let cols: Vec<String> = header.split(',').map(|c| c.trim().to_lowercase()).collect();
    let idx = |name: &str| cols.iter().position(|c| c == name);
    let need = |name: &str| -> Result<usize> {
        idx(name).with_context(|| format!("weather CSV missing required column `{name}`"))
    };
    let (iy, imo, id, ih, ig) =
        (need("year")?, need("month")?, need("day")?, need("hour")?, need("ghi")?);
    let (idni, idhi, itemp, iwind) = (idx("dni"), idx("dhi"), idx("temp_air"), idx("wind"));

    let mut out = Vec::new();
    for (n, line) in lines.enumerate() {
        let v: Vec<&str> = line.split(',').map(|s| s.trim()).collect();
        let get = |i: usize| -> Result<f64> {
            v.get(i)
                .and_then(|s| s.parse::<f64>().ok())
                .with_context(|| format!("row {}: bad number in column {i}", n + 2))
        };
        let opt = |oi: Option<usize>| -> Option<f64> {
            oi.and_then(|i| v.get(i)).and_then(|s| s.parse::<f64>().ok())
        };
        let when = DateTimeUtc::new(
            get(iy)? as i32,
            get(imo)? as u32,
            get(id)? as u32,
            get(ih)? as u32,
            0,
            0,
        )
        .map_err(|e| anyhow::anyhow!(e))?;
        out.push(WeatherRecord {
            when,
            ghi: get(ig)?,
            dni: opt(idni),
            dhi: opt(idhi),
            temp_air: opt(itemp),
            wind: opt(iwind),
        });
    }
    if out.is_empty() {
        bail!("weather CSV {path} has no data rows");
    }
    Ok(out)
}

/// Parse `YYYY-MM-DD` into a midnight-UTC datetime.
fn parse_date(s: &str) -> Result<DateTimeUtc> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 3 {
        bail!("date must be YYYY-MM-DD, got `{s}`");
    }
    let year: i32 = parts[0].parse().context("year")?;
    let month: u32 = parts[1].parse().context("month")?;
    let day: u32 = parts[2].parse().context("day")?;
    DateTimeUtc::new(year, month, day, 0, 0, 0).map_err(|e| anyhow::anyhow!(e))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let date = parse_date(&cli.date)?;
    let center = Location::new(cli.lat, cli.lon).map_err(|e| anyhow::anyhow!(e))?;

    let dem = read_geotiff::<f64, _>(&cli.dem, None)
        .with_context(|| format!("reading DEM {}", cli.dem))?;
    let (rows, cols) = dem.shape();
    eprintln!("DEM {}×{} ({} cells), cell size {} m", rows, cols, rows * cols, dem.cell_size());

    let mut system = PvSystem::reference_1kw();
    system.pdc0 = cli.pdc0;
    system.gamma_pdc = cli.gamma;
    system.pdc0_inv = cli.pdc0 / 1.1;
    system.system_losses = cli.loss;

    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = cli.step;
    cfg.sky_model = cli.sky.into();
    cfg.albedo = cli.albedo;
    cfg.system = system;
    cfg.temp_air = cli.temp_air;
    cfg.wind = cli.wind;
    cfg.horizon.radius = cli.horizon_radius;
    cfg.horizon.directions = cli.horizon_dirs;
    cfg.latitude_mode = if cli.per_cell_lat {
        LatitudeMode::PerCellGeographic
    } else {
        LatitudeMode::Center
    };
    if cli.spa {
        cfg.spa = Some(solarpv_core::spa::SpaParams { elevation: 0.0, ..Default::default() });
    }
    cfg.iam = match cli.iam {
        IamKind::None => None,
        IamKind::Ashrae => Some(solarpv_core::losses::IamModel::ashrae()),
        IamKind::MartinRuiz => Some(solarpv_core::losses::IamModel::martin_ruiz()),
        IamKind::Physical => Some(solarpv_core::losses::IamModel::physical()),
    };
    cfg.pressure_pa = cli.pressure;
    cfg.spectral = if cli.spectral { Some(SpectralLoss::c_si()) } else { None };
    cfg.apply_sky_view_factor = cli.svf;
    cfg.mount = match cli.mount {
        MountKind::Terrain => Mount::FixedTerrain,
        MountKind::Tilt => Mount::FixedTilt {
            tilt: cli.tilt,
            surface_azimuth: cli.surface_azimuth,
        },
        MountKind::Tracker => Mount::SingleAxis(SingleAxisTracker { gcr: cli.gcr, ..Default::default() }),
        MountKind::DualAxis => Mount::DualAxis(solarpv_core::tracking::DualAxisTracker { max_tilt: cli.max_tilt }),
    };

    let (res, unit) = if let Some(ref wpath) = cli.weather {
        let records = read_weather_csv(wpath)?;
        eprintln!(
            "Computing PV potential from {} weather records ({:?} sky, dt={} h)…",
            records.len(), cli.sky, cli.weather_dt
        );
        let r = pv_potential_series(&dem, &cfg, &records, cli.weather_dt)
            .map_err(|e| anyhow::anyhow!(e))?;
        eprintln!("Integrated {} daylight records.", r.sun_steps);
        (r, "period")
    } else if cli.annual {
        let sampling = if cli.day_stride == 0 {
            DaySampling::MonthlyRepresentative
        } else {
            DaySampling::EveryNDays(cli.day_stride)
        };
        eprintln!(
            "Computing ANNUAL PV potential for {} ({:?} sky, {:?})…",
            cfg.date.year, cli.sky, sampling
        );
        let r = pv_potential_annual(&dem, &cfg, sampling).map_err(|e| anyhow::anyhow!(e))?;
        eprintln!("Integrated {} representative days.", r.sun_steps);
        (r, "year")
    } else {
        eprintln!("Computing PV potential for {} ({:?} sky model)…", cli.date, cli.sky);
        let r = pv_potential(&dem, &cfg).map_err(|e| anyhow::anyhow!(e))?;
        eprintln!("Integrated {} daylight steps.", r.sun_steps);
        (r, "day")
    };

    let labels = [
        ("poa", &res.poa_wh, format!("POA insolation (Wh/m²/{unit})")),
        ("ac", &res.ac_wh, format!("AC energy (Wh/{unit})")),
        ("specific_yield", &res.specific_yield, format!("specific yield (kWh/kWp/{unit})")),
    ];
    let outputs: Vec<_> = labels.iter().map(|(s, r, l)| (*s, *r, l.as_str())).collect();
    for (suffix, raster, label) in outputs {
        let path = format!("{}_{}.tif", cli.out_prefix, suffix);
        write_geotiff(raster, &path, None)
            .with_context(|| format!("writing {path}"))?;
        eprintln!("  wrote {path} — {label}");
    }

    Ok(())
}
