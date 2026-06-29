//! `solarpv` — command-line interface for the terrain PV-potential engine.
//!
//! Reads a DEM GeoTIFF, computes the daily clear-sky PV potential over it and
//! writes per-cell POA insolation, AC energy and specific-yield rasters.

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};

use solarpv_core::grid::{pv_potential, GridConfig, LatitudeMode};
use solarpv_core::irradiance::SkyModel;
use solarpv_core::pv::PvSystem;
use solarpv_core::solpos::{DateTimeUtc, Location};
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

    /// Use per-cell latitude/longitude from the DEM transform (requires a
    /// geographic lon/lat DEM, e.g. EPSG:4326). Default: scene-centre only.
    #[arg(long, default_value_t = false)]
    per_cell_lat: bool,
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

    eprintln!("Computing PV potential for {} ({} sky model)…", cli.date, format!("{:?}", cli.sky));
    let res = pv_potential(&dem, &cfg).map_err(|e| anyhow::anyhow!(e))?;
    eprintln!("Integrated {} daylight steps.", res.sun_steps);

    let outputs = [
        ("poa", &res.poa_wh, "POA insolation (Wh/m²/day)"),
        ("ac", &res.ac_wh, "AC energy (Wh/day)"),
        ("specific_yield", &res.specific_yield, "specific yield (kWh/kWp/day)"),
    ];
    for (suffix, raster, label) in outputs {
        let path = format!("{}_{}.tif", cli.out_prefix, suffix);
        write_geotiff(raster, &path, None)
            .with_context(|| format!("writing {path}"))?;
        eprintln!("  wrote {path} — {label}");
    }

    Ok(())
}
