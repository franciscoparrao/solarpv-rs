//! Integration tests for the gridded PV potential (feature `terrain`).
//!
//! Uses small synthetic DEMs so the geometry is known a priori:
//!   * a horizontal plane should reproduce the point-model order of magnitude;
//!   * on an E–W ridge in the Southern Hemisphere the North-facing (equator-
//!     facing) slope must out-yield the symmetric South-facing slope.
#![cfg(feature = "terrain")]

use solarpv_core::grid::{
    pv_potential, pv_potential_annual, pv_potential_series, DaySampling, GridConfig, LatitudeMode,
    Mount, WeatherRecord,
};
use solarpv_core::irradiance::haurwitz_clearsky_ghi;
use solarpv_core::losses::{IamModel, SpectralLoss};
use solarpv_core::solpos::{solar_position, DateTimeUtc, Location};
use solarpv_core::tracking::{DualAxisTracker, SingleAxisTracker};
use surtgis_core::{GeoTransform, Raster};

const CELL: f64 = 30.0;

fn with_transform(mut r: Raster<f64>) -> Raster<f64> {
    // North-up: origin at top-left, negative pixel height.
    r.set_transform(GeoTransform::new(-69.0, -23.0, CELL, -CELL));
    r
}

/// A flat DEM at constant elevation.
fn flat_dem(n: usize) -> Raster<f64> {
    with_transform(Raster::from_vec(vec![100.0; n * n], n, n).unwrap())
}

/// An E–W ridge: elevation peaks at the centre row and descends towards both
/// the northern (row 0) and southern (last row) edges.
fn ridge_dem(n: usize) -> Raster<f64> {
    let mid = (n / 2) as i64;
    let mut data = vec![0.0; n * n];
    for r in 0..n {
        let z = 100.0 - 10.0 * ((r as i64 - mid).abs() as f64);
        for c in 0..n {
            data[r * n + c] = z;
        }
    }
    with_transform(Raster::from_vec(data, n, n).unwrap())
}

/// A tall E–W wall to the North: the northern third is high, the rest is flat
/// low ground. Cells in the low zone have a high northern horizon.
fn north_wall_dem(n: usize) -> Raster<f64> {
    let mut data = vec![0.0; n * n];
    for r in 0..n {
        let z = if r < n / 3 { 100.0 } else { 0.0 };
        for c in 0..n {
            data[r * n + c] = z;
        }
    }
    with_transform(Raster::from_vec(data, n, n).unwrap())
}

fn config_on(date: DateTimeUtc) -> GridConfig {
    let center = Location::new(-23.0, -69.0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = 30; // keep the test fast
    cfg
}

/// Summer solstice: sun near the zenith over the Atacama (max yield).
fn summer_config() -> GridConfig {
    config_on(DateTimeUtc::new(2026, 12, 21, 0, 0, 0).unwrap())
}

/// Winter solstice: sun firmly in the northern sky, so North-facing slopes are
/// unambiguously equator-facing.
fn winter_config() -> GridConfig {
    config_on(DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap())
}

#[test]
fn flat_plane_yields_plausible_clear_sky_potential() {
    let dem = flat_dem(5);
    let res = pv_potential(&dem, &summer_config()).unwrap();

    assert!(res.sun_steps > 0, "no daylight steps");
    let poa = res.poa_wh.get(2, 2).unwrap();
    let sy = res.specific_yield.get(2, 2).unwrap();

    // Clear-sky summer day in the Atacama: a few kWh/m² of POA and a healthy
    // specific yield, but bounded well below the 24 h ceiling.
    assert!((5_000.0..=13_000.0).contains(&poa), "poa_wh {poa}");
    assert!((4.0..=11.0).contains(&sy), "specific_yield {sy} kWh/kWp");
}

#[test]
fn north_facing_slope_outyields_south_facing_in_southern_hemisphere() {
    let n = 11;
    let dem = ridge_dem(n);
    let res = pv_potential(&dem, &winter_config()).unwrap();

    let mid = n / 2;
    let col = n / 2;
    // Symmetric cells three rows either side of the ridge: row<mid faces North,
    // row>mid faces South.
    let north = res.poa_wh.get(mid - 3, col).unwrap();
    let south = res.poa_wh.get(mid + 3, col).unwrap();

    assert!(
        north > south * 1.05,
        "north-facing POA {north} should clearly exceed south-facing {south}"
    );
}

/// A flat geographic DEM spanning many degrees of latitude (north at row 0).
fn flat_geographic_lat_span(rows: usize, cols: usize, lat_n: f64, lat_s: f64) -> Raster<f64> {
    let mut r = Raster::from_vec(vec![100.0; rows * cols], rows, cols).unwrap();
    let pix_h = (lat_s - lat_n) / (rows as f64 - 1.0); // negative (lat decreases southward)
    r.set_transform(GeoTransform::new(-69.0, lat_n, 0.01, pix_h));
    r
}

#[test]
fn per_cell_latitude_varies_yield_across_a_wide_span() {
    // Flat plane from lat -10 (row 0) to lat -40 (last row), winter solstice.
    let (rows, cols) = (31, 5);
    let dem = flat_geographic_lat_span(rows, cols, -10.0, -40.0);
    let center = Location::new(-25.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();

    let mut per_cell = GridConfig::new(center, date);
    per_cell.time_step_minutes = 30;
    per_cell.latitude_mode = LatitudeMode::PerCellGeographic;
    let res = pv_potential(&dem, &per_cell).unwrap();

    let col = cols / 2;
    let north = res.poa_wh.get(0, col).unwrap(); // lat -10, near equator
    let south = res.poa_wh.get(rows - 1, col).unwrap(); // lat -40
    // In Southern-Hemisphere winter the equatorward cell receives clearly more.
    assert!(north > south * 1.2, "equatorward {north} vs poleward {south}");

    // Center mode ignores the span: rows are ~uniform.
    let mut center_mode = per_cell.clone();
    center_mode.latitude_mode = LatitudeMode::Center;
    let res_c = pv_potential(&dem, &center_mode).unwrap();
    let n_c = res_c.poa_wh.get(0, col).unwrap();
    let s_c = res_c.poa_wh.get(rows - 1, col).unwrap();
    assert!((n_c - s_c).abs() < 1.0, "center mode should be uniform: {n_c} vs {s_c}");
}

#[test]
fn annual_integration_is_plausible_and_sampling_agnostic() {
    let dem = flat_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    // Year is what matters for annual; start date day/month are ignored.
    let date = DateTimeUtc::new(2026, 1, 1, 0, 0, 0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = 30;

    let monthly = pv_potential_annual(&dem, &cfg, DaySampling::MonthlyRepresentative).unwrap();
    let weekly = pv_potential_annual(&dem, &cfg, DaySampling::EveryNDays(7)).unwrap();

    assert_eq!(monthly.sun_steps, 12, "12 representative days");

    let sy_m = monthly.specific_yield.get(2, 2).unwrap();
    let sy_w = weekly.specific_yield.get(2, 2).unwrap();

    // Annual clear-sky specific yield for the Atacama: well over 1000 kWh/kWp,
    // below the physical ceiling.
    assert!((1200.0..=3200.0).contains(&sy_m), "annual specific yield {sy_m} kWh/kWp");
    // The two sampling schemes should agree closely.
    let rel = (sy_m - sy_w).abs() / sy_m;
    assert!(rel < 0.05, "monthly {sy_m} vs weekly {sy_w} differ {:.1}%", rel * 100.0);

    // A single summer day must be a small fraction of the annual total.
    let day = pv_potential(&dem, &{
        let mut c = cfg.clone();
        c.date = DateTimeUtc::new(2026, 12, 21, 0, 0, 0).unwrap();
        c
    })
    .unwrap();
    let sy_day = day.specific_yield.get(2, 2).unwrap();
    assert!(sy_m > sy_day * 50.0, "annual {sy_m} vs single day {sy_day}");
}

#[test]
fn series_with_clearsky_ghi_reproduces_internal_path() {
    // Feeding pv_potential_series the same Haurwitz GHI the engine generates
    // internally must reproduce pv_potential exactly (same geometry, same Erbs).
    let dem = ridge_dem(7);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    let step = 30u32;
    cfg.time_step_minutes = step;

    let internal = pv_potential(&dem, &cfg).unwrap();

    // Build the matching Haurwitz series at the scene centre.
    let mut records = Vec::new();
    for i in 0..(24 * 60 / step) {
        let total_min = i * step;
        let when = DateTimeUtc::new(2026, 6, 21, total_min / 60, total_min % 60, 0).unwrap();
        let sun = solar_position(when, center);
        records.push(WeatherRecord {
            when,
            ghi: haurwitz_clearsky_ghi(sun.apparent_zenith),
            dni: None,
            dhi: None,
            temp_air: None,
            wind: None,
        });
    }
    let series = pv_potential_series(&dem, &cfg, &records, step as f64 / 60.0).unwrap();

    for r in 0..7 {
        for c in 0..7 {
            let a = internal.specific_yield.get(r, c).unwrap();
            let b = series.specific_yield.get(r, c).unwrap();
            assert!((a - b).abs() < 1e-6, "cell ({r},{c}): internal {a} vs series {b}");
        }
    }
}

#[test]
fn cloudy_series_yields_less_than_clearsky() {
    // Halving the irradiance series must reduce the yield.
    let dem = flat_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();
    let cfg = GridConfig::new(center, date);

    let mut clear = Vec::new();
    let mut cloudy = Vec::new();
    for i in 0..48u32 {
        let when = DateTimeUtc::new(2026, 6, 21, i * 30 / 60, i * 30 % 60, 0).unwrap();
        let g = haurwitz_clearsky_ghi(solar_position(when, center).apparent_zenith);
        clear.push(WeatherRecord { when, ghi: g, dni: None, dhi: None, temp_air: None, wind: None });
        cloudy.push(WeatherRecord { when, ghi: g * 0.5, dni: None, dhi: None, temp_air: None, wind: None });
    }
    let c = pv_potential_series(&dem, &cfg, &clear, 0.5).unwrap();
    let d = pv_potential_series(&dem, &cfg, &cloudy, 0.5).unwrap();
    assert!(d.specific_yield.get(2, 2).unwrap() < c.specific_yield.get(2, 2).unwrap());
}

#[test]
fn iam_reduces_yield_but_not_reported_poa() {
    // Enabling the IAM lowers AC yield (angular reflection) while the reported
    // geometric POA insolation is unchanged.
    let dem = ridge_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();
    let mut base = GridConfig::new(center, date);
    base.time_step_minutes = 30;
    base.mount = Mount::FixedTilt { tilt: 23.0, surface_azimuth: 0.0 };

    let mut withiam = base.clone();
    withiam.iam = Some(IamModel::physical());

    let plain = pv_potential(&dem, &base).unwrap();
    let modded = pv_potential(&dem, &withiam).unwrap();

    let (r, c) = (2, 2);
    assert!(
        modded.ac_wh.get(r, c).unwrap() < plain.ac_wh.get(r, c).unwrap(),
        "IAM should reduce AC energy"
    );
    assert!(
        (modded.poa_wh.get(r, c).unwrap() - plain.poa_wh.get(r, c).unwrap()).abs() < 1e-6,
        "geometric POA should be unchanged by IAM"
    );
}

#[test]
fn single_axis_tracker_outyields_fixed_horizontal() {
    // On a flat plane, a single-axis tracker should harvest clearly more annual
    // energy than fixed horizontal modules.
    let dem = flat_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 1, 1, 0, 0, 0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = 30;

    let mut fixed = cfg.clone();
    fixed.mount = Mount::FixedTilt { tilt: 0.0, surface_azimuth: 0.0 };
    let mut tracked = cfg.clone();
    tracked.mount = Mount::SingleAxis(SingleAxisTracker::default());

    let yr_fixed = pv_potential_annual(&dem, &fixed, DaySampling::MonthlyRepresentative).unwrap();
    let yr_track = pv_potential_annual(&dem, &tracked, DaySampling::MonthlyRepresentative).unwrap();

    let sy_fixed = yr_fixed.specific_yield.get(2, 2).unwrap();
    let sy_track = yr_track.specific_yield.get(2, 2).unwrap();
    assert!(sy_track > sy_fixed * 1.15, "tracker {sy_track} vs fixed-flat {sy_fixed}");
}

#[test]
fn dual_axis_tracker_outyields_single_axis_on_flat_plane() {
    // Ideal dual-axis tracking should harvest more annual energy than a
    // single-axis tracker because it follows both azimuth and elevation.
    let dem = flat_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 1, 1, 0, 0, 0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = 30;

    let mut single = cfg.clone();
    single.mount = Mount::SingleAxis(SingleAxisTracker::default());
    let mut dual = cfg.clone();
    dual.mount = Mount::DualAxis(DualAxisTracker::default());

    let yr_single = pv_potential_annual(&dem, &single, DaySampling::MonthlyRepresentative).unwrap();
    let yr_dual = pv_potential_annual(&dem, &dual, DaySampling::MonthlyRepresentative).unwrap();

    let sy_single = yr_single.specific_yield.get(2, 2).unwrap();
    let sy_dual = yr_dual.specific_yield.get(2, 2).unwrap();
    assert!(sy_dual > sy_single * 1.05, "dual-axis {sy_dual} should exceed single-axis {sy_single}");
}

#[test]
fn ridge_specific_yield_is_finite_and_positive_everywhere() {
    let dem = ridge_dem(7);
    let res = pv_potential(&dem, &summer_config()).unwrap();
    for r in 0..7 {
        for c in 0..7 {
            let sy = res.specific_yield.get(r, c).unwrap();
            assert!(sy.is_finite() && sy >= 0.0, "bad specific yield at ({r},{c}): {sy}");
        }
    }
}

#[test]
fn svf_is_neutral_on_flat_terrain() {
    // On a flat DEM the horizon is ~0 in every direction, so SVF ≈ 1 and the
    // result should be essentially unchanged.
    let dem = flat_dem(5);
    let mut base = winter_config();
    base.time_step_minutes = 30;
    let mut with_svf = base.clone();
    with_svf.apply_sky_view_factor = true;

    let plain = pv_potential(&dem, &base).unwrap();
    let modded = pv_potential(&dem, &with_svf).unwrap();
    let (r, c) = (2, 2);
    let rel = (modded.ac_wh.get(r, c).unwrap() - plain.ac_wh.get(r, c).unwrap()).abs()
        / plain.ac_wh.get(r, c).unwrap();
    assert!(rel < 0.01, "flat SVF changed yield by {rel:.4}");
}

#[test]
fn svf_reduces_diffuse_yield_in_wall_shadow() {
    // A high E–W wall to the North blocks a large fraction of the sky dome for
    // cells in the low southern zone. Enabling SVF lowers the sky-diffuse
    // component and therefore the AC yield.
    let n = 15;
    let dem = north_wall_dem(n);
    let mut base = winter_config();
    base.time_step_minutes = 30;
    base.mount = Mount::FixedTilt { tilt: 0.0, surface_azimuth: 0.0 };

    let mut with_svf = base.clone();
    with_svf.apply_sky_view_factor = true;

    let plain = pv_potential(&dem, &base).unwrap();
    let modded = pv_potential(&dem, &with_svf).unwrap();

    // Pick a low-zone cell well south of the wall.
    let (r, c) = (n - 3, n / 2);
    assert!(
        modded.ac_wh.get(r, c).unwrap() < plain.ac_wh.get(r, c).unwrap(),
        "SVF should reduce AC in wall shadow"
    );
    assert!(
        modded.poa_wh.get(r, c).unwrap() < plain.poa_wh.get(r, c).unwrap(),
        "SVF should reduce POA in wall shadow"
    );
}

#[test]
fn spectral_loss_reduces_yield_at_altitude() {
    // Atacama pressure (~756 hPa) increases absolute airmass and the c-Si SAPM
    // spectral factor is < 1 around noon, so enabling it lowers AC yield while
    // reported geometric POA stays unchanged.
    let dem = ridge_dem(5);
    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();
    let mut base = GridConfig::new(center, date);
    base.time_step_minutes = 30;
    base.pressure_pa = 75_626.0;

    let mut with_spectral = base.clone();
    with_spectral.spectral = Some(SpectralLoss::c_si());

    let plain = pv_potential(&dem, &base).unwrap();
    let modded = pv_potential(&dem, &with_spectral).unwrap();

    let (r, c) = (2, 2);
    assert!(
        modded.ac_wh.get(r, c).unwrap() < plain.ac_wh.get(r, c).unwrap(),
        "spectral loss should reduce AC energy"
    );
    assert!(
        (modded.poa_wh.get(r, c).unwrap() - plain.poa_wh.get(r, c).unwrap()).abs() < 1e-6,
        "geometric POA should be unchanged by spectral loss"
    );
}
