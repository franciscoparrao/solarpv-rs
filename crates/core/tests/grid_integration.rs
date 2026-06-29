//! Integration tests for the gridded PV potential (feature `terrain`).
//!
//! Uses small synthetic DEMs so the geometry is known a priori:
//!   * a horizontal plane should reproduce the point-model order of magnitude;
//!   * on an E–W ridge in the Southern Hemisphere the North-facing (equator-
//!     facing) slope must out-yield the symmetric South-facing slope.
#![cfg(feature = "terrain")]

use solarpv_core::grid::{pv_potential, GridConfig};
use solarpv_core::solpos::{DateTimeUtc, Location};
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
