//! Quick scale check / usage example for the gridded PV potential.
//!
//! Run with:
//!     cargo run --release --example grid_bench --features terrain
//!
//! Builds a synthetic ridged DEM, computes the daily clear-sky PV potential over
//! the Atacama point and reports timing plus a few summary statistics.

#[cfg(not(feature = "terrain"))]
fn main() {
    eprintln!("rebuild with `--features terrain`");
}

#[cfg(feature = "terrain")]
fn main() {
    use solarpv_core::grid::{pv_potential, GridConfig};
    use solarpv_core::solpos::{DateTimeUtc, Location};
    use std::time::Instant;
    use surtgis_core::{GeoTransform, Raster};

    let n = 600;
    // A diagonal ripple so cells take a spread of slopes and aspects.
    let mut data = vec![0.0f64; n * n];
    for r in 0..n {
        for c in 0..n {
            let rr = r as f64;
            let cc = c as f64;
            data[r * n + c] = 50.0 * ((rr / 40.0).sin() + (cc / 55.0).cos()) + 0.05 * rr;
        }
    }
    let mut dem = Raster::from_vec(data, n, n).unwrap();
    dem.set_transform(GeoTransform::new(-69.0, -23.0, 30.0, -30.0));

    let center = Location::new(-23.0, -69.0).unwrap();
    let date = DateTimeUtc::new(2026, 6, 21, 0, 0, 0).unwrap();
    let mut cfg = GridConfig::new(center, date);
    cfg.time_step_minutes = 15;

    let t0 = Instant::now();
    let res = pv_potential(&dem, &cfg).unwrap();
    let dt = t0.elapsed();

    let sy: Vec<f64> = (0..n * n)
        .map(|i| res.specific_yield.get(i / n, i % n).unwrap())
        .collect();
    let mean = sy.iter().sum::<f64>() / sy.len() as f64;
    let max = sy.iter().cloned().fold(f64::MIN, f64::max);
    let min = sy.iter().cloned().fold(f64::MAX, f64::min);

    println!(
        "{n}x{n} DEM, {} daylight steps @ {}-min → {:.2?}",
        res.sun_steps, cfg.time_step_minutes, dt
    );
    println!("specific yield kWh/kWp/day: min {min:.2}, mean {mean:.2}, max {max:.2}");
}
