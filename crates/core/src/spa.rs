//! NREL Solar Position Algorithm (Reda & Andreas, 2004; NREL/TP-560-34302).
//!
//! High-accuracy (±0.0003°) solar position, the reference used by
//! `pvlib.solarposition.spa_python`. Heavier than the Michalsky model in
//! [`crate::solpos`]; use it when sub-arcminute precision matters. The periodic
//! term tables live in `spa_tables.rs`, generated verbatim from pvlib.

use crate::solpos::{DateTimeUtc, Location, SolarPosition};
use crate::spa_tables as t;

const DEG: f64 = std::f64::consts::PI / 180.0;
const RAD: f64 = 180.0 / std::f64::consts::PI;

/// Atmospheric parameters for the SPA topocentric correction and refraction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpaParams {
    /// Difference TT − UT1 in seconds (ΔT). pvlib default 67.0.
    pub delta_t: f64,
    /// Site elevation above sea level, metres.
    pub elevation: f64,
    /// Local pressure in pascals (used for refraction).
    pub pressure_pa: f64,
    /// Local air temperature in °C (used for refraction).
    pub temperature_c: f64,
    /// Atmospheric refraction at the horizon, degrees (pvlib default 0.5667).
    pub atmos_refract: f64,
}

impl Default for SpaParams {
    fn default() -> Self {
        // pvlib spa_python defaults.
        Self { delta_t: 67.0, elevation: 0.0, pressure_pa: 101_325.0, temperature_c: 12.0, atmos_refract: 0.5667 }
    }
}

/// Sum a periodic-term table `Σ A·cos(B + C·jme)`.
fn periodic_sum(table: &[[f64; 3]], jme: f64) -> f64 {
    table.iter().map(|r| r[0] * (r[1] + r[2] * jme).cos()).sum()
}

/// Heliocentric series evaluated as `(t0 + t1·jme + … )/1e8`.
fn series5(tables: &[&[[f64; 3]]], jme: f64) -> f64 {
    let mut acc = 0.0;
    let mut p = 1.0;
    for tab in tables {
        acc += periodic_sum(tab, jme) * p;
        p *= jme;
    }
    acc / 1e8
}

/// Compute the solar position with the NREL SPA.
pub fn solar_position_spa(when: DateTimeUtc, loc: Location, p: &SpaParams) -> SolarPosition {
    let jd = when.julian_day();
    let jde = jd + p.delta_t / 86400.0;
    let jc = (jd - 2_451_545.0) / 36525.0;
    let jce = (jde - 2_451_545.0) / 36525.0;
    let jme = jce / 10.0;

    // Heliocentric longitude L (rad → deg), latitude B (rad → deg), radius R (AU).
    let l_rad = series5(&[&t::L0, &t::L1, &t::L2, &t::L3, &t::L4, &t::L5], jme);
    let b_rad = series5(&[&t::B0, &t::B1], jme);
    let r = series5(&[&t::R0, &t::R1, &t::R2, &t::R3, &t::R4], jme);
    let l_deg = (l_rad * RAD).rem_euclid(360.0);
    let b_deg = b_rad * RAD;

    // Geocentric longitude Θ and latitude β.
    let theta = (l_deg + 180.0).rem_euclid(360.0);
    let beta = -b_deg;

    // Nutation in longitude (Δψ) and obliquity (Δε), degrees.
    let x = [
        297.85036 + 445267.111480 * jce - 0.0019142 * jce.powi(2) + jce.powi(3) / 189474.0,
        357.52772 + 35999.050340 * jce - 0.0001603 * jce.powi(2) - jce.powi(3) / 300000.0,
        134.96298 + 477198.867398 * jce + 0.0086972 * jce.powi(2) + jce.powi(3) / 56250.0,
        93.27191 + 483202.017538 * jce - 0.0036825 * jce.powi(2) + jce.powi(3) / 327270.0,
        125.04452 - 1934.136261 * jce + 0.0020708 * jce.powi(2) + jce.powi(3) / 450000.0,
    ];
    let (mut dpsi, mut deps) = (0.0, 0.0);
    for i in 0..63 {
        let mut arg = 0.0;
        for j in 0..5 {
            arg += x[j] * t::NUTATION_YTERM[i][j];
        }
        let arg = (arg * DEG).sin_cos();
        let (sin_arg, cos_arg) = (arg.0, arg.1);
        let abcd = t::NUTATION_ABCD[i];
        dpsi += (abcd[0] + abcd[1] * jce) * sin_arg;
        deps += (abcd[2] + abcd[3] * jce) * cos_arg;
    }
    dpsi /= 36_000_000.0;
    deps /= 36_000_000.0;

    // Mean and true obliquity of the ecliptic (degrees).
    let u = jme / 10.0;
    let eps0 = 84381.448 - 4680.93 * u - 1.55 * u.powi(2) + 1999.25 * u.powi(3)
        - 51.38 * u.powi(4) - 249.67 * u.powi(5) - 39.05 * u.powi(6) + 7.12 * u.powi(7)
        + 27.87 * u.powi(8) + 5.79 * u.powi(9) + 2.45 * u.powi(10);
    let eps = eps0 / 3600.0 + deps;

    // Aberration and apparent sun longitude.
    let dtau = -20.4898 / (3600.0 * r);
    let lambda = theta + dpsi + dtau;

    // Apparent sidereal time at Greenwich (degrees).
    let nu0 = (280.46061837 + 360.98564736629 * (jd - 2_451_545.0) + 0.000387933 * jc.powi(2)
        - jc.powi(3) / 38_710_000.0)
        .rem_euclid(360.0);
    let nu = nu0 + dpsi * (eps * DEG).cos();

    // Geocentric right ascension α and declination δ (degrees).
    let (lam, bet, ep) = (lambda * DEG, beta * DEG, eps * DEG);
    let alpha = (lam.sin() * ep.cos() - bet.tan() * ep.sin())
        .atan2(lam.cos())
        .to_degrees()
        .rem_euclid(360.0);
    let delta = (bet.sin() * ep.cos() + bet.cos() * ep.sin() * lam.sin()).asin() * RAD;

    // Observer local hour angle H (degrees).
    let h = (nu + loc.longitude - alpha).rem_euclid(360.0);

    // Topocentric corrections (equatorial horizontal parallax).
    let xi = 8.794 / (3600.0 * r); // degrees
    let lat = loc.latitude * DEG;
    let uu = (0.99664719 * lat.tan()).atan();
    let xterm = uu.cos() + (p.elevation / 6_378_140.0) * lat.cos();
    let yterm = 0.99664719 * uu.sin() + (p.elevation / 6_378_140.0) * lat.sin();
    let (xi_r, h_r, d_r) = (xi * DEG, h * DEG, delta * DEG);
    let dalpha = (-xterm * xi_r.sin() * h_r.sin())
        .atan2(d_r.cos() - xterm * xi_r.sin() * h_r.cos())
        * RAD;
    let delta_prime =
        ((d_r.sin() - yterm * xi_r.sin()) * (dalpha * DEG).cos()).atan2(d_r.cos() - xterm * xi_r.sin() * h_r.cos()) * RAD;
    let h_prime = h - dalpha;

    // Topocentric elevation without refraction, then refraction correction.
    let (hp_r, dp_r) = (h_prime * DEG, delta_prime * DEG);
    let e0 = (lat.sin() * dp_r.sin() + lat.cos() * dp_r.cos() * hp_r.cos()).asin() * RAD;
    let sun_radius = 0.26667;
    let del_e = if e0 >= -(sun_radius + p.atmos_refract) {
        let pressure_mbar = p.pressure_pa / 100.0;
        (pressure_mbar / 1010.0) * (283.0 / (273.0 + p.temperature_c)) * 1.02
            / (60.0 * ((e0 + 10.3 / (e0 + 5.11)) * DEG).tan())
    } else {
        0.0
    };
    let e = e0 + del_e;

    // Topocentric azimuth (measured clockwise from North).
    let gamma = hp_r.sin().atan2(hp_r.cos() * lat.sin() - dp_r.tan() * lat.cos()) * RAD;
    let azimuth = (gamma + 180.0).rem_euclid(360.0);

    // Equation of time (minutes).
    let m = (280.4664567 + 360007.6982779 * jme + 0.03032028 * jme.powi(2) + jme.powi(3) / 49931.0
        - jme.powi(4) / 15300.0 - jme.powi(5) / 2_000_000.0)
        .rem_euclid(360.0);
    let mut eot = 4.0 * (m - 0.0057183 - alpha + dpsi * (eps * DEG).cos());
    eot = eot.rem_euclid(1440.0);
    if eot > 20.0 {
        eot -= 1440.0;
    } else if eot < -20.0 {
        eot += 1440.0;
    }

    SolarPosition {
        zenith: 90.0 - e0,
        apparent_zenith: 90.0 - e,
        elevation: e0,
        apparent_elevation: e,
        azimuth,
        equation_of_time: eot,
        declination: delta_prime,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn agrees_with_michalsky_within_a_tenth_degree() {
        // SPA and the Michalsky model should be close (sub-0.1°) for a normal date.
        let when = DateTimeUtc::new(2026, 3, 21, 15, 30, 0).unwrap();
        let loc = Location::new(-23.0, -69.0).unwrap();
        let spa = solar_position_spa(when, loc, &SpaParams::default());
        let mich = crate::solpos::solar_position(when, loc);
        assert!((spa.apparent_zenith - mich.apparent_zenith).abs() < 0.1);
    }
}
