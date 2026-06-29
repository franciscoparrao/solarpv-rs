//! Solar geometry: position of the sun and angle of incidence on a tilted plane.
//!
//! The position algorithm is the Astronomical Almanac's low-precision formula
//! (Michalsky, 1988, *Solar Energy* 40(3):227–235), the same family used by
//! `pvlib`'s `ephemeris` method. Accuracy is on the order of 0.01° for dates near
//! the present, which is far below the sensitivity of downstream irradiance and
//! PV yield. The full NREL SPA (Reda & Andreas, 2004) is planned for v0.2 when
//! sub-arcminute precision is needed.
//!
//! Atmospheric refraction is applied with the standard piecewise model so that
//! [`SolarPosition::apparent_zenith`] matches what `pvlib` reports.

use crate::error::{Error, Result};

/// A UTC instant, expressed as civil calendar fields.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DateTimeUtc {
    /// Full year, e.g. 2026.
    pub year: i32,
    /// Month 1–12.
    pub month: u32,
    /// Day of month 1–31.
    pub day: u32,
    /// Hour 0–23 (UTC).
    pub hour: u32,
    /// Minute 0–59.
    pub minute: u32,
    /// Second 0–59 (fractional seconds not modelled).
    pub second: u32,
}

impl DateTimeUtc {
    /// Construct a UTC datetime, validating the calendar fields.
    pub fn new(year: i32, month: u32, day: u32, hour: u32, minute: u32, second: u32) -> Result<Self> {
        if !(1..=12).contains(&month) {
            return Err(Error::InvalidDate(format!("month {month} out of range")));
        }
        if !(1..=31).contains(&day) {
            return Err(Error::InvalidDate(format!("day {day} out of range")));
        }
        if hour > 23 || minute > 59 || second > 59 {
            return Err(Error::InvalidDate(format!(
                "time {hour:02}:{minute:02}:{second:02} out of range"
            )));
        }
        Ok(Self { year, month, day, hour, minute, second })
    }

    /// Day of year in `1..=366`.
    pub fn day_of_year(&self) -> u32 {
        const CUM: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
        let leap = (self.year % 4 == 0 && self.year % 100 != 0) || self.year % 400 == 0;
        let extra = if leap && self.month > 2 { 1 } else { 0 };
        CUM[(self.month - 1) as usize] + self.day + extra
    }

    /// Fractional hour of day in UTC, in `[0, 24)`.
    fn hour_fraction(&self) -> f64 {
        self.hour as f64 + self.minute as f64 / 60.0 + self.second as f64 / 3600.0
    }

    /// Julian Day (including the fractional day), valid for the Gregorian
    /// calendar. Meeus, *Astronomical Algorithms*, ch. 7.
    pub fn julian_day(&self) -> f64 {
        let (mut y, mut m) = (self.year, self.month as i32);
        if m <= 2 {
            y -= 1;
            m += 12;
        }
        let a = (y as f64 / 100.0).floor();
        let b = 2.0 - a + (a / 4.0).floor();
        let day_frac = self.day as f64 + self.hour_fraction() / 24.0;
        (365.25 * (y as f64 + 4716.0)).floor()
            + (30.6001 * (m as f64 + 1.0)).floor()
            + day_frac
            + b
            - 1524.5
    }
}

/// Geographic location on the WGS-84 ellipsoid (height ignored for geometry).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Location {
    /// Latitude in degrees, north positive (`[-90, 90]`).
    pub latitude: f64,
    /// Longitude in degrees, east positive (`[-180, 180]`).
    pub longitude: f64,
}

impl Location {
    /// Build a location, validating the angular ranges.
    pub fn new(latitude: f64, longitude: f64) -> Result<Self> {
        if !(-90.0..=90.0).contains(&latitude) {
            return Err(Error::param("latitude", latitude, "[-90, 90] degrees"));
        }
        if !(-180.0..=180.0).contains(&longitude) {
            return Err(Error::param("longitude", longitude, "[-180, 180] degrees"));
        }
        Ok(Self { latitude, longitude })
    }
}

/// Position of the sun in the local sky.
///
/// Azimuth is measured clockwise from North (0° = N, 90° = E, 180° = S). Both
/// the geometric (`*_zenith`/`elevation`) and refraction-corrected
/// (`apparent_*`) elevations are provided; downstream transposition uses the
/// apparent zenith, matching `pvlib`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SolarPosition {
    /// Geometric zenith angle in degrees (90° − true elevation).
    pub zenith: f64,
    /// Refraction-corrected (apparent) zenith angle in degrees.
    pub apparent_zenith: f64,
    /// Geometric elevation above the horizon in degrees.
    pub elevation: f64,
    /// Refraction-corrected (apparent) elevation in degrees.
    pub apparent_elevation: f64,
    /// Azimuth in degrees, clockwise from North in `[0, 360)`.
    pub azimuth: f64,
    /// Equation of time in minutes (apparent − mean solar time).
    pub equation_of_time: f64,
    /// Solar declination in degrees.
    pub declination: f64,
}

const DEG: f64 = std::f64::consts::PI / 180.0;
const RAD: f64 = 180.0 / std::f64::consts::PI;

/// Atmospheric refraction correction (degrees) added to true elevation.
///
/// Standard piecewise model (Astronomical Almanac / `pvlib` `ephemeris`),
/// returning arc-seconds internally then converting to degrees. `elev_deg` is
/// the *true* (geometric) elevation in degrees.
fn refraction_correction(elev_deg: f64) -> f64 {
    if elev_deg > 85.0 {
        return 0.0;
    }
    let te = (elev_deg * DEG).tan();
    let arcsec = if elev_deg > 5.0 {
        58.1 / te - 0.07 / te.powi(3) + 0.000_086 / te.powi(5)
    } else if elev_deg > -0.575 {
        let e = elev_deg;
        1735.0 + e * (-518.2 + e * (103.4 + e * (-12.79 + e * 0.711)))
    } else {
        -20.774 / te
    };
    arcsec / 3600.0
}

/// Compute the solar position for a UTC instant at a location.
pub fn solar_position(when: DateTimeUtc, loc: Location) -> SolarPosition {
    let jd = when.julian_day();
    // Days since the J2000.0 epoch.
    let n = jd - 2_451_545.0;

    // Mean longitude and mean anomaly of the sun (degrees).
    let mean_long = (280.460 + 0.985_647_4 * n).rem_euclid(360.0);
    let mean_anom = ((357.528 + 0.985_600_3 * n).rem_euclid(360.0)) * DEG;

    // Apparent ecliptic longitude (degrees), with equation-of-centre terms.
    let ecl_long = (mean_long + 1.915 * mean_anom.sin() + 0.020 * (2.0 * mean_anom).sin()) * DEG;
    // Obliquity of the ecliptic (degrees → radians).
    let obliquity = (23.439 - 0.000_000_4 * n) * DEG;

    // Right ascension and declination.
    let right_asc = (obliquity.cos() * ecl_long.sin()).atan2(ecl_long.cos());
    let declination = (obliquity.sin() * ecl_long.sin()).asin();

    // Greenwich mean sidereal time (hours) → local mean sidereal time (degrees).
    let gmst = (6.697_375 + 0.065_709_824_2 * n + when.hour_fraction()).rem_euclid(24.0);
    let lmst_deg = (gmst * 15.0 + loc.longitude).rem_euclid(360.0);

    // Local hour angle (radians), positive towards the west.
    let hour_angle = lmst_deg * DEG - right_asc;

    let lat = loc.latitude * DEG;
    let sin_elev =
        declination.sin() * lat.sin() + declination.cos() * lat.cos() * hour_angle.cos();
    let elevation = sin_elev.clamp(-1.0, 1.0).asin();

    // Azimuth measured clockwise from North.
    let sin_az = -declination.cos() * hour_angle.sin() / elevation.cos();
    let cos_az = (declination.sin() - lat.sin() * elevation.sin()) / (lat.cos() * elevation.cos());
    let azimuth = sin_az.atan2(cos_az).rem_euclid(2.0 * std::f64::consts::PI);

    let elevation_deg = elevation * RAD;
    let refr = refraction_correction(elevation_deg);
    let apparent_elev = elevation_deg + refr;

    // Equation of time in minutes: 4 min per degree of (mean_long − RA).
    let ra_deg = (right_asc * RAD).rem_euclid(360.0);
    let mut eot = 4.0 * (mean_long - ra_deg);
    // Fold into the conventional ±20 min window.
    if eot > 20.0 {
        eot -= 1440.0;
    } else if eot < -20.0 {
        eot += 1440.0;
    }

    SolarPosition {
        zenith: 90.0 - elevation_deg,
        apparent_zenith: 90.0 - apparent_elev,
        elevation: elevation_deg,
        apparent_elevation: apparent_elev,
        azimuth: azimuth * RAD,
        equation_of_time: eot,
        declination: declination * RAD,
    }
}

/// Angle of incidence (degrees) of the beam on a tilted plane.
///
/// `surface_tilt` is measured from the horizontal (0° = horizontal, 90° =
/// vertical); `surface_azimuth` is the direction the plane faces, clockwise from
/// North (180° = facing the equator in the Southern Hemisphere is *North-facing*,
/// i.e. `surface_azimuth = 0`). `solar_zenith`/`solar_azimuth` are the apparent
/// position from [`solar_position`].
///
/// Mirrors `pvlib.irradiance.aoi`.
pub fn angle_of_incidence(
    surface_tilt: f64,
    surface_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
) -> f64 {
    let st = surface_tilt * DEG;
    let sz = solar_zenith * DEG;
    let cos_aoi = sz.cos() * st.cos()
        + sz.sin() * st.sin() * ((solar_azimuth - surface_azimuth) * DEG).cos();
    cos_aoi.clamp(-1.0, 1.0).acos() * RAD
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn julian_day_reference_epoch() {
        // 2000-01-01 12:00:00 UTC is JD 2451545.0 by definition.
        let dt = DateTimeUtc::new(2000, 1, 1, 12, 0, 0).unwrap();
        assert!((dt.julian_day() - 2_451_545.0).abs() < 1e-6);
    }

    #[test]
    fn solar_noon_is_high_and_southward_in_nh_summer() {
        // Greenwich, ~solar noon on the summer solstice: sun high and due south.
        let dt = DateTimeUtc::new(2026, 6, 21, 12, 0, 0).unwrap();
        let loc = Location::new(51.48, 0.0).unwrap();
        let p = solar_position(dt, loc);
        assert!(p.elevation > 60.0, "elevation {}", p.elevation);
        assert!((p.azimuth - 180.0).abs() < 5.0, "azimuth {}", p.azimuth);
    }

    #[test]
    fn aoi_zero_when_plane_faces_sun() {
        // A plane pointed straight at the sun has zero angle of incidence.
        let aoi = angle_of_incidence(30.0, 137.0, 30.0, 137.0);
        assert!(aoi < 1e-6, "aoi {aoi}");
    }
}
