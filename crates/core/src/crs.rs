//! Coordinate helpers for the gridded engine.
//!
//! The terrain step needs a **projected** DEM in metres (so that slope, aspect
//! and horizon distances are meaningful), but per-cell solar geometry needs
//! **geographic** latitude/longitude. When the DEM is a UTM grid we bridge the
//! two analytically, without pulling in PROJ: [`utm_zone_from_epsg`] recognises
//! the standard WGS84 UTM EPSG codes and [`inverse_utm`] maps easting/northing
//! back to latitude/longitude on the WGS84 ellipsoid (Snyder's inverse
//! transverse-Mercator series, sub-metre over a UTM zone).
//!
//! This exists because reading a UTM northing (millions of metres) as if it
//! were a latitude in degrees silently produces a per-row-striped, physically
//! meaningless yield map — see the regression test in `grid_integration.rs`.

// WGS84 ellipsoid.
const A: f64 = 6_378_137.0; // semi-major axis (m)
const F: f64 = 1.0 / 298.257_223_563; // flattening
const K0: f64 = 0.9996; // UTM scale factor on the central meridian
const FALSE_EASTING: f64 = 500_000.0;
const FALSE_NORTHING_SOUTH: f64 = 10_000_000.0;

/// Recognise a WGS84 UTM CRS from its EPSG code.
///
/// Returns `(zone, is_northern_hemisphere)` for the standard WGS84 UTM codes:
/// `326zz` for the northern hemisphere and `327zz` for the southern, with
/// `zz` the zone 1–60 (e.g. `32719` → zone 19 South, the Chilean north).
/// Any other code returns `None`.
pub fn utm_zone_from_epsg(epsg: u32) -> Option<(u8, bool)> {
    let (base, north) = match epsg / 100 {
        326 => (32_600, true),
        327 => (32_700, false),
        _ => return None,
    };
    let zone = epsg - base;
    if (1..=60).contains(&zone) {
        Some((zone as u8, north))
    } else {
        None
    }
}

/// Central meridian of a UTM zone, in degrees.
fn central_meridian(zone: u8) -> f64 {
    zone as f64 * 6.0 - 183.0
}

/// Inverse UTM: map easting/northing (metres) to `(latitude, longitude)` in
/// degrees on the WGS84 ellipsoid.
///
/// `zone` is the UTM zone (1–60) and `north` selects the hemisphere (which sets
/// the false northing). Uses Snyder's inverse transverse-Mercator series
/// (USGS Professional Paper 1395), accurate to well under a metre anywhere
/// inside a UTM zone — far tighter than the DEM resolution.
pub fn inverse_utm(easting: f64, northing: f64, zone: u8, north: bool) -> (f64, f64) {
    let e2 = F * (2.0 - F); // first eccentricity squared
    let ep2 = e2 / (1.0 - e2); // second eccentricity squared

    let x = easting - FALSE_EASTING;
    let y = if north { northing } else { northing - FALSE_NORTHING_SOUTH };

    let m = y / K0;
    let mu = m / (A * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0));

    let e1 = (1.0 - (1.0 - e2).sqrt()) / (1.0 + (1.0 - e2).sqrt());
    let e1_2 = e1 * e1;
    let e1_3 = e1_2 * e1;
    let e1_4 = e1_3 * e1;

    // Footpoint latitude.
    let phi1 = mu
        + (3.0 * e1 / 2.0 - 27.0 * e1_3 / 32.0) * (2.0 * mu).sin()
        + (21.0 * e1_2 / 16.0 - 55.0 * e1_4 / 32.0) * (4.0 * mu).sin()
        + (151.0 * e1_3 / 96.0) * (6.0 * mu).sin()
        + (1097.0 * e1_4 / 512.0) * (8.0 * mu).sin();

    let sin_phi1 = phi1.sin();
    let cos_phi1 = phi1.cos();
    let tan_phi1 = phi1.tan();

    let c1 = ep2 * cos_phi1 * cos_phi1;
    let t1 = tan_phi1 * tan_phi1;
    let n1 = A / (1.0 - e2 * sin_phi1 * sin_phi1).sqrt();
    let r1 = A * (1.0 - e2) / (1.0 - e2 * sin_phi1 * sin_phi1).powf(1.5);
    let d = x / (n1 * K0);

    let d2 = d * d;
    let d3 = d2 * d;
    let d4 = d3 * d;
    let d5 = d4 * d;
    let d6 = d5 * d;

    let lat = phi1
        - (n1 * tan_phi1 / r1)
            * (d2 / 2.0
                - (5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 * c1 * c1 - 9.0 * ep2) * d4 / 24.0
                + (61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1 * t1 - 252.0 * ep2 - 3.0 * c1 * c1)
                    * d6
                    / 720.0);

    let lon_rad = (d - (1.0 + 2.0 * t1 + c1) * d3 / 6.0
        + (5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1 * c1 + 8.0 * ep2 + 24.0 * t1 * t1) * d5 / 120.0)
        / cos_phi1;

    let lon = central_meridian(zone) + lon_rad.to_degrees();
    (lat.to_degrees(), lon)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognises_wgs84_utm_codes() {
        assert_eq!(utm_zone_from_epsg(32719), Some((19, false))); // Chilean north
        assert_eq!(utm_zone_from_epsg(32619), Some((19, true))); // northern mirror
        assert_eq!(utm_zone_from_epsg(32601), Some((1, true)));
        assert_eq!(utm_zone_from_epsg(32760), Some((60, false)));
        assert_eq!(utm_zone_from_epsg(4326), None); // geographic
        assert_eq!(utm_zone_from_epsg(3857), None); // web mercator
        assert_eq!(utm_zone_from_epsg(32761), None); // zone 61 does not exist
        assert_eq!(utm_zone_from_epsg(32700), None); // zone 0 does not exist
    }

    #[test]
    fn central_meridian_of_easting_maps_to_exact_longitude() {
        // On the central meridian (easting = false easting) the longitude is
        // exactly the zone's central meridian, independent of northing.
        let (_, lon) = inverse_utm(FALSE_EASTING, 7_400_000.0, 19, false);
        assert!((lon - (-69.0)).abs() < 1e-9, "lon {lon} should be -69.0 exactly");
    }

    #[test]
    fn inverts_a_known_atacama_point() {
        // Zone 19S, points in the Chilean north. Reference latitude/longitude
        // from pyproj (EPSG:32719 → EPSG:4326); the series matches to well
        // under a metre, so 1e-5° (~1 m) is a tight check.
        let (lat, lon) = inverse_utm(500_000.0, 7_400_000.0, 19, false);
        assert!((lat - (-23.510195)).abs() < 1e-5, "lat {lat}");
        assert!((lon - (-69.0)).abs() < 1e-9, "lon {lon}");

        // Off the central meridian: E=360000 (140 km west), N=7390000.
        let (lat2, lon2) = inverse_utm(360_000.0, 7_390_000.0, 19, false);
        assert!((lat2 - (-23.594461)).abs() < 1e-5, "lat {lat2}");
        assert!((lon2 - (-70.372093)).abs() < 1e-5, "lon {lon2}");
    }
}
