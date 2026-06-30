//! Single-axis tracker geometry.
//!
//! Computes the rotation and resulting surface orientation of a single-axis
//! tracker following the sun, with optional backtracking to avoid row-to-row
//! shading. Matches `pvlib.tracking.singleaxis` (Lorenzo et al., 2011); the
//! dominant configuration for utility-scale PV in the Atacama is the horizontal
//! North–South axis (`axis_tilt = 0`, `axis_azimuth = 0`).

use crate::solpos::angle_of_incidence;

#[inline]
fn cosd(d: f64) -> f64 {
    d.to_radians().cos()
}
#[inline]
fn sind(d: f64) -> f64 {
    d.to_radians().sin()
}

/// Configuration of a single-axis tracker.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SingleAxisTracker {
    /// Tilt of the rotation axis from horizontal, degrees (0 = horizontal).
    pub axis_tilt: f64,
    /// Azimuth of the rotation axis, degrees clockwise from North
    /// (0 = North–South axis, panels swinging East↔West).
    pub axis_azimuth: f64,
    /// Maximum rotation from horizontal, degrees (the limit is symmetric, ±).
    pub max_angle: f64,
    /// Whether to backtrack to avoid row-to-row shading at low sun angles.
    pub backtrack: bool,
    /// Ground coverage ratio (collector width ÷ row pitch).
    pub gcr: f64,
    /// Cross-axis slope tilt, degrees (0 for level ground).
    pub cross_axis_tilt: f64,
}

impl Default for SingleAxisTracker {
    fn default() -> Self {
        // pvlib defaults: horizontal N–S axis, ±90°, backtracking, gcr 2/7.
        Self {
            axis_tilt: 0.0,
            axis_azimuth: 0.0,
            max_angle: 90.0,
            backtrack: true,
            gcr: 2.0 / 7.0,
            cross_axis_tilt: 0.0,
        }
    }
}

/// Orientation produced by a tracker at a given sun position.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TrackerOrientation {
    /// Tracker rotation from horizontal, degrees (clockwise positive: a
    /// rotation to the West is positive for a South-pointing axis).
    pub tracker_theta: f64,
    /// Resulting surface tilt from horizontal, degrees.
    pub surface_tilt: f64,
    /// Resulting surface azimuth, degrees clockwise from North.
    pub surface_azimuth: f64,
    /// Angle of incidence of the beam on the rotated surface, degrees.
    pub aoi: f64,
}

/// Ideal tracker rotation (no backtracking): the angle that places the sun in
/// the plane normal to the panel and containing the axis. `pvlib`'s
/// `shading.projected_solar_zenith_angle`.
fn projected_solar_zenith(
    axis_tilt: f64,
    axis_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
) -> f64 {
    let sz = sind(solar_zenith);
    let sx = sz * sind(solar_azimuth);
    let sy = sz * cosd(solar_azimuth);
    let szz = cosd(solar_zenith);
    let sx_prime = sx * cosd(axis_azimuth) - sy * sind(axis_azimuth);
    let sz_prime = sx * sind(axis_azimuth) * sind(axis_tilt)
        + sy * sind(axis_tilt) * cosd(axis_azimuth)
        + szz * cosd(axis_tilt);
    sx_prime.atan2(sz_prime).to_degrees()
}

/// Surface tilt & azimuth for a tracker rotated by `theta`. `pvlib`'s
/// `tracking.calc_surface_orientation`.
fn surface_orientation(tracker_theta: f64, axis_tilt: f64, axis_azimuth: f64) -> (f64, f64) {
    let surface_tilt = (cosd(tracker_theta) * cosd(axis_tilt)).clamp(-1.0, 1.0).acos().to_degrees();

    // Unit normal R*(0,0,1) with R = Rz(-azimuth) Rx(-tilt) Ry(theta).
    let (ca, sa) = (cosd(-axis_azimuth), sind(-axis_azimuth));
    let st = sind(-axis_tilt); // only the x/y normal components drive azimuth
    let (cth, sth) = (cosd(tracker_theta), sind(tracker_theta));
    let nx = sa * st * cth + ca * sth;
    let ny = sa * sth - ca * st * cth;

    let mut surface_azimuth = if surface_tilt == 0.0 {
        axis_azimuth - 90.0
    } else {
        nx.atan2(ny).to_degrees()
    };
    surface_azimuth = surface_azimuth.rem_euclid(360.0);
    (surface_tilt, surface_azimuth)
}

/// Compute the tracker orientation for a sun position.
///
/// Returns `None` when the sun is at or below the horizon (`solar_zenith ≥ 90`).
/// `solar_zenith`/`solar_azimuth` are the apparent position in degrees.
pub fn single_axis(
    tracker: &SingleAxisTracker,
    solar_zenith: f64,
    solar_azimuth: f64,
) -> Option<TrackerOrientation> {
    if solar_zenith >= 90.0 {
        return None;
    }

    let omega_ideal =
        projected_solar_zenith(tracker.axis_tilt, tracker.axis_azimuth, solar_zenith, solar_azimuth);

    let mut theta = omega_ideal;
    if tracker.backtrack {
        let axes_distance = 1.0 / (tracker.gcr * cosd(tracker.cross_axis_tilt));
        let temp = (axes_distance * cosd(omega_ideal - tracker.cross_axis_tilt)).abs();
        // arccos only defined for temp < 1; otherwise no shading → no correction.
        if temp < 1.0 {
            let correction = -omega_ideal.signum() * temp.acos().to_degrees();
            theta = omega_ideal + correction;
        }
    }

    theta = theta.clamp(-tracker.max_angle, tracker.max_angle);

    let (surface_tilt, surface_azimuth) =
        surface_orientation(theta, tracker.axis_tilt, tracker.axis_azimuth);
    let aoi = angle_of_incidence(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);

    Some(TrackerOrientation { tracker_theta: theta, surface_tilt, surface_azimuth, aoi })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn night_returns_none() {
        assert!(single_axis(&SingleAxisTracker::default(), 95.0, 180.0).is_none());
    }

    #[test]
    fn morning_sun_tilts_tracker_east() {
        // N–S horizontal axis, sun in the East (azimuth ~90°): panel faces East,
        // so its azimuth is ~90° and it is tilted off horizontal.
        let t = single_axis(&SingleAxisTracker { backtrack: false, ..Default::default() }, 60.0, 90.0)
            .unwrap();
        assert!(t.surface_tilt > 10.0, "tilt {}", t.surface_tilt);
        assert!((t.surface_azimuth - 90.0).abs() < 1.0, "azimuth {}", t.surface_azimuth);
    }

    #[test]
    fn noon_sun_overhead_is_flat() {
        // Sun near zenith → tracker nearly horizontal.
        let t = single_axis(&SingleAxisTracker::default(), 2.0, 180.0).unwrap();
        assert!(t.surface_tilt < 3.0, "tilt {}", t.surface_tilt);
    }

    #[test]
    fn backtracking_reduces_low_sun_rotation() {
        // At a low sun the backtracking tracker rotates less than the ideal one.
        let ideal = single_axis(
            &SingleAxisTracker { backtrack: false, gcr: 0.5, ..Default::default() },
            80.0,
            90.0,
        )
        .unwrap();
        let bt = single_axis(
            &SingleAxisTracker { backtrack: true, gcr: 0.5, ..Default::default() },
            80.0,
            90.0,
        )
        .unwrap();
        assert!(bt.tracker_theta.abs() < ideal.tracker_theta.abs(), "bt {} ideal {}", bt.tracker_theta, ideal.tracker_theta);
    }
}
