//! Detailed PV losses: incidence-angle modifiers (IAM) and the PVWatts loss
//! breakdown.
//!
//! The IAM models give the angular reflection loss on the module cover as a
//! function of the angle of incidence (AOI); they mirror `pvlib.iam`. The
//! [`PvLosses`] breakdown replaces a single flat derate with named components,
//! combined as in `pvlib.pvsystem.pvwatts_losses`.

#[inline]
fn cosd(d: f64) -> f64 {
    d.to_radians().cos()
}

/// Incidence-angle-modifier model. `iam(aoi)` returns the transmittance relative
/// to normal incidence (1.0 at AOI = 0, falling towards 0 near grazing).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum IamModel {
    /// ASHRAE model: `1 − b·(1/cos(aoi) − 1)`.
    Ashrae { b: f64 },
    /// Martín & Ruiz exponential model with angular-loss coefficient `a_r`.
    MartinRuiz { a_r: f64 },
    /// Physical (Fresnel + glass absorption) model: refractive index `n`,
    /// extinction coefficient `k` (1/m) and cover thickness `l` (m).
    Physical { n: f64, k: f64, l: f64 },
}

impl IamModel {
    /// ASHRAE with the common `b = 0.05`.
    pub const fn ashrae() -> Self {
        IamModel::Ashrae { b: 0.05 }
    }
    /// Martín–Ruiz with the common `a_r = 0.16`.
    pub const fn martin_ruiz() -> Self {
        IamModel::MartinRuiz { a_r: 0.16 }
    }
    /// Physical model with pvlib's defaults (`n = 1.526`, `k = 4`, `l = 0.002`).
    pub const fn physical() -> Self {
        IamModel::Physical { n: 1.526, k: 4.0, l: 0.002 }
    }

    /// Incidence-angle modifier at angle `aoi` (degrees).
    pub fn iam(&self, aoi: f64) -> f64 {
        match *self {
            IamModel::Ashrae { b } => {
                if aoi.abs() >= 90.0 {
                    0.0
                } else {
                    (1.0 - b * (1.0 / cosd(aoi) - 1.0)).max(0.0)
                }
            }
            IamModel::MartinRuiz { a_r } => {
                ((1.0 - (-cosd(aoi) / a_r).exp()) / (1.0 - (-1.0 / a_r).exp())).max(0.0)
            }
            IamModel::Physical { n, k, l } => physical_iam(aoi, n, k, l),
        }
    }
}

/// Physical IAM (no anti-reflective coating, `n_ar = n`), per `pvlib.iam.physical`.
fn physical_iam(aoi: f64, n: f64, k: f64, l: f64) -> f64 {
    let (n1, n2) = (1.0_f64, n);
    let costheta1 = cosd(aoi).max(0.0);
    let sintheta1 = (1.0 - costheta1 * costheta1).sqrt();
    // Refraction at the air→glass interface.
    let sintheta2 = n1 / n2 * sintheta1;
    let costheta2 = (1.0 - sintheta2 * sintheta2).sqrt();

    let (n1c1, n2c1) = (n1 * costheta1, n2 * costheta1);
    let (n1c2, n2c2) = (n1 * costheta2, n2 * costheta2);
    let rho_s = ((n1c1 - n2c2) / (n1c1 + n2c2)).powi(2);
    let rho_p = ((n1c2 - n2c1) / (n1c2 + n2c1)).powi(2);
    let rho_0 = ((n1 - n2) / (n1 + n2)).powi(2);

    // Transmittance through the interface and glass absorption.
    let tau_s = (1.0 - rho_s) * (-k * l / costheta2).exp();
    let tau_p = (1.0 - rho_p) * (-k * l / costheta2).exp();
    let tau_0 = (1.0 - rho_0) * (-k * l).exp();

    let iam = (tau_s + tau_p) / 2.0 / tau_0;
    if aoi >= 90.0 {
        0.0
    } else {
        iam.max(0.0)
    }
}

/// PVWatts system loss breakdown (each value is a percentage, 0–100).
///
/// The components combine multiplicatively, exactly as
/// `pvlib.pvsystem.pvwatts_losses`:
/// `total = 1 − Π(1 − lᵢ/100)`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PvLosses {
    pub soiling: f64,
    pub shading: f64,
    pub snow: f64,
    pub mismatch: f64,
    pub wiring: f64,
    pub connections: f64,
    pub lid: f64,
    pub nameplate_rating: f64,
    pub age: f64,
    pub availability: f64,
}

impl Default for PvLosses {
    fn default() -> Self {
        // pvlib.pvsystem.pvwatts_losses defaults → ~14.08 % total.
        Self {
            soiling: 2.0,
            shading: 3.0,
            snow: 0.0,
            mismatch: 2.0,
            wiring: 2.0,
            connections: 0.5,
            lid: 1.5,
            nameplate_rating: 1.0,
            age: 0.0,
            availability: 3.0,
        }
    }
}

impl PvLosses {
    /// Combined loss as a fraction in `[0, 1]` (e.g. 0.1408 for the defaults).
    pub fn total_fraction(&self) -> f64 {
        let comps = [
            self.soiling,
            self.shading,
            self.snow,
            self.mismatch,
            self.wiring,
            self.connections,
            self.lid,
            self.nameplate_rating,
            self.age,
            self.availability,
        ];
        1.0 - comps.iter().map(|l| 1.0 - l / 100.0).product::<f64>()
    }

    /// Combined loss as a percentage (matches `pvlib.pvsystem.pvwatts_losses`).
    pub fn total_percent(&self) -> f64 {
        self.total_fraction() * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iam_unity_at_normal_incidence() {
        for m in [IamModel::ashrae(), IamModel::martin_ruiz(), IamModel::physical()] {
            assert!((m.iam(0.0) - 1.0).abs() < 0.02, "{m:?} at 0° = {}", m.iam(0.0));
        }
    }

    #[test]
    fn iam_decreases_with_angle_and_zero_past_grazing() {
        for m in [IamModel::ashrae(), IamModel::martin_ruiz(), IamModel::physical()] {
            assert!(m.iam(60.0) < m.iam(20.0), "{m:?} not decreasing");
            assert!(m.iam(95.0) <= 1e-9, "{m:?} should vanish past 90°");
        }
    }

    #[test]
    fn pvwatts_losses_default_is_about_14_percent() {
        assert!((PvLosses::default().total_percent() - 14.08).abs() < 0.05);
    }
}
