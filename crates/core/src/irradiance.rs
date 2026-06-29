//! Irradiance: decomposition of global horizontal irradiance and transposition
//! to the plane of array (POA).
//!
//! Implements the Erbs (1982) diffuse-fraction model and three sky-diffuse
//! transposition models — isotropic (Liu & Jordan), Hay & Davies (1980) and
//! Perez et al. (1990) — matching the corresponding `pvlib.irradiance` routines.
//!
//! Symbols: GHI = global horizontal, DHI = diffuse horizontal, DNI = direct
//! normal. All irradiances in W/m²; all angles in degrees.

use crate::solpos::angle_of_incidence;

const DEG: f64 = std::f64::consts::PI / 180.0;

/// Solar constant (W/m²). Uses 1366.1, the default in
/// `pvlib.irradiance.get_extra_radiation`, so the Spencer expansion below is
/// numerically identical to pvlib (the modern WMO value is 1361, but Spencer's
/// fit is referenced to 1366.1).
pub const SOLAR_CONSTANT: f64 = 1366.1;

/// Extraterrestrial normal irradiance (W/m²) for a day of year.
///
/// Spencer (1971) Fourier expansion of the Earth–Sun distance correction, as in
/// `pvlib.irradiance.get_extra_radiation(method="spencer")`.
pub fn extra_radiation(day_of_year: u32) -> f64 {
    let b = 2.0 * std::f64::consts::PI * (day_of_year as f64 - 1.0) / 365.0;
    let r0 = 1.000_110
        + 0.034_221 * b.cos()
        + 0.001_280 * b.sin()
        + 0.000_719 * (2.0 * b).cos()
        + 0.000_077 * (2.0 * b).sin();
    SOLAR_CONSTANT * r0
}

/// Kasten & Young (1989) relative (unitless) air mass from apparent zenith.
///
/// Returns `f64::INFINITY` for the sun at or below the horizon. Mirrors
/// `pvlib.atmosphere.get_relative_airmass(model="kastenyoung1989")`.
pub fn relative_airmass(apparent_zenith_deg: f64) -> f64 {
    if apparent_zenith_deg >= 90.0 {
        return f64::INFINITY;
    }
    let z = apparent_zenith_deg;
    1.0 / ((z * DEG).cos() + 0.505_72 * (96.079_95 - z).powf(-1.636_4))
}

/// Haurwitz (1945) clear-sky global horizontal irradiance (W/m²).
///
/// A single-parameter clear-sky model depending only on the apparent solar
/// zenith — useful for a self-contained "clear-sky PV potential" map when no
/// measured or TMY irradiance series is available. Mirrors
/// `pvlib.clearsky.haurwitz`. Returns 0 with the sun at or below the horizon.
pub fn haurwitz_clearsky_ghi(apparent_zenith_deg: f64) -> f64 {
    if apparent_zenith_deg >= 90.0 {
        return 0.0;
    }
    let cos_z = (apparent_zenith_deg * DEG).cos();
    (1098.0 * cos_z * (-0.059 / cos_z).exp()).max(0.0)
}

/// Global irradiance split into its direct-normal and diffuse-horizontal parts.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Decomposition {
    /// Global horizontal irradiance (input), W/m².
    pub ghi: f64,
    /// Direct normal irradiance, W/m².
    pub dni: f64,
    /// Diffuse horizontal irradiance, W/m².
    pub dhi: f64,
}

/// Erbs (1982) decomposition: estimate DNI and DHI from GHI.
///
/// `zenith_deg` is the apparent solar zenith; `day_of_year` drives the
/// extraterrestrial reference. Matches `pvlib.irradiance.erbs`.
pub fn erbs(ghi: f64, zenith_deg: f64, day_of_year: u32) -> Decomposition {
    let cos_z = (zenith_deg * DEG).cos();
    // Below the horizon (or numerically near it) there is no usable beam.
    if zenith_deg >= 90.0 || ghi <= 0.0 || cos_z <= 0.0 {
        return Decomposition { ghi: ghi.max(0.0), dni: 0.0, dhi: ghi.max(0.0) };
    }
    let i0h = extra_radiation(day_of_year) * cos_z;
    let kt = (ghi / i0h).clamp(0.0, 1.0);
    let df = if kt <= 0.22 {
        1.0 - 0.09 * kt
    } else if kt <= 0.80 {
        0.9511 - 0.1604 * kt + 4.388 * kt.powi(2) - 16.638 * kt.powi(3) + 12.336 * kt.powi(4)
    } else {
        0.165
    };
    let dhi = ghi * df;
    let dni = (ghi - dhi) / cos_z;
    Decomposition { ghi, dni, dhi }
}

/// Sky-diffuse transposition model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SkyModel {
    /// Liu & Jordan isotropic sky.
    Isotropic,
    /// Hay & Davies (1980) anisotropic (circumsolar) model.
    HayDavies,
    /// Perez et al. (1990) anisotropic model (`allsitescomposite1990`).
    Perez,
}

/// Plane-of-array irradiance components (W/m²).
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PoaIrradiance {
    /// Total POA irradiance reaching the module plane.
    pub global: f64,
    /// Beam (direct) component on the plane.
    pub direct: f64,
    /// Sky-diffuse component on the plane.
    pub sky_diffuse: f64,
    /// Ground-reflected component on the plane.
    pub ground_diffuse: f64,
}

/// Ground-reflected diffuse on a tilted plane (isotropic ground).
fn ground_diffuse(ghi: f64, surface_tilt: f64, albedo: f64) -> f64 {
    ghi * albedo * (1.0 - (surface_tilt * DEG).cos()) / 2.0
}

/// Perez (1990) F1/F2 coefficient table, set `allsitescomposite1990`
/// (the `pvlib` default). Rows are the 8 clearness bins; columns
/// `[F11, F12, F13, F21, F22, F23]`.
const PEREZ_COEFFS: [[f64; 6]; 8] = [
    [-0.008, 0.588, -0.062, -0.060, 0.072, -0.022],
    [0.130, 0.683, -0.151, -0.019, 0.066, -0.029],
    [0.330, 0.487, -0.221, 0.055, -0.064, -0.026],
    [0.568, 0.187, -0.295, 0.109, -0.152, -0.014],
    [0.873, -0.392, -0.362, 0.226, -0.462, 0.001],
    [1.132, -1.237, -0.412, 0.288, -0.823, 0.056],
    [1.060, -1.600, -0.359, 0.264, -1.127, 0.131],
    [0.678, -0.327, -0.250, 0.156, -1.377, 0.251],
];

/// Select the Perez clearness bin (0-based) from the clearness index ε.
fn perez_bin(eps: f64) -> usize {
    match eps {
        e if e < 1.065 => 0,
        e if e < 1.23 => 1,
        e if e < 1.5 => 2,
        e if e < 1.95 => 3,
        e if e < 2.8 => 4,
        e if e < 4.5 => 5,
        e if e < 6.2 => 6,
        _ => 7,
    }
}

/// Perez (1990) sky-diffuse irradiance on a tilted plane (W/m²).
///
/// `air_mass` is the relative air mass (see [`relative_airmass`]); `dni_extra`
/// is the extraterrestrial normal irradiance (see [`extra_radiation`]).
/// Mirrors `pvlib.irradiance.perez`.
#[allow(clippy::too_many_arguments)]
pub fn perez_sky_diffuse(
    surface_tilt: f64,
    surface_azimuth: f64,
    dhi: f64,
    dni: f64,
    dni_extra: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    air_mass: f64,
) -> f64 {
    if dhi <= 0.0 {
        return 0.0;
    }
    let z = solar_zenith * DEG;
    let kappa = 1.041;
    let kz3 = kappa * z.powi(3);
    let eps = ((dhi + dni) / dhi + kz3) / (1.0 + kz3);
    let c = PEREZ_COEFFS[perez_bin(eps)];

    // Sky brightness.
    let delta = dhi * air_mass / dni_extra;
    let f1 = (c[0] + c[1] * delta + c[2] * z).max(0.0);
    let f2 = c[3] + c[4] * delta + c[5] * z;

    let aoi = angle_of_incidence(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);
    let a = (aoi * DEG).cos().max(0.0);
    let b = (solar_zenith * DEG).cos().max((85.0_f64 * DEG).cos());

    dhi * ((1.0 - f1) * (1.0 + (surface_tilt * DEG).cos()) / 2.0
        + f1 * a / b
        + f2 * (surface_tilt * DEG).sin())
}

/// Hay & Davies (1980) sky-diffuse irradiance on a tilted plane (W/m²).
fn hay_davies_sky_diffuse(
    surface_tilt: f64,
    surface_azimuth: f64,
    dhi: f64,
    dni: f64,
    dni_extra: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
) -> f64 {
    if dhi <= 0.0 {
        return 0.0;
    }
    // Anisotropy index and geometric (beam) ratio Rb = cos(aoi)/cos(zenith).
    let ai = dni / dni_extra;
    let aoi = angle_of_incidence(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);
    let cos_z = (solar_zenith * DEG).cos();
    let rb = if cos_z > 0.0 {
        ((aoi * DEG).cos() / cos_z).max(0.0)
    } else {
        0.0
    };
    dhi * (ai * rb + (1.0 - ai) * (1.0 + (surface_tilt * DEG).cos()) / 2.0)
}

/// Compute plane-of-array irradiance from decomposed components.
///
/// `albedo` is the ground reflectance (e.g. 0.2 vegetation, 0.7 snow). For the
/// Perez model, `air_mass` and `dni_extra` must be supplied; pass the values
/// from [`relative_airmass`] and [`extra_radiation`].
#[allow(clippy::too_many_arguments)]
pub fn poa_irradiance(
    model: SkyModel,
    surface_tilt: f64,
    surface_azimuth: f64,
    decomp: Decomposition,
    solar_zenith: f64,
    solar_azimuth: f64,
    albedo: f64,
    dni_extra: f64,
    air_mass: f64,
) -> PoaIrradiance {
    let aoi = angle_of_incidence(surface_tilt, surface_azimuth, solar_zenith, solar_azimuth);
    let direct = (decomp.dni * (aoi * DEG).cos()).max(0.0);

    let sky_diffuse = match model {
        SkyModel::Isotropic => {
            decomp.dhi * (1.0 + (surface_tilt * DEG).cos()) / 2.0
        }
        SkyModel::HayDavies => hay_davies_sky_diffuse(
            surface_tilt,
            surface_azimuth,
            decomp.dhi,
            decomp.dni,
            dni_extra,
            solar_zenith,
            solar_azimuth,
        ),
        SkyModel::Perez => perez_sky_diffuse(
            surface_tilt,
            surface_azimuth,
            decomp.dhi,
            decomp.dni,
            dni_extra,
            solar_zenith,
            solar_azimuth,
            air_mass,
        ),
    };

    let ground = ground_diffuse(decomp.ghi, surface_tilt, albedo);
    PoaIrradiance {
        global: direct + sky_diffuse + ground,
        direct,
        sky_diffuse,
        ground_diffuse: ground,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extra_radiation_in_range() {
        // Annual variation stays within ±3.4% of the solar constant.
        for doy in [1, 80, 172, 264, 355] {
            let e = extra_radiation(doy);
            assert!((1310.0..=1415.0).contains(&e), "doy {doy}: {e}");
        }
    }

    #[test]
    fn erbs_clear_sky_is_mostly_beam() {
        // High GHI near solar noon → low diffuse fraction, large DNI.
        let d = erbs(900.0, 25.0, 172);
        assert!(d.dni > d.dhi, "dni {} dhi {}", d.dni, d.dhi);
        assert!((d.dni * (25.0_f64.to_radians()).cos() + d.dhi - 900.0).abs() < 1e-6);
    }

    #[test]
    fn isotropic_poa_horizontal_equals_ghi() {
        // On a horizontal plane, POA must reconstruct GHI = DNI·cos(z) + DHI.
        let zenith = 30.0_f64;
        let dni = 700.0;
        let dhi = 140.0;
        let ghi = dni * zenith.to_radians().cos() + dhi;
        let decomp = Decomposition { ghi, dni, dhi };
        let poa = poa_irradiance(
            SkyModel::Isotropic, 0.0, 180.0, decomp, zenith, 150.0, 0.2,
            extra_radiation(172), relative_airmass(zenith),
        );
        assert!((poa.global - ghi).abs() < 1e-6, "global {} ghi {}", poa.global, ghi);
        assert!(poa.ground_diffuse.abs() < 1e-9);
    }
}
