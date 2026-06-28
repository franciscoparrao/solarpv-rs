//! PV conversion: cell temperature, DC power and AC energy yield.
//!
//! Uses the PVWatts (Dobos, 2014) DC and inverter models together with the SAPM
//! cell-temperature model, matching `pvlib.temperature.sapm_cell`,
//! `pvlib.pvsystem.pvwatts_dc` and `pvlib.inverter.pvwatts`.

/// Parameters of the SAPM module/cell temperature model.
///
/// `Tcell = POA·exp(a + b·wind) + Tair + POA/1000·ΔT`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TempModel {
    /// Empirical coefficient `a` (upper temperature limit at low wind).
    pub a: f64,
    /// Empirical coefficient `b` (wind-speed dependence).
    pub b: f64,
    /// Temperature difference between cell and module back at 1000 W/m² (°C).
    pub delta_t: f64,
}

impl TempModel {
    /// `open_rack_glass_glass` — the `pvlib` default mounting preset.
    pub const OPEN_RACK_GLASS_GLASS: TempModel = TempModel { a: -3.47, b: -0.0594, delta_t: 3.0 };
    /// `close_mount_glass_glass` — roof-parallel close mounting.
    pub const CLOSE_MOUNT_GLASS_GLASS: TempModel = TempModel { a: -2.98, b: -0.0471, delta_t: 1.0 };
    /// `open_rack_glass_polymer` — common utility-scale preset.
    pub const OPEN_RACK_GLASS_POLYMER: TempModel = TempModel { a: -3.56, b: -0.0750, delta_t: 3.0 };
}

impl Default for TempModel {
    fn default() -> Self {
        TempModel::OPEN_RACK_GLASS_GLASS
    }
}

/// SAPM cell temperature (°C).
///
/// `poa_global` in W/m², `temp_air` in °C, `wind_speed` in m/s (at 10 m).
/// Mirrors `pvlib.temperature.sapm_cell`.
pub fn sapm_cell_temperature(
    poa_global: f64,
    temp_air: f64,
    wind_speed: f64,
    model: TempModel,
) -> f64 {
    let module_temp = poa_global * (model.a + model.b * wind_speed).exp() + temp_air;
    module_temp + (poa_global / 1000.0) * model.delta_t
}

/// PVWatts DC power (W).
///
/// `pdc0` is the nameplate DC power at STC (W); `gamma_pdc` the power
/// temperature coefficient (fraction per °C, e.g. −0.004). Mirrors
/// `pvlib.pvsystem.pvwatts_dc` with `temp_ref = 25 °C`.
pub fn pvwatts_dc(poa_effective: f64, temp_cell: f64, pdc0: f64, gamma_pdc: f64) -> f64 {
    if poa_effective <= 0.0 {
        return 0.0;
    }
    (poa_effective / 1000.0) * pdc0 * (1.0 + gamma_pdc * (temp_cell - 25.0))
}

/// PVWatts inverter model: DC → AC power (W).
///
/// `pdc0_inv` is the inverter's rated DC input (W); `eta_nom`/`eta_ref` are the
/// nominal and reference efficiencies (defaults 0.96 / 0.9637). Mirrors
/// `pvlib.inverter.pvwatts`.
pub fn pvwatts_ac(pdc: f64, pdc0_inv: f64, eta_nom: f64, eta_ref: f64) -> f64 {
    if pdc <= 0.0 {
        return 0.0;
    }
    let pac0 = eta_nom * pdc0_inv;
    let zeta = pdc / pdc0_inv;
    let eta = (eta_nom / eta_ref) * (-0.0162 * zeta - 0.0059 / zeta + 0.9858);
    (eta * pdc).clamp(0.0, pac0)
}

/// Defaults for [`pvwatts_ac`] (`eta_nom = 0.96`, `eta_ref = 0.9637`).
pub const ETA_INV_NOM: f64 = 0.96;
/// Reference inverter efficiency used to normalise the PVWatts curve.
pub const ETA_INV_REF: f64 = 0.9637;

/// Configuration of a simple grid-tied PV system at a point.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PvSystem {
    /// Nameplate DC power at STC (W).
    pub pdc0: f64,
    /// Power temperature coefficient (fraction per °C), typically negative.
    pub gamma_pdc: f64,
    /// Inverter rated DC input (W).
    pub pdc0_inv: f64,
    /// Fractional system losses applied to DC (soiling, wiring, mismatch, …).
    /// PVWatts default is 0.14.
    pub system_losses: f64,
    /// Cell temperature model.
    pub temp_model: TempModel,
}

impl PvSystem {
    /// A 1 kW reference system with PVWatts defaults: −0.4 %/°C, DC/AC ≈ 1.1,
    /// 14 % system losses, open-rack glass/glass.
    pub fn reference_1kw() -> Self {
        Self {
            pdc0: 1000.0,
            gamma_pdc: -0.004,
            pdc0_inv: 1000.0 / 1.1,
            system_losses: 0.14,
            temp_model: TempModel::OPEN_RACK_GLASS_GLASS,
        }
    }
}

/// AC power (W) from POA irradiance and weather for a [`PvSystem`].
///
/// Applies the SAPM cell temperature, PVWatts DC with the fractional system
/// losses, then the PVWatts inverter.
pub fn ac_power(
    system: &PvSystem,
    poa_global: f64,
    temp_air: f64,
    wind_speed: f64,
) -> f64 {
    let t_cell = sapm_cell_temperature(poa_global, temp_air, wind_speed, system.temp_model);
    let dc = pvwatts_dc(poa_global, t_cell, system.pdc0, system.gamma_pdc)
        * (1.0 - system.system_losses);
    pvwatts_ac(dc, system.pdc0_inv, ETA_INV_NOM, ETA_INV_REF)
}

/// Integrate a series of AC power samples (W) into energy (Wh).
///
/// `time_step_hours` is the spacing between consecutive samples (e.g. 1.0 for
/// hourly). Uses a left-Riemann sum, the convention for instantaneous-rate
/// weather series.
pub fn integrate_energy_wh(power_w: &[f64], time_step_hours: f64) -> f64 {
    power_w.iter().sum::<f64>() * time_step_hours
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_hotter_than_air_under_sun() {
        let t = sapm_cell_temperature(800.0, 20.0, 1.0, TempModel::OPEN_RACK_GLASS_GLASS);
        assert!(t > 20.0 && t < 70.0, "tcell {t}");
    }

    #[test]
    fn dc_at_stc_equals_nameplate() {
        // 1000 W/m², 25 °C cell → exactly nameplate power.
        let dc = pvwatts_dc(1000.0, 25.0, 1000.0, -0.004);
        assert!((dc - 1000.0).abs() < 1e-9, "dc {dc}");
    }

    #[test]
    fn hot_cell_reduces_dc() {
        let cold = pvwatts_dc(1000.0, 25.0, 1000.0, -0.004);
        let hot = pvwatts_dc(1000.0, 50.0, 1000.0, -0.004);
        assert!(hot < cold);
        // −0.4 %/°C over 25 °C → −10 %.
        assert!((hot - 900.0).abs() < 1e-6, "hot {hot}");
    }

    #[test]
    fn inverter_clips_at_rated() {
        let pac = pvwatts_ac(2000.0, 1000.0, ETA_INV_NOM, ETA_INV_REF);
        assert!((pac - ETA_INV_NOM * 1000.0).abs() < 1e-9, "pac {pac}");
    }

    #[test]
    fn energy_integration_hourly() {
        let p = [0.0, 100.0, 200.0, 100.0, 0.0];
        assert!((integrate_energy_wh(&p, 1.0) - 400.0).abs() < 1e-9);
    }
}
