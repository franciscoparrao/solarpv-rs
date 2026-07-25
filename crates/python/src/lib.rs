//! Python bindings for `solarpv-core` via PyO3.
//!
//! Exposes the validated point model chain (solar position, POA irradiance,
//! PV power) and, with the `terrain` feature, the gridded engine over a DEM.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use solarpv_core::irradiance::{poa_irradiance, Decomposition, SkyModel};
use solarpv_core::pv::{ac_power, pvwatts_ac, pvwatts_dc, sapm_cell_temperature, PvSystem, TempModel};
use solarpv_core::solpos::{solar_position, DateTimeUtc, Location};

/// Solar position at a given instant.
///
/// Returns a dict with keys `zenith`, `azimuth`, `apparent_zenith`,
/// `apparent_elevation`, `elevation` (all in degrees).
#[pyfunction(name = "solar_position")]
#[pyo3(signature = (lat, lon, year, month, day, hour, minute=0, second=0))]
fn solar_position_py(
    lat: f64,
    lon: f64,
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: u32,
) -> PyResult<PyObject> {
    Python::with_gil(|py| {
        let when = DateTimeUtc::new(year, month, day, hour, minute, second)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let loc = Location::new(lat, lon)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let s = solar_position(when, loc);

        let dict = PyDict::new(py);
        dict.set_item("zenith", s.zenith)?;
        dict.set_item("azimuth", s.azimuth)?;
        dict.set_item("apparent_zenith", s.apparent_zenith)?;
        dict.set_item("apparent_elevation", s.apparent_elevation)?;
        dict.set_item("elevation", s.elevation)?;
        Ok(dict.into())
    })
}

/// Plane-of-array irradiance components (W/m²).
///
/// `model` may be `"isotropic"`, `"haydavies"` or `"perez"`.
#[pyfunction(name = "poa_irradiance")]
#[pyo3(signature = (
    ghi, dni, dhi,
    tilt, surface_azimuth,
    solar_zenith, solar_azimuth,
    albedo, dni_extra, airmass,
    model="perez"
))]
fn poa_irradiance_py(
    ghi: f64,
    dni: f64,
    dhi: f64,
    tilt: f64,
    surface_azimuth: f64,
    solar_zenith: f64,
    solar_azimuth: f64,
    albedo: f64,
    dni_extra: f64,
    airmass: f64,
    model: &str,
) -> PyResult<PyObject> {
    Python::with_gil(|py| {
        let sky_model = parse_sky_model(model)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
        let decomp = Decomposition { ghi, dni, dhi };
        let poa = poa_irradiance(
            sky_model, tilt, surface_azimuth, decomp,
            solar_zenith, solar_azimuth, albedo, dni_extra, airmass,
        );

        let dict = PyDict::new(py);
        dict.set_item("global", poa.global)?;
        dict.set_item("direct", poa.direct)?;
        dict.set_item("sky_diffuse", poa.sky_diffuse)?;
        dict.set_item("ground_diffuse", poa.ground_diffuse)?;
        Ok(dict.into())
    })
}

/// Run the PV chain for a single time step.
///
/// Returns a dict with `tcell` (°C), `dc` (W) and `ac` (W).
#[pyfunction(name = "ac_power")]
#[pyo3(signature = (
    poa_global, temp_air, wind,
    pdc0=1000.0, gamma_pdc=-0.004, system_losses=0.14,
    pdc0_inv=None, eta_nom=None, eta_ref=None,
    temp_model="open_rack_glass_glass"
))]
fn ac_power_py(
    poa_global: f64,
    temp_air: f64,
    wind: f64,
    pdc0: f64,
    gamma_pdc: f64,
    system_losses: f64,
    pdc0_inv: Option<f64>,
    eta_nom: Option<f64>,
    eta_ref: Option<f64>,
    temp_model: &str,
) -> PyResult<PyObject> {
    Python::with_gil(|py| {
        let temp_model = parse_temp_model(temp_model)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
        let system = PvSystem {
            pdc0,
            gamma_pdc,
            pdc0_inv: pdc0_inv.unwrap_or(pdc0 / 1.1),
            system_losses,
            temp_model,
        };

        let tcell = sapm_cell_temperature(poa_global, temp_air, wind, temp_model);
        let dc = pvwatts_dc(poa_global, tcell, system.pdc0, system.gamma_pdc)
            * (1.0 - system.system_losses);
        let ac = if let (Some(eta_nom), Some(eta_ref)) = (eta_nom, eta_ref) {
            pvwatts_ac(dc, system.pdc0_inv, eta_nom, eta_ref)
        } else {
            ac_power(&system, poa_global, temp_air, wind)
        };

        let dict = PyDict::new(py);
        dict.set_item("tcell", tcell)?;
        dict.set_item("dc", dc)?;
        dict.set_item("ac", ac)?;
        Ok(dict.into())
    })
}

fn parse_sky_model(s: &str) -> Result<SkyModel, String> {
    match s.to_lowercase().as_str() {
        "isotropic" => Ok(SkyModel::Isotropic),
        "haydavies" | "hay-davies" => Ok(SkyModel::HayDavies),
        "perez" => Ok(SkyModel::Perez),
        _ => Err(format!("unknown sky model `{s}`; expected isotropic/haydavies/perez")),
    }
}

fn parse_temp_model(s: &str) -> Result<TempModel, String> {
    match s.to_lowercase().as_str() {
        "open_rack_glass_glass" => Ok(TempModel::OPEN_RACK_GLASS_GLASS),
        "close_mount_glass_glass" => Ok(TempModel::CLOSE_MOUNT_GLASS_GLASS),
        "open_rack_glass_polymer" => Ok(TempModel::OPEN_RACK_GLASS_POLYMER),
        _ => Err(format!("unknown temperature model `{s}`")),
    }
}

#[cfg(feature = "terrain")]
mod terrain {
    use super::*;
    use solarpv_core::grid::{pv_potential, GridConfig, LatitudeMode, Mount};
    use solarpv_core::pv::PvSystem;
    use solarpv_core::solpos::{DateTimeUtc, Location};
    use solarpv_core::tracking::{DualAxisTracker, SingleAxisTracker};
    use surtgis_core::io::read_geotiff;

    fn raster_to_vec2d(r: &surtgis_core::Raster<f64>) -> Vec<Vec<f64>> {
        let (rows, cols) = r.shape();
        (0..rows)
            .map(|i| (0..cols).map(|j| r.get(i, j).unwrap_or(f64::NAN)).collect())
            .collect()
    }

    /// Gridded PV potential over a DEM.
    ///
    /// Returns a dict with `poa` (Wh/m²), `ac` (Wh), `specific_yield` (kWh/kWp),
    /// `sun_steps` and `shape` (rows, cols).
    #[pyfunction(name = "pv_potential")]
    #[pyo3(signature = (
        dem_path,
        lat, lon,
        year, month, day,
        step=60,
        sky="perez",
        albedo=0.25,
        pdc0=1000.0,
        gamma_pdc=-0.004,
        system_losses=0.14,
        temp_air=18.0,
        wind=2.0,
        pressure_pa=101325.0,
        horizon_radius=100,
        horizon_dirs=36,
        per_cell_lat=false,
        mount="terrain",
        tilt=0.0,
        surface_azimuth=0.0,
        gcr=2.0/7.0,
        max_tilt=90.0,
        iam=None,
        spectral=false,
    ))]
    #[allow(clippy::too_many_arguments)]
    pub fn pv_potential_py(
        dem_path: &str,
        lat: f64,
        lon: f64,
        year: i32,
        month: u32,
        day: u32,
        step: u32,
        sky: &str,
        albedo: f64,
        pdc0: f64,
        gamma_pdc: f64,
        system_losses: f64,
        temp_air: f64,
        wind: f64,
        pressure_pa: f64,
        horizon_radius: usize,
        horizon_dirs: usize,
        per_cell_lat: bool,
        mount: &str,
        tilt: f64,
        surface_azimuth: f64,
        gcr: f64,
        max_tilt: f64,
        iam: Option<&str>,
        spectral: bool,
    ) -> PyResult<PyObject> {
        Python::with_gil(|py| {
            let dem = read_geotiff::<f64, _>(dem_path, None)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

            let date = DateTimeUtc::new(year, month, day, 0, 0, 0)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
            let center = Location::new(lat, lon)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

            let sky_model = parse_sky_model(sky)
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
            let mut cfg = GridConfig::new(center, date);
            cfg.time_step_minutes = step;
            cfg.sky_model = sky_model;
            cfg.albedo = albedo;
            cfg.system = PvSystem {
                pdc0,
                gamma_pdc,
                pdc0_inv: pdc0 / 1.1,
                system_losses,
                temp_model: TempModel::OPEN_RACK_GLASS_GLASS,
            };
            cfg.temp_air = temp_air;
            cfg.wind = wind;
            cfg.pressure_pa = pressure_pa;
            cfg.horizon.radius = horizon_radius;
            cfg.horizon.directions = horizon_dirs;
            cfg.latitude_mode = if per_cell_lat {
                LatitudeMode::PerCellGeographic
            } else {
                LatitudeMode::Center
            };
            cfg.mount = match mount {
                "terrain" => Mount::FixedTerrain,
                "tilt" => Mount::FixedTilt { tilt, surface_azimuth },
                "tracker" => Mount::SingleAxis(SingleAxisTracker { gcr, ..Default::default() }),
                "dual-axis" | "dual_axis" => Mount::DualAxis(DualAxisTracker { max_tilt }),
                _ => return Err(pyo3::exceptions::PyValueError::new_err(
                    format!("unknown mount `{mount}`; expected terrain/tilt/tracker/dual-axis"))),
            };
            cfg.iam = iam.map(|kind| match kind {
                "ashrae" => solarpv_core::losses::IamModel::ashrae(),
                "martin_ruiz" | "martin-ruiz" => solarpv_core::losses::IamModel::martin_ruiz(),
                "physical" => solarpv_core::losses::IamModel::physical(),
                _ => panic!("unknown IAM model `{kind}`"),
            });
            cfg.spectral = if spectral {
                Some(solarpv_core::losses::SpectralLoss::c_si())
            } else {
                None
            };

            let res = pv_potential(&dem, &cfg)
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

            let dict = PyDict::new(py);
            dict.set_item("poa", raster_to_vec2d(&res.poa_wh))?;
            dict.set_item("ac", raster_to_vec2d(&res.ac_wh))?;
            dict.set_item("specific_yield", raster_to_vec2d(&res.specific_yield))?;
            dict.set_item("sun_steps", res.sun_steps)?;
            dict.set_item("shape", (res.poa_wh.shape().0, res.poa_wh.shape().1))?;
            Ok(dict.into())
        })
    }
}

/// The solarpv Python module.
#[pymodule]
fn solarpv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(solar_position_py, m)?)?;
    m.add_function(wrap_pyfunction!(poa_irradiance_py, m)?)?;
    m.add_function(wrap_pyfunction!(ac_power_py, m)?)?;
    #[cfg(feature = "terrain")]
    m.add_function(wrap_pyfunction!(terrain::pv_potential_py, m)?)?;
    Ok(())
}
