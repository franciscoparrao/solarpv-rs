"""Python-side tests for the solarpv PyO3 bindings.

Run with:
    maturin develop --features terrain
    pytest crates/python/tests
"""
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin

import solarpv


def test_version_available():
    assert isinstance(solarpv.__version__, str)


def test_solar_position_noon_southern_hemisphere_winter():
    # Winter solstice, noon UTC at Atacama: sun north of zenith.
    s = solarpv.solar_position(-23.0, -69.0, 2026, 6, 21, 12)
    assert set(s.keys()) == {"zenith", "azimuth", "apparent_zenith", "apparent_elevation", "elevation"}
    assert 0.0 < s["elevation"] < 90.0
    assert 0.0 <= s["azimuth"] <= 360.0


def test_poa_irradiance_perez_components_sum():
    poa = solarpv.poa_irradiance(
        ghi=800.0, dni=600.0, dhi=150.0,
        tilt=23.0, surface_azimuth=0.0,
        solar_zenith=20.0, solar_azimuth=0.0,
        albedo=0.25, dni_extra=1366.0, airmass=1.2,
        model="perez",
    )
    assert set(poa.keys()) == {"global", "direct", "sky_diffuse", "ground_diffuse"}
    assert poa["global"] > 0.0
    assert poa["direct"] > 0.0
    total = poa["direct"] + poa["sky_diffuse"] + poa["ground_diffuse"]
    assert poa["global"] == pytest.approx(total, rel=1e-12)


def test_ac_power_at_stc():
    # 1000 W/m², cool air → output is positive and below nameplate after losses.
    pv = solarpv.ac_power(1000.0, 10.0, 5.0)
    assert set(pv.keys()) == {"tcell", "dc", "ac"}
    assert pv["tcell"] > 10.0
    assert 0.0 < pv["dc"] < 1000.0
    assert pv["ac"] > 0.0


def test_ac_power_hot_cell_reduces_output():
    cold = solarpv.ac_power(1000.0, 10.0, 2.0)
    hot = solarpv.ac_power(1000.0, 40.0, 2.0)
    assert hot["ac"] < cold["ac"]


def _make_flat_dem(path, n=9, cell_deg=0.001, lat=-23.0, lon=-69.0, elev=2400.0):
    data = np.full((n, n), elev, dtype=np.float64)
    transform = from_origin(
        lon - (n / 2) * cell_deg,
        lat + (n / 2) * cell_deg,
        cell_deg,
        cell_deg,
    )
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=n,
        width=n,
        count=1,
        dtype=data.dtype,
        crs="EPSG:4326",
        transform=transform,
    ) as dst:
        dst.write(data, 1)


def test_pv_potential_flat_dem(tmp_path):
    dem = tmp_path / "flat.tif"
    _make_flat_dem(dem)
    res = solarpv.pv_potential(
        str(dem), -23.0, -69.0, 2026, 6, 21, step=30,
    )
    assert set(res.keys()) == {"poa", "ac", "specific_yield", "sun_steps", "shape"}
    rows, cols = res["shape"]
    assert rows == cols == 9
    assert res["sun_steps"] > 0
    center = rows // 2
    assert res["poa"][center][center] > 0.0
    assert res["ac"][center][center] > 0.0
    assert res["specific_yield"][center][center] > 0.0


def test_pv_potential_mounts_differ(tmp_path):
    dem = tmp_path / "flat.tif"
    _make_flat_dem(dem)
    terrain = solarpv.pv_potential(str(dem), -23.0, -69.0, 2026, 6, 21, step=30, mount="terrain")
    tilt = solarpv.pv_potential(str(dem), -23.0, -69.0, 2026, 6, 21, step=30, mount="tilt", tilt=23.0)
    # A flat DEM with terrain mount is horizontal; a 23° tilt changes yield.
    c = 4
    assert terrain["specific_yield"][c][c] != tilt["specific_yield"][c][c]
