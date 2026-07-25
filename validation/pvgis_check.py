#!/usr/bin/env python3
"""Cross-check solarpv-rs against the PVGIS v5.2 API.

Queries PVGIS for the Atacama reference point, writes the hourly ERA5 series to
a solarpv-cli weather CSV, runs the Rust engine on a small flat DEM, and compares
annual energy yield.

Run inside the project venv:
    python3 validation/pvgis_check.py
"""
import csv
import json
import pathlib
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pvlib
import rasterio
from rasterio.transform import from_origin
import requests

LAT, LON = -23.0, -69.0
YEAR = 2020
SYSTEM_LOSS = 14.0
TILT, AZIMUTH = 23.0, 0.0
PDC0 = 1000.0
GAMMA_PDC = -0.004

PVGIS_URL = "https://re.jrc.ec.europa.eu/api/v5_2/seriescalc"


def pvgis_hourly(pvcalculation: int):
    params = {
        "lat": LAT,
        "lon": LON,
        "startyear": YEAR,
        "endyear": YEAR,
        "pvcalculation": pvcalculation,
        "peakpower": 1.0,
        "loss": SYSTEM_LOSS,
        "mountingplace": "free",
        "angle": TILT,
        "aspect": AZIMUTH,
        "outputformat": "json",
    }
    if pvcalculation == 0:
        # Drop PV-specific params for the irradiance-only query.
        params.pop("peakpower", None)
        params.pop("loss", None)
        params.pop("mountingplace", None)
        params.pop("angle", None)
        params.pop("aspect", None)

    r = requests.get(PVGIS_URL, params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def parse_time(t: str):
    """PVGIS time format 'YYYYMMDD:HHMM' -> (year, month, day, hour)."""
    return int(t[:4]), int(t[4:6]), int(t[6:8]), int(t[9:11])


def write_weather_csv(path: pathlib.Path, records):
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["year", "month", "day", "hour", "ghi", "temp_air", "wind"])
        for rec in records:
            writer.writerow(
                [
                    rec["year"],
                    rec["month"],
                    rec["day"],
                    rec["hour"],
                    rec["ghi"],
                    rec["temp_air"],
                    rec["wind"],
                ]
            )


def make_flat_dem(path: pathlib.Path, n: int = 9, cell_deg: float = 0.001):
    """Small flat geographic DEM centered on (LON, LAT)."""
    data = np.full((n, n), 2400.0, dtype=np.float64)
    transform = from_origin(LON - (n / 2) * cell_deg, LAT + (n / 2) * cell_deg, cell_deg, cell_deg)
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


def run_solarpv(weather_csv: pathlib.Path, dem: pathlib.Path, out_prefix: pathlib.Path):
    cmd = [
        "cargo",
        "run",
        "--release",
        "-p",
        "solarpv-cli",
        "--",
        "--dem",
        str(dem),
        "--lat",
        str(LAT),
        "--lon",
        str(LON),
        "--date",
        f"{YEAR}-01-01",
        "--out-prefix",
        str(out_prefix),
        "--step",
        "60",
        "--annual",
        "--day-stride",
        "1",
        "--mount",
        "tilt",
        "--tilt",
        str(TILT),
        "--surface-azimuth",
        str(AZIMUTH),
        "--pdc0",
        str(PDC0),
        "--gamma",
        str(GAMMA_PDC),
        "--loss",
        str(SYSTEM_LOSS / 100.0),
        "--weather",
        str(weather_csv),
        "--weather-dt",
        "1.0",
    ]
    subprocess.run(cmd, check=True, cwd=pathlib.Path(__file__).parent.parent)


def read_raster_value(path: pathlib.Path) -> float:
    with rasterio.open(path) as src:
        data = src.read(1)
        r, c = data.shape[0] // 2, data.shape[1] // 2
        return float(data[r, c])


def pvlib_chain_from_pvgis_ghi(records, elevation=2400.0):
    """Reproduce the solarpv-rs point chain using pvlib and the same GHI input.

    Returns annual POA Wh/m² and AC Wh for a 1 kW system.
    """
    loc = pvlib.location.Location(LAT, LON, altitude=elevation)
    times = pd.to_datetime(
        [f"{r['year']}-{r['month']:02d}-{r['day']:02d} {r['hour']:02d}:00:00" for r in records]
    ).tz_localize("UTC")
    ghi = np.array([r["ghi"] for r in records])
    temp_air = np.array([r["temp_air"] for r in records])
    wind = np.array([r["wind"] for r in records])

    solpos = loc.get_solarposition(times)
    app_zen = solpos["apparent_zenith"].to_numpy()
    azimuth = solpos["azimuth"].to_numpy()
    airmass = pvlib.atmosphere.get_relative_airmass(app_zen)
    dni_extra = pvlib.irradiance.get_extra_radiation(times, method="spencer")
    doy = times.dayofyear.to_numpy()

    erbs = pvlib.irradiance.erbs(ghi, app_zen, doy)
    dni = np.nan_to_num(np.asarray(erbs["dni"]))
    dhi = np.nan_to_num(np.asarray(erbs["dhi"]))

    poa = pvlib.irradiance.get_total_irradiance(
        TILT, AZIMUTH, app_zen, azimuth, dni, ghi, dhi,
        dni_extra=dni_extra, airmass=airmass, albedo=0.25, model="perez",
    )
    poa_global = np.nan_to_num(poa["poa_global"].to_numpy())

    tparams = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"]["open_rack_glass_glass"]
    tcell = pvlib.temperature.sapm_cell(poa_global, temp_air, wind, **tparams)
    dc = pvlib.pvsystem.pvwatts_dc(poa_global, tcell, PDC0, GAMMA_PDC) * (1.0 - SYSTEM_LOSS / 100.0)
    ac = np.nan_to_num(pvlib.inverter.pvwatts(dc, PDC0 / 1.1))

    return float(poa_global.sum()), float(ac.sum())


def main():
    base = pathlib.Path(__file__).parent
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="pvgis_check_"))

    print("Fetching PVGIS irradiance series ...")
    rad = pvgis_hourly(pvcalculation=0)
    records = []
    for h in rad["outputs"]["hourly"]:
        y, mo, d, hr = parse_time(h["time"])
        records.append(
            {
                "year": y,
                "month": mo,
                "day": d,
                "hour": hr,
                "ghi": float(h["G(i)"]),
                "temp_air": float(h["T2m"]),
                "wind": float(h["WS10m"]),
            }
        )
    weather_csv = tmp / "weather.csv"
    write_weather_csv(weather_csv, records)
    print(f"  wrote {len(records)} hourly records to {weather_csv}")

    print("Fetching PVGIS PV calculation ...")
    pv = pvgis_hourly(pvcalculation=1)
    pvgis_p_wh = sum(float(h["P"]) for h in pv["outputs"]["hourly"])
    pvgis_poa_wh = sum(float(h["G(i)"]) for h in pv["outputs"]["hourly"])
    print(f"  PVGIS annual AC energy: {pvgis_p_wh/1000:.2f} kWh")
    print(f"  PVGIS annual POA:       {pvgis_poa_wh/1000:.2f} kWh/m²")

    print("Running pvlib chain on PVGIS GHI ...")
    pvlib_poa_wh, pvlib_ac_wh = pvlib_chain_from_pvgis_ghi(records)
    print(f"  pvlib annual AC energy: {pvlib_ac_wh/1000:.2f} kWh")
    print(f"  pvlib annual POA:       {pvlib_poa_wh/1000:.2f} kWh/m²")

    dem = tmp / "dem.tif"
    make_flat_dem(dem)
    out_prefix = tmp / "out"
    print("Running solarpv-cli ...")
    run_solarpv(weather_csv, dem, out_prefix)

    sol_poa_wh = read_raster_value(pathlib.Path(f"{out_prefix}_poa.tif"))
    sol_ac_wh = read_raster_value(pathlib.Path(f"{out_prefix}_ac.tif"))
    sy = read_raster_value(pathlib.Path(f"{out_prefix}_specific_yield.tif"))
    print(f"  solarpv annual AC energy: {sol_ac_wh/1000:.2f} kWh")
    print(f"  solarpv annual POA:       {sol_poa_wh/1000:.2f} kWh/m²")
    print(f"  solarpv specific yield:   {sy:.2f} kWh/kWp")

    print()
    print(f"solarpv vs pvlib (same GHI):")
    print(f"  AC:  {(sol_ac_wh/pvlib_ac_wh - 1)*100:+.2f} %")
    print(f"  POA: {(sol_poa_wh/pvlib_poa_wh - 1)*100:+.2f} %")
    print(f"solarpv vs PVGIS:")
    print(f"  AC:  {(sol_ac_wh/pvgis_p_wh - 1)*100:+.2f} %")
    print(f"  POA: {(sol_poa_wh/pvgis_poa_wh - 1)*100:+.2f} %")

    report = {
        "location": {"lat": LAT, "lon": LON, "year": YEAR},
        "pvgis": {"annual_ac_wh": pvgis_p_wh, "annual_poa_wh": pvgis_poa_wh},
        "pvlib_same_ghi": {"annual_ac_wh": pvlib_ac_wh, "annual_poa_wh": pvlib_poa_wh},
        "solarpv": {"annual_ac_wh": sol_ac_wh, "annual_poa_wh": sol_poa_wh},
    }
    report_path = base / "pvgis_check_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\nReport written to {report_path}")


if __name__ == "__main__":
    main()
