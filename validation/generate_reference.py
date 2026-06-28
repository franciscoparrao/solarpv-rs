#!/usr/bin/env python3
"""Generate pvlib reference cases for validating solarpv-core.

Produces ``reference.json`` consumed by ``crates/core/tests/pvlib_validation.rs``.
Each sample is a clear-sky instant over a point in the Atacama desert (northern
Chile), with every intermediate quantity pvlib computes so that each Rust module
can be validated against the exact same inputs.

Run inside a venv with pvlib installed:
    python3 validation/generate_reference.py
"""
import json
import pathlib

import numpy as np
import pandas as pd
import pvlib

# Atacama desert point — Chile's solar heartland (near Crucero, ~2400 m).
LAT, LON, ALT = -23.0, -69.0, 2400.0
ALBEDO = 0.25            # bright desert ground
SURF_TILT = 23.0         # tilt ≈ |latitude|, equator-facing
SURF_AZIM = 0.0          # North-facing in the Southern Hemisphere (0° = N)

# Two clear-sky days (summer + winter), hourly, in UTC. Chile std time = UTC-4
# (-3 in summer DST); we work in UTC throughout and let pvlib handle geometry.
DAYS = ["2026-01-15", "2026-06-21"]


def build_samples():
    loc = pvlib.location.Location(LAT, LON, altitude=ALT)
    samples = []
    for day in DAYS:
        times = pd.date_range(f"{day} 00:00", f"{day} 23:00", freq="1h", tz="UTC")
        solpos = loc.get_solarposition(times)        # NREL SPA
        cs = loc.get_clearsky(times, model="ineichen")  # ghi, dni, dhi
        dni_extra = pvlib.irradiance.get_extra_radiation(times, method="spencer")
        airmass = pvlib.atmosphere.get_relative_airmass(solpos["apparent_zenith"])

        for t in times:
            app_zen = float(solpos["apparent_zenith"][t])
            # Only keep daylight samples (sun above horizon).
            if app_zen >= 90.0:
                continue
            ghi = float(cs["ghi"][t])
            dni = float(cs["dni"][t])
            dhi = float(cs["dhi"][t])
            dnie = float(dni_extra[t])
            am = float(airmass[t])
            doy = int(t.dayofyear)

            # Erbs decomposition from GHI (independent of the clear-sky DNI/DHI).
            erbs = pvlib.irradiance.erbs(ghi, app_zen, doy)

            poa = {}
            for model in ("isotropic", "haydavies", "perez"):
                r = pvlib.irradiance.get_total_irradiance(
                    SURF_TILT, SURF_AZIM,
                    app_zen, float(solpos["azimuth"][t]),
                    dni, ghi, dhi,
                    dni_extra=dnie, airmass=am, albedo=ALBEDO, model=model,
                )
                poa[model] = {
                    "global": float(r["poa_global"]),
                    "direct": float(r["poa_direct"]),
                    "sky_diffuse": float(r["poa_diffuse"]) - float(r["poa_ground_diffuse"]),
                    "ground_diffuse": float(r["poa_ground_diffuse"]),
                }

            # PV chain on the realistic Perez POA.
            poa_g = poa["perez"]["global"]
            temp_air, wind = 18.0, 2.0
            tparams = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][
                "open_rack_glass_glass"
            ]
            tcell = float(pvlib.temperature.sapm_cell(poa_g, temp_air, wind, **tparams))
            dc = float(pvlib.pvsystem.pvwatts_dc(poa_g, tcell, 1000.0, -0.004))
            ac = float(pvlib.inverter.pvwatts(dc, 1000.0 / 1.1))

            samples.append({
                "utc": {"year": t.year, "month": t.month, "day": t.day,
                        "hour": t.hour, "minute": t.minute, "second": t.second},
                "doy": doy,
                "lat": LAT, "lon": LON,
                "tilt": SURF_TILT, "surface_azimuth": SURF_AZIM, "albedo": ALBEDO,
                # pvlib solar position (NREL SPA) — the oracle for solpos.
                "pv_apparent_zenith": app_zen,
                "pv_azimuth": float(solpos["azimuth"][t]),
                "pv_elevation": float(solpos["apparent_elevation"][t]),
                # irradiance inputs / oracles.
                "ghi": ghi, "dni": dni, "dhi": dhi,
                "dni_extra": dnie, "airmass": am,
                "erbs_dni": float(erbs["dni"]), "erbs_dhi": float(erbs["dhi"]),
                "poa": poa,
                # pv chain.
                "temp_air": temp_air, "wind": wind,
                "tcell": tcell, "dc": dc, "ac": ac,
            })
    return samples


def main():
    samples = build_samples()
    out = {
        "meta": {
            "source": "pvlib " + pvlib.__version__,
            "location": {"lat": LAT, "lon": LON, "altitude": ALT,
                         "name": "Atacama desert, Chile"},
            "n_samples": len(samples),
        },
        "samples": samples,
    }
    path = pathlib.Path(__file__).with_name("reference.json")
    path.write_text(json.dumps(out, indent=2))
    print(f"wrote {len(samples)} samples to {path}")


if __name__ == "__main__":
    main()
