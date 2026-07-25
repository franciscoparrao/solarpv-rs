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

            # Haurwitz clear-sky GHI for this zenith.
            hz = pvlib.clearsky.haurwitz(pd.Series([app_zen]))
            haurwitz_ghi = float(np.asarray(hz).ravel()[0])

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
                "haurwitz_ghi": haurwitz_ghi,
                "poa": poa,
                # pv chain.
                "temp_air": temp_air, "wind": wind,
                "tcell": tcell, "dc": dc, "ac": ac,
            })
    return samples


def build_annual_point(year=2026, tilt=23.0, surf_azim=0.0,
                       temp_air=18.0, wind=2.0, losses=0.14):
    """Full-year HOURLY clear-sky → fixed-tilt POA → PVWatts integration at the
    Atacama point, using exactly solarpv-core's model chain (Haurwitz GHI, Erbs,
    Perez, SAPM cell temp, PVWatts DC with system losses, PVWatts inverter).

    This is the end-to-end oracle: the Rust engine should reproduce the annual
    POA insolation and AC energy from the same 8760-hour integration.
    """
    loc = pvlib.location.Location(LAT, LON, altitude=ALT)
    times = pd.date_range(f"{year}-01-01 00:00", f"{year}-12-31 23:00",
                          freq="1h", tz="UTC")
    solpos = loc.get_solarposition(times)
    app_zen = solpos["apparent_zenith"]
    azimuth = solpos["azimuth"]
    doy = times.dayofyear

    ghi = pvlib.clearsky.haurwitz(app_zen)
    ghi = np.asarray(ghi).ravel()
    app_zen = app_zen.to_numpy()
    azimuth = azimuth.to_numpy()

    erbs = pvlib.irradiance.erbs(ghi, app_zen, np.asarray(doy))
    dni = np.nan_to_num(np.asarray(erbs["dni"]))
    dhi = np.nan_to_num(np.asarray(erbs["dhi"]))
    dni_extra = np.asarray(pvlib.irradiance.get_extra_radiation(times, method="spencer"))
    airmass = np.nan_to_num(np.asarray(pvlib.atmosphere.get_relative_airmass(app_zen)), nan=0.0)

    poa = pvlib.irradiance.get_total_irradiance(
        tilt, surf_azim, app_zen, azimuth, dni, ghi, dhi,
        dni_extra=dni_extra, airmass=airmass, albedo=ALBEDO, model="perez",
    )
    poa_global = np.nan_to_num(np.asarray(poa["poa_global"]))

    tparams = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][
        "open_rack_glass_glass"]
    tcell = pvlib.temperature.sapm_cell(poa_global, temp_air, wind, **tparams)
    dc = pvlib.pvsystem.pvwatts_dc(poa_global, tcell, 1000.0, -0.004) * (1.0 - losses)
    ac = np.nan_to_num(np.asarray(pvlib.inverter.pvwatts(dc, 1000.0 / 1.1)))

    return {
        "year": year, "tilt": tilt, "surface_azimuth": surf_azim,
        "temp_air": temp_air, "wind": wind, "system_losses": losses,
        "albedo": ALBEDO, "pdc0": 1000.0, "gamma_pdc": -0.004,
        "annual_poa_wh": float(poa_global.sum()),   # dt = 1 h
        "annual_ac_wh": float(ac.sum()),
        "annual_specific_yield": float(ac.sum() / 1000.0),  # kWh/kWp
    }


def build_tracking():
    """Single-axis tracker cases across a range of sun positions and configs,
    for validating tracking::single_axis against pvlib.tracking.singleaxis."""
    cases = []
    configs = [
        # (axis_tilt, axis_azimuth, max_angle, backtrack, gcr)
        (0.0, 0.0, 90.0, False, 2.0 / 7.0),
        (0.0, 0.0, 60.0, True, 0.5),
        (0.0, 0.0, 45.0, True, 0.35),
        (20.0, 180.0, 90.0, False, 2.0 / 7.0),
        (10.0, 200.0, 90.0, True, 0.4),
    ]
    suns = [(20.0, 90.0), (45.0, 120.0), (60.0, 80.0), (75.0, 95.0),
            (85.0, 270.0), (50.0, 200.0), (30.0, 180.0)]
    for (at, aa, ma, bt, gcr) in configs:
        for (zen, az) in suns:
            r = pvlib.tracking.singleaxis(zen, az, axis_tilt=at, axis_azimuth=aa,
                                          max_angle=ma, backtrack=bt, gcr=gcr)
            cases.append({
                "axis_tilt": at, "axis_azimuth": aa, "max_angle": ma,
                "backtrack": bt, "gcr": gcr, "zenith": zen, "azimuth": az,
                "tracker_theta": float(np.asarray(r["tracker_theta"]).ravel()[0]),
                "surface_tilt": float(np.asarray(r["surface_tilt"]).ravel()[0]),
                "surface_azimuth": float(np.asarray(r["surface_azimuth"]).ravel()[0]),
                "aoi": float(np.asarray(r["aoi"]).ravel()[0]),
            })
    return cases


def build_spa():
    """NREL SPA reference (pvlib.solarposition.spa_python, default delta_t=67,
    pressure=101325 Pa, temperature=12 C) across the year, for validating the
    Rust spa module to sub-arcminute accuracy."""
    cases = []
    times = pd.date_range("2026-01-01 00:00", "2026-12-28 21:00", freq="37h", tz="UTC")
    sp = pvlib.solarposition.spa_python(
        times, LAT, LON, altitude=ALT, pressure=101325.0, temperature=12.0, delta_t=67.0,
    )
    for ts in times:
        zen = float(sp["apparent_zenith"][ts])
        if zen >= 90.0:
            continue  # below horizon: azimuth ill-defined
        cases.append({
            "utc": {"year": ts.year, "month": ts.month, "day": ts.day,
                    "hour": ts.hour, "minute": ts.minute, "second": ts.second},
            "lat": LAT, "lon": LON, "altitude": ALT,
            "apparent_zenith": zen,
            "zenith": float(sp["zenith"][ts]),
            "azimuth": float(sp["azimuth"][ts]),
            "apparent_elevation": float(sp["apparent_elevation"][ts]),
        })
    return cases


def build_losses():
    """IAM curves (ashrae, martin_ruiz, physical) and the PVWatts loss
    breakdown, for validating the losses module against pvlib."""
    aois = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 85.0]
    iam = {
        "ashrae": [float(pvlib.iam.ashrae(a, b=0.05)) for a in aois],
        "martin_ruiz": [float(pvlib.iam.martin_ruiz(a, a_r=0.16)) for a in aois],
        "physical": [float(pvlib.iam.physical(a, n=1.526, K=4.0, L=0.002)) for a in aois],
    }
    pvwatts_losses = float(pvlib.pvsystem.pvwatts_losses())  # defaults → ~14.08 %
    return {"aois": aois, "iam": iam, "pvwatts_losses_pct": pvwatts_losses}


def build_spectral():
    """SAPM spectral mismatch factor f1 vs absolute airmass, for validating
    losses::SpectralLoss against pvlib.spectrum.spectral_factor_sapm.
    """
    # Use the same crystalline-silicon module the Rust side defaults to.
    module = pvlib.pvsystem.retrieve_sam("SandiaMod")["Canadian_Solar_CS5P_220M___2009_"]
    ams = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
    # Absolute airmass at the Atacama point (2400 m, ~775 hPa).
    pressure = pvlib.atmosphere.alt2pres(ALT)
    am_abs = [float(pvlib.atmosphere.get_absolute_airmass(am, pressure)) for am in ams]
    factors = [
        float(pvlib.spectrum.spectral_factor_sapm(am, module))
        for am in am_abs
    ]
    return {
        "module": "Canadian_Solar_CS5P_220M___2009_",
        "coefficients": {
            "A0": float(module["A0"]), "A1": float(module["A1"]),
            "A2": float(module["A2"]), "A3": float(module["A3"]),
            "A4": float(module["A4"]),
        },
        "pressure_pa": float(pressure),
        "airmass_absolute": am_abs,
        "factors": factors,
    }


def main():
    samples = build_samples()
    annual_point = build_annual_point()
    tracking = build_tracking()
    spa_cases = build_spa()
    losses = build_losses()
    spectral = build_spectral()
    out = {
        "meta": {
            "source": "pvlib " + pvlib.__version__,
            "location": {"lat": LAT, "lon": LON, "altitude": ALT,
                         "name": "Atacama desert, Chile"},
            "n_samples": len(samples),
        },
        "samples": samples,
        "annual_point": annual_point,
        "tracking": tracking,
        "spa": spa_cases,
        "losses": losses,
        "spectral": spectral,
    }
    path = pathlib.Path(__file__).with_name("reference.json")
    path.write_text(json.dumps(out, indent=2))
    print(f"wrote {len(samples)} samples to {path}")


if __name__ == "__main__":
    main()
