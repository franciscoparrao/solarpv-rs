# solarpv-rs

Terrain-aware **photovoltaic potential** engine in Rust — a "PVGIS lite" that
extends topographic solar radiation to actual PV production: solar geometry,
plane-of-array irradiance, and energy yield.

Part of the author's Rust geospatial engine family (SurtGIS, Hydroflux, Smelt,
Anvil, Cantus, Criterium). It reuses SurtGIS's terrain / horizon rasters for the
gridded step and fills the gap SurtGIS leaves open: **irradiance → PV energy**.

## Status — v0.1 (point model, validated)

The validatable point chain is implemented in `solarpv-core` and checked for
numerical parity against [`pvlib`](https://pvlib-python.readthedocs.io) at a
point in the Atacama desert (northern Chile):

| Module | What it does | pvlib oracle |
|--------|--------------|--------------|
| `solpos` | Sun position (Michalsky) + angle of incidence | `solarposition.spa_python` (< 0.1° zenith) |
| `irradiance` | Erbs decomposition, POA transposition (isotropic / Hay-Davies / Perez) | `irradiance.erbs`, `irradiance.get_total_irradiance` |
| `pv` | SAPM cell temperature, PVWatts DC + inverter, energy integration | `temperature.sapm_cell`, `pvsystem.pvwatts_dc`, `inverter.pvwatts` |

All transposition and PV-chain quantities match pvlib to ≤ 0.2 % relative.

## Layout

```
crates/core/            solarpv-core: solar geometry, irradiance, PV models
crates/core/tests/      pvlib parity tests
validation/             pvlib reference generator + reference.json
```

## Running

```bash
cargo test -p solarpv-core            # unit + pvlib parity tests

# regenerate the pvlib reference (needs a venv with pvlib):
python3 validation/generate_reference.py
```

## Roadmap

- **v0.1** ✅ point PV chain validated against pvlib.
- **v0.1 (next)** gridded step: map POA / yield over a DEM reusing SurtGIS
  `horizon_angles` + `slope`/`aspect` (feature `terrain`); CLI.
- **v0.2** TMY time series, tracking, detailed losses; NREL SPA for sub-arcminute
  solar position; PV potential map of a northern-Chile zone.

## Conventions

Angles in degrees (azimuth clockwise from North: 0°=N, 90°=E, 180°=S). Irradiance
in W/m², energy in Wh. Southern-Hemisphere arrays face North (`surface_azimuth = 0`).

## License

MIT OR Apache-2.0.
