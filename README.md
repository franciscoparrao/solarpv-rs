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
| `irradiance` | Erbs decomposition, POA transposition (isotropic / Hay-Davies / Perez), Haurwitz clear-sky | `irradiance.erbs`, `irradiance.get_total_irradiance` |
| `pv` | SAPM cell temperature, PVWatts DC + inverter, energy integration | `temperature.sapm_cell`, `pvsystem.pvwatts_dc`, `inverter.pvwatts` |

All transposition and PV-chain quantities match pvlib to ≤ 0.2 % relative.

### Gridded step (feature `terrain`)

`grid::pv_potential` maps the point chain over a DEM, reusing SurtGIS
`slope`/`aspect`/`horizon_angles`: each cell is a ground-following surface with
beam shading from the surrounding topography. It produces per-cell POA insolation
(Wh/m²/day), AC energy (Wh/day) and specific yield (kWh/kWp/day). Enable with
`--features terrain`.

## Layout

```
crates/core/            solarpv-core: solar geometry, irradiance, PV models, grid
crates/core/tests/      pvlib parity + grid integration tests
crates/cli/             solarpv-cli: the `solarpv` binary (DEM in, rasters out)
validation/             pvlib reference generator + reference.json
```

## CLI

```bash
cargo run --release -p solarpv-cli -- \
  --dem dem.tif --lat -23.0 --lon -69.0 --date 2026-06-21 \
  --out-prefix out/atacama --step 15 --sky perez

# writes out/atacama_{poa,ac,specific_yield}.tif (GeoTIFF, georeferenced)
```

Reads a single-band elevation GeoTIFF (native reader, no GDAL needed) and writes
per-cell POA insolation, AC energy and specific yield. `--help` lists all flags
(albedo, system nameplate, temperature coefficient, weather, horizon resolution).

## Running

```bash
cargo test -p solarpv-core            # unit + pvlib parity tests

# regenerate the pvlib reference (needs a venv with pvlib):
python3 validation/generate_reference.py
```

## Roadmap

- **v0.1** ✅ point PV chain validated against pvlib.
- **v0.1** ✅ gridded step over a DEM reusing SurtGIS `horizon_angles` +
  `slope`/`aspect` (feature `terrain`), parallelised with rayon.
- **v0.1** ✅ CLI: read DEM GeoTIFF → write POA / AC / specific-yield rasters.
- **v0.1 (next)** per-cell latitude for large scenes; annual integration
  (multi-day); GDAL feature for broader format support.
- **v0.2** TMY / measured irradiance series, tracking, detailed losses; NREL SPA
  for sub-arcminute solar position; PV potential map of a northern-Chile zone.

## Conventions

Angles in degrees (azimuth clockwise from North: 0°=N, 90°=E, 180°=S). Irradiance
in W/m², energy in Wh. Southern-Hemisphere arrays face North (`surface_azimuth = 0`).

## License

MIT OR Apache-2.0.
