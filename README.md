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
| `tracking` | Single-axis tracker geometry with backtracking | `tracking.singleaxis` (≤ 1e-4°) |

All transposition and PV-chain quantities match pvlib to ≤ 0.2 % relative.
Beyond the per-component checks, an **end-to-end** test reproduces pvlib's full
8760-hour clear-sky → fixed-tilt POA → PVWatts annual integration at the Atacama
point (2054 kWh/kWp·yr) to within **0.5 %** — the only residual being Michalsky
vs NREL SPA solar position.

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

# annual potential (one representative day per month), per-cell latitude:
cargo run --release -p solarpv-cli -- \
  --dem dem.tif --lat -23.0 --lon -69.0 --date 2026-01-01 \
  --out-prefix out/atacama_year --annual --per-cell-lat
```

Reads a single-band elevation GeoTIFF (native reader, no GDAL needed) and writes
per-cell POA insolation, AC energy and specific yield. `--help` lists all flags
(albedo, system nameplate, temperature coefficient, weather, horizon resolution).

Mounting is selectable: `--mount terrain` (default, ground-following),
`--mount tilt --tilt 25 --surface-azimuth 0` (fixed racks), or
`--mount tracker --gcr 0.3` (horizontal single-axis tracker with backtracking).

To use real irradiance instead of clear-sky, pass a CSV with
`--weather tmy.csv --weather-dt 1.0` (header
`year,month,day,hour,ghi[,dni,dhi,temp_air,wind]`, UTC).

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
- **v0.1** ✅ per-cell latitude/longitude for large scenes (`--per-cell-lat`,
  geographic DEMs): solar ephemeris precomputed per step, sky position per cell.
- **v0.1** ✅ annual integration (`--annual`): terrain computed once, year
  sampled by monthly representative days or an N-day stride.
- **v0.1** ✅ end-to-end annual validation vs pvlib (Atacama point, < 0.5 %);
  demonstrated on a real 637×570 Chilean DEM (UTM 19S, ~−32.9°): yield 325–2018
  kWh/kWp·yr, with North-facing slopes (1854) > flat (1726) > South-facing (1509).
- **v0.2** ✅ single-axis tracking (`--mount tracker`, validated vs pvlib): on a
  flat Atacama site the tracker yields ~+36 % over fixed horizontal modules.
  Also adds fixed-tilt mounting (`--mount tilt`).
- **v0.2** ✅ measured / TMY irradiance series (`--weather tmy.csv`): drive the
  engine from real GHI (+ optional DNI/DHI/temp/wind) instead of clear-sky.
  Feeding back the internal Haurwitz series reproduces the clear-sky path to 1e-6.
- **v0.2 (next)** detailed losses; NREL SPA for sub-arcminute solar position;
  GDAL feature; cross-check against PVGIS.

## Conventions

Angles in degrees (azimuth clockwise from North: 0°=N, 90°=E, 180°=S). Irradiance
in W/m², energy in Wh. Southern-Hemisphere arrays face North (`surface_azimuth = 0`).

## License

MIT OR Apache-2.0.
