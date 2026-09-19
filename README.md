# solarpv-rs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22845645.svg)](https://doi.org/10.5281/zenodo.22845645)

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
| `spa` | NREL SPA (Reda & Andreas 2004), high accuracy | `solarposition.spa_python` (< 2 arcsec) |
| `irradiance` | Erbs decomposition, POA transposition (isotropic / Hay-Davies / Perez), Haurwitz clear-sky | `irradiance.erbs`, `irradiance.get_total_irradiance` |
| `pv` | SAPM cell temperature, PVWatts DC + inverter, energy integration | `temperature.sapm_cell`, `pvsystem.pvwatts_dc`, `inverter.pvwatts` |
| `tracking` | Single- and dual-axis tracker geometry | `tracking.singleaxis` (≤ 1e-4°); dual-axis ideal sun-following |
| `losses` | IAM (ASHRAE/Martín-Ruiz/physical), PVWatts loss breakdown, SAPM spectral mismatch | `iam.*`, `pvsystem.pvwatts_losses`, `spectrum.spectral_factor_sapm` (≤ 1e-6) |

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

Optional `--svf` applies a **topographic sky-view factor** (`SVF = 1 − mean(sin²(horizon))`)
to the diffuse-sky component, reducing the diffuse irradiance received in valleys
or near ridges. It is off by default so the default clear-sky path remains
numerically identical to the pvlib validation.

Optional `--features gdal` switches the GeoTIFF I/O backend to GDAL through
SurtGIS (`surtgis-core/gdal`) for broader format support; the default native
backend keeps the engine free of any GDAL dependency.

## Layout

```
crates/core/            solarpv-core: solar geometry, irradiance, PV models, grid
crates/core/tests/      pvlib parity + grid integration tests
crates/cli/             solarpv-cli: the `solarpv` binary (DEM in, rasters out)
crates/python/          solarpv-py: Python bindings via PyO3 / maturin
validation/             pvlib reference generator + reference.json + PVGIS cross-check
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
`--mount tilt --tilt 25 --surface-azimuth 0` (fixed racks),
`--mount tracker --gcr 0.3` (horizontal single-axis tracker with backtracking), or
`--mount dual-axis --max-tilt 90` (full hemispherical dual-axis tracker).

For terrain-aware diffuse reduction add `--svf`.

To use real irradiance instead of clear-sky, pass a CSV with
`--weather tmy.csv --weather-dt 1.0` (header
`year,month,day,hour,ghi[,dni,dhi,temp_air,wind]`, UTC).

## Python bindings

A PyO3 package (`solarpv`) exposes the validated point chain and the gridded
engine to Python. Build and install it with [maturin](https://www.maturin.rs):

```bash
cd crates/python
maturin develop --features terrain      # editable install
pytest tests -v                         # Python-side tests
```

Usage example:

```python
import solarpv

# Solar position
sol = solarpv.solar_position(-23.0, -69.0, 2026, 6, 21, 12)

# Plane-of-array irradiance
poa = solarpv.poa_irradiance(
    ghi=800, dni=600, dhi=150,
    tilt=23, surface_azimuth=0,
    solar_zenith=20, solar_azimuth=0,
    albedo=0.25, dni_extra=1366, airmass=1.2,
    model="perez",
)

# PV power for one time step
pv = solarpv.ac_power(poa["global"], temp_air=18, wind=2)

# Gridded potential over a DEM (requires --features terrain)
grid = solarpv.pv_potential(
    "dem.tif", lat=-23.0, lon=-69.0,
    year=2026, month=6, day=21,
    step=30, mount="tilt", tilt=23.0,
)
```

## Running

```bash
cargo test --workspace                   # Rust unit + integration tests

# regenerate the pvlib reference (needs a venv with pvlib):
python3 validation/generate_reference.py

# PVGIS cross-check:
python3 validation/pvgis_check.py

# Python bindings tests:
cd crates/python && maturin develop --features terrain && pytest tests -v
```

## Cross-check vs PVGIS

A direct comparison against the [PVGIS](https://re.jrc.ec.europa.eu/pvgis.html)
v5.2 API for the Atacama reference point (ERA5, 2020, 1 kWp fixed tilt 23°,
North-facing, c-Si, 14 % system loss) shows that `solarpv-rs` tracks **pvlib**
rather than PVGIS when fed the same GHI:

| Source | Annual POA (kWh/m²) | Annual AC (kWh/kWp) |
|--------|--------------------:|--------------------:|
| PVGIS v5.2 | 2089 | 1608 |
| pvlib (same PVGIS GHI, Erbs + Perez + PVWatts) | 2831 | 2180 |
| **solarpv-rs** (same PVGIS GHI) | **2869** | **2196** |

`solarpv-rs` differs from pvlib by **< 2 %** on the same GHI input, confirming
internal consistency. PVGIS itself reports ~37 % lower yield because it uses its
own satellite-derived decomposition/transposition and additional operational
loss/attenuation modelling not reproduced in the open pvlib chain.

Run the cross-check with:

```bash
python3 validation/pvgis_check.py
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
  Also adds fixed-tilt mounting (`--mount tilt`) and dual-axis tracking
  (`--mount dual-axis`).
- **v0.2** ✅ measured / TMY irradiance series (`--weather tmy.csv`): drive the
  engine from real GHI (+ optional DNI/DHI/temp/wind) instead of clear-sky.
  Feeding back the internal Haurwitz series reproduces the clear-sky path to 1e-6.
- **v0.2** ✅ NREL SPA solar position (`--spa`, validated vs pvlib to < 2 arcsec):
  the high-accuracy alternative to Michalsky, optional per run.
- **v0.2** ✅ detailed losses (validated vs pvlib to 1e-6): incidence-angle
  modifiers (`--iam ashrae|martin-ruiz|physical`, beam reflection loss ≈ −2 %/yr
  fixed tilt), the PVWatts loss breakdown, and SAPM spectral mismatch (`--spectral`);
  `--loss` overrides the DC derate.
- **v0.2** ✅ GDAL I/O backend as an opt-in feature (`--features gdal`) delegated to
  `surtgis-core/gdal`; default build stays GDAL-free.
- **v0.2** ✅ cross-check against PVGIS v5.2 documented.
- **v0.2** ✅ Python bindings (PyO3 / maturin) for the point chain and gridded engine.
- **v0.3** ✅ topographic sky-view factor (`--svf`) reducing diffuse-sky irradiance
  in obstructed terrain.
- **v0.3 (next)** WASM bindings for an interactive browser estimator; publication prep
  (Zenodo + Renewable Energy manuscript).

## Conventions

Angles in degrees (azimuth clockwise from North: 0°=N, 90°=E, 180°=S). Irradiance
in W/m², energy in Wh. Southern-Hemisphere arrays face North (`surface_azimuth = 0`).

## Reproducibility & citation

**Toolchain.** Rust (edition 2024); build with a recent stable `cargo`. The
gridded engine needs the `terrain` feature, which reuses the sibling SurtGIS
crates via path dependencies (`../surtgis`), no system GDAL required.

```sh
# build the CLI and run the test suite (point chain + gridded + pvlib parity)
cargo build --release -p solarpv-cli --features terrain
cargo test --features terrain
```

**Minimal run** (clear-sky annual PV potential over a DEM, fixed racks at 23°):

```sh
solarpv --dem dem_utm19s.tif --annual --mount tilt --tilt 23         --lat -23.6 --lon -69.5 --date 2026-01-01 --out-prefix out
# large scenes: add --tile 1024 to bound memory; --per-cell-lat for wide spans;
# rescale to observed irradiance with --ghi ghi_annual.tif (Explorador Solar).
```

**Determinism.** The physics is deterministic: no random state, and the rayon
parallelism only sums independent per-cell contributions, so results are
reproducible across runs and core counts. Numerical parity against
[`pvlib`](https://pvlib-python.readthedocs.io) is checked in
`crates/core/tests/pvlib_validation.rs` from a reference generated by
`validation/generate_reference.py`.

**Cite.** This release is archived at [doi:10.5281/zenodo.22845645](https://doi.org/10.5281/zenodo.22845645); see `CITATION.cff`. Cite the specific version you used.

## License

MIT OR Apache-2.0.
