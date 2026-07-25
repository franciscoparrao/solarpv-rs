# solarpv-rs — Potencial solar fotovoltaico sobre terreno (Rust, "PVGIS lite")

> **Estado:** FUNCIONAL v0.1/v0.3 (actualizado 2026-07-24). ~4.000 LOC, 50 tests Rust + 7 tests Python. Implementado: posición solar Michalsky + NREL SPA, descomposición Erbs, transposición POA (isotropic/Hay-Davies/Perez), cadena PV (SAPM, PVWatts DC+inversor), tracking single-axis con backtracking + dual-axis ideal, pérdidas IAM + spectral mismatch SAPM, paso grillado sobre DEM (feature `terrain`) con **sky-view factor topográfico** (`--svf`), TMY/medidas, bindings Python (PyO3/maturin). Feature `gdal` delegado a `surtgis-core/gdal` (I/O con librerías GDAL del sistema, sin dependencia directa de GDAL en este crate). Sin stubs. Validado: paridad con pvlib en Atacama (<0.5% anual); cross-check vs PVGIS v5.2 muestra que solarpv-rs sigue a pvlib (<2% mismo GHI) mientras PVGIS reporta ~37% menos por modelos propietarios de descomposición/atenuación. GAP respecto a la familia: **falta WASM y publicación**. Sin paper (venue: Renewable Energy). Próximo: WASM o preparar publicación (Zenodo + Renewable Energy).
> Familia de motores Rust del autor: SurtGIS, Hydroflux, Smelt, Anvil, Cantus, Criterium.
> Doc madre: `~/proyectos/ideas-motores-rust.md` (idea K1).

## Qué es
Motor que extiende la radiación solar topográfica a producción fotovoltaica:
horizonte, sombreado, series temporales y yield energético.

## El gap que llena
SurtGIS tiene radiación solar pero no llega a **producción PV**. El campo es
**PVGIS** (servicio web), pvlib (Python). Chile = capital solar mundial → caso
de uso fuerte.

## Alcance MVP (v0.1)
- [ ] Geometría solar (posición, ángulo de incidencia) sobre DEM.
- [ ] Horizonte y sombreado topográfico (reusa viewshed/openness de SurtGIS).
- [ ] Irradiancia POA (plane-of-array) con descomposición difusa.
- [ ] Yield PV (modelo simple temperatura-eficiencia) y energía anual.
- [x] (v0.2) Series TMY; tracking; pérdidas detalladas.
- [x] (v0.2+) Tracking dual-axis; feature GDAL.
- [x] (v0.3) Sky-view factor topográfico para reducción de difusa en terreno obstruido.

## Arquitectura tentativa
- `solarpv-core`: geometría solar, modelos de irradiancia y PV.
- Targets: native + Python (PyO3) + CLI; WASM para estimador interactivo.
- Reusa terrain (slope/aspect/horizonte) de SurtGIS.

## Validación / paridad numérica
Cross-check contra **pvlib** y datos de estaciones solares (ej. red chilena).

## Venue objetivo
**Renewable Energy** o **Solar Energy**.

## Conexiones con tu ecosistema
- **SurtGIS**: slope/aspect/horizonte/radiación como base.
- Proyecto `solar_grupo_4` como posible caso.

## Próximos pasos al retomar
1. Implementar geometría solar + sombreado de horizonte sobre DEM.
2. Calcular POA + yield anual; validar contra pvlib en un punto.
3. Mapear potencial PV de una zona del norte de Chile.
