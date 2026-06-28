# solarpv-rs — Potencial solar fotovoltaico sobre terreno (Rust, "PVGIS lite")

> **Estado:** IDEA (sin código). Creado 2026-06-10.
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
- [ ] (v0.2) Series TMY; tracking; pérdidas detalladas.

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
