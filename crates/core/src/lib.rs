//! # solarpv-core
//!
//! Core physics for terrain-aware photovoltaic potential ("PVGIS lite").
//!
//! This crate implements the *point* model — the validatable chain that turns a
//! location, time and sky condition into PV energy:
//!
//! 1. [`solpos`] — solar geometry (position on the sky, angle of incidence).
//! 2. [`irradiance`] — decomposition of global irradiance and transposition to
//!    the plane of array (POA).
//! 3. [`pv`] — cell temperature, DC power and AC energy yield.
//!
//! The gridded step (mapping POA / yield over a DEM, reusing SurtGIS terrain and
//! horizon rasters) is gated behind the `terrain` feature and built on top of
//! this point core, so the physics can be validated against `pvlib` in isolation.
//!
//! ## Conventions
//!
//! Angles in public APIs are in **degrees** to match `pvlib`. Azimuth follows the
//! solar/`pvlib` convention: measured clockwise from North, so 0° = N, 90° = E,
//! 180° = S, 270° = W. Irradiance is in W/m²; energy in Wh.

pub mod error;
pub mod irradiance;
pub mod pv;
pub mod solpos;

/// Gridded PV potential over a DEM, reusing SurtGIS terrain rasters.
#[cfg(feature = "terrain")]
pub mod grid;

pub use error::{Error, Result};
