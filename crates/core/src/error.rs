//! Error types for `solarpv-core`.

use thiserror::Error;

/// Result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can arise from invalid inputs to the physical models.
#[derive(Debug, Error)]
pub enum Error {
    /// A parameter was outside its physically meaningful range.
    #[error("invalid parameter `{name}`: {value} (expected {expected})")]
    InvalidParameter {
        /// Parameter name.
        name: &'static str,
        /// Offending value, formatted.
        value: String,
        /// Human-readable description of the valid range.
        expected: &'static str,
    },

    /// A calendar date was not valid (e.g. month 13).
    #[error("invalid date: {0}")]
    InvalidDate(String),

    /// A terrain/raster computation (slope, aspect, horizon) failed.
    #[error("terrain computation failed: {0}")]
    Terrain(String),
}

impl Error {
    /// Helper to build an [`Error::InvalidParameter`].
    pub(crate) fn param(name: &'static str, value: impl std::fmt::Display, expected: &'static str) -> Self {
        Error::InvalidParameter {
            name,
            value: value.to_string(),
            expected,
        }
    }
}
