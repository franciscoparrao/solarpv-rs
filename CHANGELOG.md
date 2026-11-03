# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/).

## [0.4.1] - 2026-11-03

### Fixed
- `--tile` combined with `--ghi` / `--ghi-monthly` read the observed GHI with
  tile-local indices, so every tile used the GHI of the DEM's north-west corner.
  On a coastal DEM whose corner is sea the whole yield map came back as zero;
  elsewhere each tile silently used another area's clearness index. Each tile now
  reads the observed GHI at its position in the full DEM. Untiled and clear-sky
  runs are unchanged bit for bit. Reported and fixed by Braulio Bravo and Alonso
  Henríquez (#1), with a regression test.

### Changed
- Documentation no longer claims that tiled runs are identical to untiled ones.
  Horizon rays whose step lands exactly on a half cell round by absolute index,
  so a few horizon angles near tile seams can differ (≲0.1 % in specific yield
  on rugged test terrain).
- Python package version aligned with the crate (0.4.1).

## [0.4.0] - 2026-09-19

### Added
- `--tile` to process large DEMs in memory-bounded tiles with a horizon halo.
- `--ghi` / `--ghi-monthly` (`--ghi-unit kwh|wh`) to rescale the clear-sky model
  per cell to an observed mean-daily GHI raster before decomposition.
- Tiling and observed-GHI options in the Python bindings.

### Fixed
- `--per-cell-lat` read UTM northings as degrees of latitude. The engine now
  inverts UTM per cell and validates the DEM's CRS.

[0.4.1]: https://github.com/franciscoparrao/solarpv-rs/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/franciscoparrao/solarpv-rs/releases/tag/v0.4.0
