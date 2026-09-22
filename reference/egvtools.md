# egvtools: High-resolution Ecogeographical Variable Workflows

`egvtools` provides a coherent set of wrappers and utilities that make
large-scale **EGV** creation reproducible and pleasant on real datasets.
The package leans on robust building blocks—`terra`, `sf`, `sfarrow`,
`exactextractr`, and `whitebox`—and standardizes I/O, naming
conventions, and multi-scale zonal statistics so your pipelines are
repeatable across machines and projects.

## Core workflow (generalizable)

Functions that form the backbone of your analyses (thin wrappers over
existing libs, but opinionated for stability, speed and consistency):

- [`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md)
  — rasterize polygons to template, handle background/mask.

- [`downscale2egv()`](https://aavotins.github.io/egvtools/reference/downscale2egv.md)
  — downscale coarse rasters to template grid and optionally smooth with
  IDW.

- [`distance2egv()`](https://aavotins.github.io/egvtools/reference/distance2egv.md)
  — distances to features with optional gap filling at the edges.

- [`input2egv()`](https://aavotins.github.io/egvtools/reference/input2egv.md)
  — normalize/align inputs to EGV outputs with guards.

- [`landscape_function()`](https://aavotins.github.io/egvtools/reference/landscape_function.md)
  — landscapemetrics landscape-level per-zone metrics, tiled.

- [`radius_function()`](https://aavotins.github.io/egvtools/reference/radius_function.md)
  — multi-scale zonal statistics (dense/sparse).

- [`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)
  — tile template grids for chunked processing.

- [`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md)
  — precompute buffered tiles for multiple radii.

- [`create_backgrounds()`](https://aavotins.github.io/egvtools/reference/create_backgrounds.md)
  — build consistent background rasters/values.

## Reproducibility helpers

Utilities that set up inputs and structure so our (project HiQBioDiv)
results can be reproduced:

- [`download_raster_templates()`](https://aavotins.github.io/egvtools/reference/download_raster_templates.md)
  — fetch template rasters (Zenodo) to canonical paths.

- [`download_vector_templates()`](https://aavotins.github.io/egvtools/reference/download_vector_templates.md)
  — fetch template vector grids/points.

## Design principles

- **Tiled, RAM-aware I/O**; workers open data by path (avoid big
  globals) and caches them if necessary (avoid big I/O).

- **Deterministic filenames** and strict layer/radius naming.

- **Guards** for empty/invalid geometries, coordinate reference systems,
  naming and all-NA joins.

- **Cross-platform parallel** via `{future}/{furrr}` when enabled.

## Options

Package defaults are set on load (see `zzz-options.R`). Users may
override:

- `options(egvtools.future_plan = "sequential")`

- `options(egvtools.progress = TRUE)`

&nbsp;

- `egvtools.future_plan`: default parallel plan name (e.g.
  `"sequential"`, `"multisession"`).

- `egvtools.progress`: logical, show progress bars (`TRUE`/`FALSE`).

## Getting started

1.  Run `download_*_templates()` to fetch canonical inputs.

2.  Use
    [`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)
    /
    [`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md)
    for scalable chunks.

3.  Produce EGVs at site scale with
    [`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md),
    [`downscale2egv()`](https://aavotins.github.io/egvtools/reference/downscale2egv.md),
    [`distance2egv()`](https://aavotins.github.io/egvtools/reference/distance2egv.md),
    [`input2egv()`](https://aavotins.github.io/egvtools/reference/input2egv.md),
    and summarise from larger scales with
    [`landscape_function()`](https://aavotins.github.io/egvtools/reference/landscape_function.md)
    /
    [`radius_function()`](https://aavotins.github.io/egvtools/reference/radius_function.md).

## See also

Useful links:

- GitHub repo: <https://github.com/aavotins/egvtools>

- Package site: <https://aavotins.github.io/egvtools/>

- Issues / bugs: <https://github.com/aavotins/egvtools/issues>

## Author

**Maintainer**: Andris Avotiņš <andris.avotins@lu.lv>

Authors:

- Andris Avotiņš <andris.avotins@lu.lv>
