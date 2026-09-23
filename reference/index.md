# Package index

## Core workflow

Generalizable EGV creators/wrappers

- [`polygon2input()`](https://aavotins.github.io/egvtools/reference/polygon2input.md)
  : Rasterize polygons to a template grid, optionally restrict & cover
  gaps, then write GeoTIFF
- [`downscale2egv()`](https://aavotins.github.io/egvtools/reference/downscale2egv.md)
  : Downscale & Align a Raster to a Template, with optional gap fill &
  IDW smoothing
- [`distance2egv()`](https://aavotins.github.io/egvtools/reference/distance2egv.md)
  : Distance-to-class on the EGV grid
- [`input2egv()`](https://aavotins.github.io/egvtools/reference/input2egv.md)
  : Convert an input raster to an EGV-aligned raster
- [`landscape_function()`](https://aavotins.github.io/egvtools/reference/landscape_function.md)
  : Compute landscape-level metrics per zone (tiled), merge, analyze
  gaps, optionally fill with IDW, and write rasters
- [`radius_function()`](https://aavotins.github.io/egvtools/reference/radius_function.md)
  : Extract and Rasterize Summary Statistics from Buffered Radii Using
  exactextractr
- [`generalized_radius_function()`](https://aavotins.github.io/egvtools/reference/generalized_radius_function.md)
  : Extract and Rasterize Buffered-Radius Summary Statistics

## Reproducibility helpers

Inputs, tiling, and reproducible scaffolding

- [`download_raster_templates()`](https://aavotins.github.io/egvtools/reference/download_raster_templates.md)
  : Download and unpack raster templates
- [`download_vector_templates()`](https://aavotins.github.io/egvtools/reference/download_vector_templates.md)
  : Download and unpack vector templates (points/grids/gpkg)
- [`tile_vector_grid()`](https://aavotins.github.io/egvtools/reference/tile_vector_grid.md)
  : Tile a vector grid into parquet tiles
- [`tiled_buffers()`](https://aavotins.github.io/egvtools/reference/tiled_buffers.md)
  : Create buffered tiles from point layers
- [`create_backgrounds()`](https://aavotins.github.io/egvtools/reference/create_backgrounds.md)
  : Create constant-background rasters from a directory of GeoTIFFs

## Internals

- [`egvtools`](https://aavotins.github.io/egvtools/reference/egvtools.md)
  [`egvtools-package`](https://aavotins.github.io/egvtools/reference/egvtools.md)
  : egvtools: High-resolution Ecogeographical Variable Workflows

## Datasets

Spatial reference datasets used by egvtools for examples, testing, and
spatial processing workflows.

- [`clc18`](https://aavotins.github.io/egvtools/reference/clc18.md) :
  Example Corine Land Cover (2018) polygons
- [`pts1000_sauzeme`](https://aavotins.github.io/egvtools/reference/pts1000_sauzeme.md)
  : Example 1000-m grid points
- [`pts100_sauzeme`](https://aavotins.github.io/egvtools/reference/pts100_sauzeme.md)
  : Example 100-m grid points
- [`pts300_sauzeme`](https://aavotins.github.io/egvtools/reference/pts300_sauzeme.md)
  : Example 300-m grid points
- [`pts500_sauzeme`](https://aavotins.github.io/egvtools/reference/pts500_sauzeme.md)
  : Example 500-m grid points
- [`tikls100_sauzeme`](https://aavotins.github.io/egvtools/reference/tikls100_sauzeme.md)
  : Example 100-m grid polygons
- [`tikls1km_sauzeme`](https://aavotins.github.io/egvtools/reference/tikls1km_sauzeme.md)
  : Example 1000-m grid polygons
- [`tikls300_sauzeme`](https://aavotins.github.io/egvtools/reference/tikls300_sauzeme.md)
  : Example 300-m grid polygons
- [`tikls500_sauzeme`](https://aavotins.github.io/egvtools/reference/tikls500_sauzeme.md)
  : Example 500-m grid points
- [`tks93_50km`](https://aavotins.github.io/egvtools/reference/tks93_50km.md)
  : Example 50-km map tiles
