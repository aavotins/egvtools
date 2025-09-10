
<!-- README.md is generated from README.Rmd. Please edit that file -->

# egvtools: High-resolution Ecogeographical Variable Workflows

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/aavotins/egvtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aavotins/egvtools/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/aavotins/egvtools/graph/badge.svg)](https://app.codecov.io/gh/aavotins/egvtools)
<!-- badges: end -->

`egvtools` provides a coherent set of wrappers and utilities that make
large-scale **EGV** creation reproducible and pleasant on real datasets.
The package leans on robust building blocks—`terra`, `sf`, `sfarrow`,
`exactextractr`, and `whitebox`—and standardizes I/O, naming
conventions, and multi-scale zonal statistics so your pipelines are
repeatable across machines and projects.

The package was developed to simply our work in project “HiQBioDiv:
High-resolution quantification of biodiversity for conservation and
management” funded by the Latvian Council of Science (Ref.
No. VPP-VARAM-DABA-2024/1-0002) and to ease reproduction of our work.
Five of the functions are strictly for replication, while others are
useful for wider audience (see documentation and
[articles](https://aavotins.github.io/egvtools/)).

## Installation

You can install the development version of egvtools from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("aavotins/egvtools")
```

## Usage

Functions in this package can be devided in two parts with extra
intermediate wrapper:

1)  helper functions to prepare analysis templates for reproduction:

- `download_raster_templates()` — fetch template rasters (from Zenodo
  repository);

- `download_vector_templates()` — fetch template vector grids/points
  (from Zenodo repository);

- `tile_vector_grid()` — tile template (vector) grid for chunked
  processing;

- `tiled_buffers()` — precompute buffered tiles for multiple radii;

2)  intermediate wrapper around `terra::ifel()`:

- `create_backgrounds()` — build consistent background rasters, guarding
  spatial cover, resolution, coordinate reference system, exact pixel
  matching, etc.;

3)  core analysis functions - small workflows, that are easily
    generalizable to other areas or usecases. Every function guards
    spatial cover, resolution, coordinate reference system, exact pixel
    matching, etc.:

- `polygon2input()` — rasterize polygons to ultra-high-resolution
  template, handle background/mask;

- `input2egv()` — normalize/align ultra-high-resolution inputs to
  broader-resolution EGV output rasters with guards to template;

- `downscale2egv()` — downscale coarse rasters to template grid and
  optionally smooth with IDW;

- `distance2egv()` — distances to features in inputs, summarised to EGV
  resolution with optional gap filling at the edges;

- `landscape_function()` — landscape-level per-zone metrics, tiled;

- `radius_function()` — multi-scale zonal statistics (dense/sparse).

In this package we use various geodata. Vector data need to be
polygonised before `polygon2input()`. Multiple outputs of this function
can be combined before creating EGVs.

We use term *input* for raster layers of higher resolution (exact
multiple) than EGV used for species distribution analysis. These layers
are for geodata harmonisation and standartisation in a *much* faster and
memory-friendly approach.

Every other function ending with `*egv()` and `landscape_*()` and
`radius_*()` functions create standartised and harmonised EGVs.

Functions ending with `*egv()` and `landscape_*()` function operate at
EGV cell resolution. While `radius_function()` creates output matching
EGV template with cell values representing aggregated information from
larger scales (at specified radius around every EGV-cells center)

## Code of Conduct

Please note that the egvtools project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
