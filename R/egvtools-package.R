#' egvtools: High-resolution Ecogeographical Variable Workflows
#'
#' `egvtools` provides a coherent set of wrappers and utilities that make
#' large-scale **EGV** creation reproducible and pleasant on real datasets.
#' The package leans on robust building blocks—`terra`, `sf`, `sfarrow`,
#' `exactextractr`, and `whitebox`—and standardizes I/O, naming conventions,
#' and multi-scale zonal statistics so your pipelines are repeatable across machines
#' and projects.
#'
#' @section Core workflow (generalizable):
#' Functions that form the backbone of your analyses (thin wrappers over
#' existing libs, but opinionated for stability, speed and consistency):
#' - `polygon2input()` — rasterize polygons to template, handle background/mask.
#' - `downscale2egv()` — downscale coarse rasters to template grid and optionally smooth with IDW.
#' - `distance2egv()` — distances to features with optional gap filling at the edges.
#' - `input2egv()` — normalize/align inputs to EGV outputs with guards.
#' - `landscape_function()` — \pkg{landscapemetrics} landscape-level per-zone metrics, tiled.
#' - `radius_function()` — multi-scale zonal statistics (dense/sparse).
#'
#' @section Reproducibility helpers:
#' Utilities that set up inputs and structure so our (project HiQBioDiv) results
#' can be reproduced:
#' - `download_raster_templates()` — fetch template rasters (Zenodo) to
#'   canonical paths.
#' - `download_vector_templates()` — fetch template vector grids/points.
#' - `tile_vector_grid()` — tile template grids for chunked processing.
#' - `tiled_buffers()` — precompute buffered tiles for multiple radii.
#' - `create_backgrounds()` — build consistent background rasters/values.
#'
#' @section Design principles:
#' - **Tiled, RAM-aware I/O**; workers open data by path (avoid big globals) and
#' caches them if necessary (avoid big I/O).
#' - **Deterministic filenames** and strict layer/radius naming.
#' - **Guards** for empty/invalid geometries, coordinate reference systems, naming and all-NA joins.
#' - **Cross-platform parallel** via `{future}/{furrr}` when enabled.
#'
#' @section Options:
#' Package defaults are set on load (see `zzz-options.R`). Users may override:
#' - `options(egvtools.future_plan = "sequential")`
#' - `options(egvtools.progress    = TRUE)`
#'
#' @section Getting started:
#' 1. Run `download_*_templates()` to fetch canonical inputs.
#' 2. Use `tile_vector_grid()` / `tiled_buffers()` for scalable chunks.
#' 3. Produce EGVs at site scale with `polygon2input()`, `downscale2egv()`, `distance2egv()`,
#'    `input2egv()`, and summarise from larger scales with `landscape_function()` / `radius_function()`.
#'
#' @seealso
#'   Useful links:
#'   \itemize{
#'     \item GitHub repo: \url{https://github.com/aavotins/egvtools}
#'     \item Package site: \url{https://aavotins.github.io/egvtools/}
#'     \item Issues / bugs: \url{https://github.com/aavotins/egvtools/issues}
#'   }
#'
#' @docType package
#' @name egvtools
#' @aliases egvtools-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
