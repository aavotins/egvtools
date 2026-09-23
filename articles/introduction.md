# Introduction to egvtools

``` r

library(egvtools)
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(terra)
#> terra 1.9.50
library(sfarrow)
```

## Overview

`egvtools` provides a coherent set of wrappers and utilities that make
large-scale EGV creation reproducible and pleasant on real datasets. The
package leans on robust building blocks — `terra`, `sf`, `sfarrow`,
`exactextractr`, and `whitebox` — and standardizes I/O, naming
conventions, and multi-scale zonal statistics so your pipelines are
repeatable across machines and projects.

The package was developed in project “HiQBioDiv: High-resolution
quantification of biodiversity for conservation and management” funded
by the Latvian Council of Science (Ref. No. VPP-VARAM-DABA-2024/1-0002).

The development version can be installed from GitHub with:

``` r

# install.packages("remotes")
remotes::install_github("aavotins/egvtools")
```

## Terminology

Athough all georeferenced data can be considered geodata, in this
material we use the following terms in the order listed below in our
workflows:

- **raw geodata** - considered as raw data obtained for a harmonised
  description of the environment. This may include tables with
  coordinates, raster or vector data. It can be anything that has been
  or can be used to create ecogeographical variables, with or without
  slight processing.

- **geodata product** - processed raw geodata that have undegone heavy
  modifications, e.g. spatial overlays and combinations of different
  sets of raw geodata, and are used as input data. In this document,
  geodata products are categorical raster layers that match the CRS and
  the pixel locations of input data. When split by categories, they
  become input data. The processing step of creating geodata products is
  necessary when decisions about the order of spatial overlays are
  important. For example, in a high-resolution pixel, there can only be
  water or forest, if the edge between water and forest need to be
  calculated.

- **input data or input layers** - very-high resolution (multiple times
  higher than that used for ecogeographical variables) raster data that
  are the direct input for the creation of most of the ecogeographical
  variables. The creation of such layers is particularly useful
  alongside geodata products, as dealing with border misalignment or
  decisions regarding the order of spatial o verlays, as well as simple
  geoprocessing, is much faster with raster data.

- **ecogeographical variables (EGVs)** - this is the final product of
  the workflow describing environment for statistical analysis
  (e.g. species distribution modelling). They are suitable also for
  publishing due to standadisation of the values. In other words, these
  are standardised landscape ecological variables in the form of
  high-resolution raster layers.

## Example data

Several example datasets are included in the package. The examples in
this vignette use harmonized raster and vector grids, grid centroids as
point layers, climate and Land Use Land Cover (LULC) data represented by
(Corine Land Cover (CLC) for the year
2018)\[<https://land.copernicus.eu/en/products/corine-land-cover/clc2018>\].

``` r

data("aves_small", package = "sdmhelpers")
data("bryophyta_one", package = "sdmhelpers")
data("habitat_inventories", package = "sdmhelpers")
data("points100", package = "sdmhelpers")
data("example_grid", package = "sdmhelpers")
example_grid=terra::unwrap(example_grid)
data("example_presences", package = "sdmhelpers")
data("example_background", package = "sdmhelpers")
data("example_independent_blocks", package = "sdmhelpers")
data("egv_names", package = "sdmhelpers")
```

The raster data are stored in `inst/extdata`.

``` r

egv_file <- system.file("extdata",
                        "example_egv.tif",
                        package = "sdmhelpers")

egv_file
```

Harmonization data can also be downloaded from Zenodo repository:

- (vector data)\[<https://doi.org/10.5281/zenodo.22912420>\];

- (raster data)\[<https://doi.org/10.5281/zenodo.22912629>\].
