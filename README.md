
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
useful for wider audience (see documentation and [articles]()).

## Installation

You can install the development version of egvtools from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("aavotins/egvtools")
```

## Example

There will be example here. Soon.

``` r
library(egvtools)
## basic example code
```

## Code of Conduct

Please note that the egvtools project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
