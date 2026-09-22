test_that("landscape_function validates rasterization engine before doing work", {
  expect_error(
    landscape_function(
      landscape = NULL,
      out_dir = tempdir(),
      out_filename = "x.tif",
      out_layername = "x",
      rasterize_engine = "invalid"
    )
  )
})

test_that("landscape_function validates zone identifiers", {
  skip_if_spatial_packages_missing(c(
    "terra", "sf", "raster", "landscapemetrics", "stringr", "dplyr",
    "tibble", "purrr", "fs", "sfarrow", "fasterize", "future", "furrr", "whitebox"
  ))

  r <- make_test_raster(values = rep(c(1, 2), 8))
  zones <- make_rect_sf(0, 4, 0, 4, tile = "A")

  expect_error(
    landscape_function(
      landscape = r,
      zones = zones,
      id_field = "missing_id",
      tile_field = "tile",
      template = r,
      out_dir = tempfile("landscape-") ,
      out_filename = "x.tif",
      out_layername = "x",
      n_workers = 1,
      quiet = TRUE
    ),
    "id_field 'missing_id'"
  )
})

test_that("landscape_function can complete a tiny deterministic mocked metric run", {
  skip_if_spatial_packages_missing(c(
    "terra", "sf", "raster", "landscapemetrics", "stringr", "dplyr",
    "tibble", "purrr", "fs", "sfarrow", "fasterize", "future", "furrr", "whitebox", "sp"
  ))
  skip_if_not_installed("testthat", minimum_version = "3.1.0")

  local_mocked_bindings(
    sample_lsm = function(landscape, y, plot_id, what, ...) {
      data.frame(
        layer = 1L,
        plot_id = as.character(plot_id),
        value = seq_along(plot_id),
        stringsAsFactors = FALSE
      )
    },
    .package = "landscapemetrics"
  )

  r <- make_test_raster(values = rep(c(1, 2), 8))
  zones <- make_two_zone_sf()
  out_dir <- tempfile("landscape-")

  ans <- landscape_function(
    landscape = r,
    zones = zones,
    id_field = "id",
    tile_field = "tile",
    template = r,
    out_dir = out_dir,
    out_filename = "metric.tif",
    out_layername = "metric",
    rasterize_engine = "terra",
    n_workers = 1,
    skip_existing = FALSE,
    report_gaps = TRUE,
    report_gap_size = FALSE,
    fill_gaps = FALSE,
    quiet = TRUE
  )

  expect_s3_class(ans, "data.frame")
  expect_equal(nrow(ans), 1L)
  expect_equal(ans$n_tiles, 1L)
  expect_equal(ans$n_zones, 2L)
  expect_equal(ans$n_layers, 1L)
  expect_true(file.exists(ans$output_path))

  out <- terra::rast(ans$output_path)
  expect_same_grid(out, r)
  expect_equal(names(out), "metric")
})
