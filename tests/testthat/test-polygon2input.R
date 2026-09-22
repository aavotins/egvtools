test_that("polygon2input validates its principal inputs", {
  skip_if_spatial_packages_missing(c("terra", "sf", "raster", "fasterize", "fs"))

  expect_error(
    polygon2input(
      vector_data = data.frame(x = 1),
      template_path = "does-not-exist.tif",
      file_name = "x.tif",
      quiet = TRUE
    ),
    "sf object"
  )
})

test_that("polygon2input rasterizes polygons, covers gaps, and preserves template NA mask", {
  skip_if_spatial_packages_missing(c("terra", "sf", "raster", "fasterize", "fs"))

  template <- make_test_raster(
    nrows = 4, ncols = 4,
    xmin = 0, xmax = 4, ymin = 0, ymax = 4,
    values = c(rep(1, 15), NA)
  )
  template_file <- tempfile(fileext = ".tif")
  write_test_raster(template, template_file)

  poly <- make_rect_sf(0, 2, 0, 4, value = 7)
  out_dir <- tempfile("p2i-")

  ans <- polygon2input(
    vector_data = poly,
    template_path = template_file,
    out_path = out_dir,
    file_name = "poly.tif",
    value_field = "value",
    background_value = 0,
    check_na = TRUE,
    overwrite = TRUE,
    quiet = TRUE
  )

  expect_true(file.exists(ans$out_file))
  expect_equal(ans$n_cells, 16L)
  expect_equal(ans$n_na_final, 0L)

  out <- terra::rast(ans$out_file)
  expect_same_grid(out, template)

  left_cell <- terra::cellFromXY(out, matrix(c(0.5, 3.5), ncol = 2))
  right_cell <- terra::cellFromXY(out, matrix(c(3.5, 3.5), ncol = 2))
  masked_cell <- 16L

  expect_equal(unname(terra::values(out)[left_cell, 1]), 7)
  expect_equal(unname(terra::values(out)[right_cell, 1]), 0)
  expect_true(is.na(terra::values(out)[masked_cell, 1]))
})

test_that("polygon2input skips an existing output unless overwrite is requested", {
  skip_if_spatial_packages_missing(c("terra", "sf", "raster", "fasterize", "fs"))

  template <- make_test_raster(values = rep(1, 16))
  template_file <- tempfile(fileext = ".tif")
  write_test_raster(template, template_file)
  poly <- make_rect_sf(0, 2, 0, 4)
  out_dir <- tempfile("p2i-existing-")
  dir.create(out_dir)
  existing <- file.path(out_dir, "x.tif")
  write_test_raster(template, existing)

  ans <- polygon2input(
    poly,
    template_file,
    out_path = out_dir,
    file_name = "x.tif",
    overwrite = FALSE,
    quiet = TRUE
  )

  expect_equal(ans$out_file, existing)
  expect_true(is.na(ans$n_cells))
  expect_true(is.na(ans$n_na_initial))
  expect_true(is.na(ans$n_na_final))
})
