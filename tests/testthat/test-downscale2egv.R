test_that("downscale2egv validates output names", {
  skip_if_spatial_packages_missing(c("terra", "sf", "sfarrow", "whitebox", "fs"))

  r <- make_test_raster()
  grid <- make_rect_sf(0, 4, 0, 4)
  out <- tempfile("dwn-")

  expect_error(
    downscale2egv(r, grid, r, out_path = out, layer_name = "x", fill_gaps = FALSE, quiet = TRUE),
    "file_name"
  )
})

test_that("downscale2egv produces a template-aligned raster without Whitebox filling", {
  skip_if_spatial_packages_missing(c("terra", "sf", "sfarrow", "whitebox", "fs"))

  template <- make_test_raster(
    nrows = 4, ncols = 4,
    xmin = 0, xmax = 4, ymin = 0, ymax = 4,
    values = rep(1, 16)
  )
  raw <- make_test_raster(
    nrows = 8, ncols = 8,
    xmin = 0, xmax = 4, ymin = 0, ymax = 4,
    values = seq_len(64)
  )
  grid <- make_rect_sf(0, 4, 0, 4)
  out_dir <- tempfile("dwn-")

  ans <- downscale2egv(
    template_path = template,
    grid_path = grid,
    rawfile_path = raw,
    out_path = out_dir,
    file_name = "downscaled.tif",
    layer_name = "downscaled",
    buffer_m = 0,
    check_na = TRUE,
    fill_gaps = FALSE,
    smooth = FALSE,
    return_visible = TRUE,
    quiet = TRUE
  )

  expect_s3_class(ans, "data.frame")
  expect_equal(nrow(ans), 1L)
  expect_true(file.exists(ans$output))
  expect_equal(ans$gap_count, 0)
  expect_false(ans$smoothed)

  out <- terra::rast(ans$output)
  expect_same_grid(out, template)
  expect_equal(names(out), "downscaled")
})

test_that("downscale2egv rejects unsupported input object types", {
  skip_if_spatial_packages_missing(c("terra", "sf", "sfarrow", "whitebox", "fs"))

  grid <- make_rect_sf(0, 4, 0, 4)
  out <- tempfile("dwn-errors-")

  expect_error(
    downscale2egv(
      template_path = list(a = 1),
      grid_path = grid,
      rawfile_path = make_test_raster(),
      out_path = out,
      file_name = "x.tif",
      layer_name = "x",
      fill_gaps = FALSE,
      quiet = TRUE
    ),
    "template raster must be"
  )
})
