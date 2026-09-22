test_that("distance2egv requires output filename and layer name", {
  skip_if_spatial_packages_missing(c("terra", "whitebox"))

  r <- make_test_raster()
  tmpl <- make_test_raster(nrows = 2, ncols = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4, values = 1)

  expect_error(
    distance2egv(r, tmpl, use_whitebox = FALSE, quiet = TRUE),
    "outfilename"
  )
  expect_error(
    distance2egv(r, tmpl, use_whitebox = FALSE, outfilename = "x.tif", quiet = TRUE),
    "layername"
  )
})

test_that("distance2egv handles numeric and interval source classes", {
  skip_if_spatial_packages_missing(c("terra", "whitebox"))

  vals <- rep(0, 36)
  vals[c(1, 10, 36)] <- c(5, 15, 20)
  r <- make_test_raster(
    nrows = 6, ncols = 6,
    xmin = 0, xmax = 60, ymin = 0, ymax = 60,
    values = vals
  )
  tmpl <- make_test_raster(
    nrows = 3, ncols = 3,
    xmin = 0, xmax = 60, ymin = 0, ymax = 60,
    values = rep(1, 9)
  )

  out_dir <- tempfile("d2e-")
  dir.create(out_dir)

  ans <- distance2egv(
    input = r,
    template_egv = tmpl,
    values_as_one = c("5", "[15,20]"),
    use_whitebox = FALSE,
    fill_gaps = FALSE,
    outlocation = out_dir,
    outfilename = "dist.tif",
    layername = "dist_test",
    check_na = TRUE,
    quiet = TRUE
  )

  expect_equal(ans$n_sources, 3L)
  expect_true(file.exists(ans$path))
  expect_true(is.finite(ans$min_dist))
  expect_true(is.finite(ans$max_dist))
  expect_gte(ans$max_dist, ans$min_dist)

  out <- terra::rast(ans$path)
  expect_same_grid(out, tmpl)
  expect_equal(names(out), "dist_test")
  expect_equal(ans$n_na_final, 0L)
})

test_that("distance2egv requires template_input when initial projection is requested", {
  skip_if_spatial_packages_missing(c("terra", "whitebox"))

  r <- make_test_raster()
  tmpl <- make_test_raster()
  out_dir <- tempfile("d2e-")
  dir.create(out_dir)

  expect_error(
    distance2egv(
      r, tmpl,
      project_to_template_input = TRUE,
      use_whitebox = FALSE,
      outlocation = out_dir,
      outfilename = "x.tif",
      layername = "x",
      quiet = TRUE
    ),
    "template_input"
  )
})
