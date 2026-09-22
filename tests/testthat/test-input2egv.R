test_that("input2egv maps mean to average and aligns to the template", {
  skip_if_spatial_packages_missing(c("terra", "fs", "whitebox"))

  input <- make_test_raster(
    nrows = 4, ncols = 4,
    xmin = 0, xmax = 4, ymin = 0, ymax = 4,
    values = 1:16
  )
  tmpl <- make_test_raster(
    nrows = 2, ncols = 2,
    xmin = 0, xmax = 4, ymin = 0, ymax = 4,
    values = rep(1, 4)
  )
  out_dir <- tempfile("i2e-")

  ans <- input2egv(
    input = input,
    egv_template = tmpl,
    summary_function = "mean",
    missing_job = "none",
    check_alignment = FALSE,
    outlocation = out_dir,
    outfilename = "aligned.tif",
    layername = "aligned",
    return_visible = TRUE,
    quiet = TRUE
  )

  expect_equal(ans$method, "average")
  expect_equal(ans$missing_job, "none")
  expect_true(file.exists(ans$path))

  out <- terra::rast(ans$path)
  expect_same_grid(out, tmpl)
  expect_equal(names(out), "aligned")
  expect_equal(
    terra::global(out, "mean", na.rm = TRUE)[[1]],
    mean(1:16),
    tolerance = 1e-6
  )
})

test_that("input2egv forces nearest-neighbour for categorical data", {
  skip_if_spatial_packages_missing(c("terra", "fs", "whitebox"))

  input <- make_test_raster(values = rep(c(1, 2), 8))
  tmpl <- make_test_raster(nrows = 2, ncols = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4, values = 1)
  out_dir <- tempfile("i2e-cat-")

  ans <- input2egv(
    input,
    tmpl,
    summary_function = "average",
    missing_job = "none",
    is_categorical = TRUE,
    check_alignment = FALSE,
    outlocation = out_dir,
    outfilename = "cat.tif",
    layername = "cat",
    return_visible = TRUE,
    quiet = TRUE
  )

  expect_equal(ans$method, "near")
})

test_that("input2egv CoverOutput fills remaining template gaps", {
  skip_if_spatial_packages_missing(c("terra", "fs", "whitebox"))

  input <- make_test_raster(nrows = 2, ncols = 2, values = c(1, NA, 3, NA))
  tmpl <- make_test_raster(nrows = 2, ncols = 2, values = 1)
  out_dir <- tempfile("i2e-cover-")

  ans <- input2egv(
    input,
    tmpl,
    missing_job = "CoverOutput",
    output_bg = 9,
    check_alignment = FALSE,
    outlocation = out_dir,
    outfilename = "covered.tif",
    layername = "covered",
    return_visible = TRUE,
    quiet = TRUE
  )

  expect_equal(ans$n_gaps_initial, 2L)
  expect_equal(ans$n_gaps_final, 0L)
  expect_equal(
    terra::values(terra::rast(ans$path), mat = FALSE),
    c(1, 9, 3, 9)
  )
})

test_that("input2egv validates missing-value workflows", {
  skip_if_spatial_packages_missing(c("terra", "fs", "whitebox"))

  r <- make_test_raster()
  out_dir <- tempfile("i2e-errors-")

  expect_error(
    input2egv(
      r, r,
      missing_job = "CoverInput",
      outlocation = out_dir,
      outfilename = "x.tif",
      layername = "x",
      quiet = TRUE
    ),
    "input_bg"
  )

  expect_error(
    input2egv(
      r, r,
      missing_job = "CoverInput",
      input_bg = 0,
      outlocation = out_dir,
      outfilename = "x.tif",
      layername = "x",
      quiet = TRUE
    ),
    "template is missing"
  )

  expect_error(
    input2egv(
      r, r,
      missing_job = "not-a-job",
      outlocation = out_dir,
      outfilename = "x.tif",
      layername = "x",
      quiet = TRUE
    )
  )
})
