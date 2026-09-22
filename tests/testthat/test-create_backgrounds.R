test_that("create_backgrounds validates input and handles an empty directory", {
  skip_if_spatial_packages_missing(c("terra", "fs"))

  expect_error(
    create_backgrounds(file.path(tempdir(), "definitely-not-an-egvtools-dir"), quiet = TRUE),
    "existing directory"
  )

  in_dir <- tempfile("cbg-empty-")
  dir.create(in_dir)
  out_dir <- tempfile("cbg-out-")

  ans <- create_backgrounds(in_dir, out_dir = out_dir, quiet = TRUE)

  expect_s3_class(ans, "data.frame")
  expect_named(ans, c("in_file", "out_file", "n_cells", "elapsed_sec"))
  expect_equal(nrow(ans), 0L)
  expect_true(dir.exists(out_dir))
})

test_that("create_backgrounds replaces non-NA cells and preserves NA cells", {
  skip_if_spatial_packages_missing(c("terra", "fs"))

  in_dir <- tempfile("cbg-in-")
  out_dir <- tempfile("cbg-out-")
  dir.create(in_dir)

  r <- make_test_raster(nrows = 2, ncols = 3, values = c(1, NA, 3, 4, NA, 6))
  write_test_raster(r, file.path(in_dir, "source.tif"))

  ans <- create_backgrounds(
    in_dir = in_dir,
    out_dir = out_dir,
    background_value = 0,
    quiet = TRUE
  )

  expect_equal(nrow(ans), 1L)
  expect_equal(basename(ans$out_file), "nulls_source.tif")
  expect_equal(ans$n_cells, 6L)
  expect_true(file.exists(ans$out_file))

  out <- terra::rast(ans$out_file)
  expect_equal(
    terra::values(out, mat = FALSE),
    c(0, NA, 0, 0, NA, 0)
  )
  expect_equal(terra::crs(out), terra::crs(r))
})

test_that("create_backgrounds uses the non-zero prefix and respects overwrite", {
  skip_if_spatial_packages_missing(c("terra", "fs"))

  in_dir <- tempfile("cbg-in-")
  out_dir <- tempfile("cbg-out-")
  dir.create(in_dir)

  r <- make_test_raster(nrows = 2, ncols = 2, values = c(1, 2, NA, 4))
  write_test_raster(r, file.path(in_dir, "a.tiff"))

  first <- create_backgrounds(
    in_dir,
    out_dir = out_dir,
    background_value = 0.5,
    quiet = TRUE
  )
  expect_equal(basename(first$out_file), "bg0.5_a.tif")

  mtime <- file.info(first$out_file)$mtime
  second <- create_backgrounds(
    in_dir,
    out_dir = out_dir,
    background_value = 0.5,
    overwrite = FALSE,
    quiet = TRUE
  )

  expect_false(isTRUE(second$n_cells[1] > 0))
  expect_true(is.na(second$n_cells[1]))
  expect_equal(file.info(first$out_file)$mtime, mtime)
})

test_that("create_backgrounds filters by recursion and pattern", {
  skip_if_spatial_packages_missing(c("terra", "fs"))

  in_dir <- tempfile("cbg-filter-")
  out_dir <- tempfile("cbg-out-")
  dir.create(in_dir)
  dir.create(file.path(in_dir, "sub"))

  r <- make_test_raster(nrows = 2, ncols = 2, values = 1:4)
  write_test_raster(r, file.path(in_dir, "keep_main.tif"))
  write_test_raster(r, file.path(in_dir, "drop_main.tif"))
  write_test_raster(r, file.path(in_dir, "sub", "keep_nested.tif"))
  writeLines("not a raster", file.path(in_dir, "keep.txt"))

  ans <- create_backgrounds(
    in_dir,
    out_dir = out_dir,
    pattern = "^keep",
    recursive = TRUE,
    quiet = TRUE
  )

  expect_equal(nrow(ans), 2L)
  expect_setequal(
    basename(ans$out_file),
    c("nulls_keep_main.tif", "nulls_keep_nested.tif")
  )
})
