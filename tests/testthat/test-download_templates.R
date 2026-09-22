test_that("download_raster_templates is idempotent when files already exist", {
  skip_if_spatial_packages_missing(c("fs", "curl"))

  out_dir <- tempfile("raster-template-")
  dir.create(out_dir)
  writeLines("already here", file.path(out_dir, "marker.txt"))

  ans <- download_raster_templates(out_dir = out_dir, overwrite = FALSE, quiet = TRUE)

  expect_equal(ans$out_dir, out_dir)
  expect_length(ans$files_written, 0L)
})

test_that("download_raster_templates unpacks mocked archive contents", {
  skip_if_spatial_packages_missing(c("fs", "curl"))
  skip_if_not_installed("testthat", minimum_version = "3.1.0")

  local_mocked_bindings(
    curl_download = function(url, destfile, mode = "wb", quiet = FALSE, ...) {
      writeBin(as.raw(c(1, 2, 3)), destfile)
      destfile
    },
    .package = "curl"
  )
  local_mocked_bindings(
    unzip = function(zipfile, exdir, ...) {
      dir.create(file.path(exdir, "archive"), recursive = TRUE)
      writeLines("a", file.path(exdir, "archive", "A.tif"))
      writeLines("b", file.path(exdir, "B.tif"))
      invisible(c("archive/A.tif", "B.tif"))
    },
    .package = "utils"
  )

  out_dir <- tempfile("raster-template-")
  ans <- download_raster_templates(
    url = "https://example.invalid/archive.zip",
    out_dir = out_dir,
    overwrite = TRUE,
    quiet = TRUE
  )

  expect_true(file.exists(file.path(out_dir, "A.tif")))
  expect_true(file.exists(file.path(out_dir, "B.tif")))
  expect_setequal(basename(ans$files_written), c("A.tif", "B.tif"))
})

test_that("download_vector_templates classifies mocked archive contents", {
  skip_if_spatial_packages_missing(c("fs", "curl"))
  skip_if_not_installed("testthat", minimum_version = "3.1.0")

  local_mocked_bindings(
    curl_download = function(url, destfile, mode = "wb", quiet = FALSE, ...) {
      writeBin(as.raw(c(1, 2, 3)), destfile)
      destfile
    },
    .package = "curl"
  )
  local_mocked_bindings(
    unzip = function(zipfile, exdir, ...) {
      writeLines("grid", file.path(exdir, "tikls100_sauzeme.parquet"))
      writeLines("points", file.path(exdir, "pts100_sauzeme.parquet"))
      writeLines("gpkg", file.path(exdir, "vector_grids.gpkg"))
      writeLines("ignore", file.path(exdir, "README.txt"))
      invisible(character())
    },
    .package = "utils"
  )

  grid_dir <- tempfile("grid-")
  points_dir <- tempfile("points-")
  gpkg_dir <- tempfile("gpkg-")

  ans <- download_vector_templates(
    url = "https://example.invalid/archive.zip",
    grid_dir = grid_dir,
    points_dir = points_dir,
    gpkg_dir = gpkg_dir,
    overwrite = TRUE,
    quiet = TRUE
  )

  expect_equal(ans$grid_dir, grid_dir)
  expect_equal(ans$points_dir, points_dir)
  expect_equal(ans$gpkg_dir, gpkg_dir)
  expect_true(file.exists(file.path(grid_dir, "tikls100_sauzeme.parquet")))
  expect_true(file.exists(file.path(points_dir, "pts100_sauzeme.parquet")))
  expect_true(file.exists(file.path(gpkg_dir, "vector_grids.gpkg")))
  expect_false(file.exists(file.path(grid_dir, "README.txt")))
})
