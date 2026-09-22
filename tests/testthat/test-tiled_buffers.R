test_that("tiled_buffers creates constant-radius buffers per tile", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs", "future", "furrr"))

  in_dir <- tempfile("buf-in-")
  out_dir <- tempfile("buf-out-")
  dir.create(in_dir)

  pts <- make_test_points_sf()
  src <- file.path(in_dir, "pts100_sauzeme.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(pts, src))

  ans <- suppress_sfarrow_geo_warning(tiled_buffers(
    in_dir = in_dir,
    out_dir = out_dir,
    buffer_mode = "dense",
    radii_dense = c(10, 20),
    split_field = "tks50km",
    n_workers = 1,
    quiet = TRUE
  ))

  expect_equal(nrow(ans), 4L)
  expect_true(all(ans$wrote))
  expect_setequal(ans$tileid, c("A", "B"))
  expect_setequal(ans$radius_m, c(10, 20))
  expect_true(all(file.exists(ans$out_file)))

  one <- sfarrow::st_read_parquet(ans$out_file[1])
  expect_true(all(sf::st_geometry_type(one) %in% c("POLYGON", "MULTIPOLYGON")))
  expect_equal(nrow(one), 2L)

  again <- suppress_sfarrow_geo_warning(tiled_buffers(
    in_dir = in_dir,
    out_dir = out_dir,
    buffer_mode = "dense",
    radii_dense = c(10, 20),
    split_field = "tks50km",
    n_workers = 1,
    overwrite = FALSE,
    quiet = TRUE
  ))
  expect_false(any(again$wrote))
})

test_that("tiled_buffers supports per-feature radii in specified mode", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs", "future", "furrr"))

  in_dir <- tempfile("buf-in-")
  out_dir <- tempfile("buf-out-")
  dir.create(in_dir)
  pts <- make_test_points_sf()
  src <- file.path(in_dir, "custom_points.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(pts, src))

  ans <- suppress_sfarrow_geo_warning(tiled_buffers(
    out_dir = out_dir,
    buffer_mode = "specified",
    points_path = src,
    radius_field = "radius_m",
    split_field = "tks50km",
    n_workers = 1,
    quiet = TRUE
  ))

  expect_equal(nrow(ans), 2L)
  expect_true(all(ans$mode == "field"))
  expect_true(all(ans$radius_field == "radius_m"))
  expect_true(all(file.exists(ans$out_file)))
})

test_that("tiled_buffers rejects unknown modes", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs", "future", "furrr"))

  expect_error(
    tiled_buffers(
      out_dir = tempfile("buf-out-"),
      buffer_mode = "other",
      n_workers = 1,
      quiet = TRUE
    ),
    "buffer_mode must be one of"
  )
})
