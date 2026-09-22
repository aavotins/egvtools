skip_radius_test_deps <- function() {
  pkgs <- c(
    "terra", "sf", "raster", "dplyr", "fs", "sfarrow",
    "fasterize", "exactextractr", "future", "furrr", "glue"
  )

  for (pkg in pkgs) {
    testthat::skip_if_not_installed(pkg)
  }
}

write_sf_parquet_silently <- function(x, path) {
  suppressWarnings(
    sfarrow::st_write_parquet(obj = x, dsn = path)
  )
  invisible(path)
}

make_radius_argument_fixture <- function() {
  root <- tempfile("generalized-radius-args-")
  dir.create(root, recursive = TRUE)

  grid_dir <- file.path(root, "grid")
  buffer_dir <- file.path(root, "buffers")
  output_dir <- file.path(root, "output")
  dir.create(grid_dir)
  dir.create(buffer_dir)
  dir.create(output_dir)

  template_path <- file.path(root, "template.tif")
  input_path <- file.path(root, "input.tif")

  r <- terra::rast(
    xmin = 0, xmax = 1,
    ymin = 0, ymax = 1,
    resolution = 1,
    crs = "EPSG:3059"
  )
  terra::values(r) <- 1

  terra::writeRaster(r, template_path, overwrite = TRUE)
  terra::writeRaster(r, input_path, overwrite = TRUE)

  list(
    root = root,
    grid_dir = grid_dir,
    buffer_dir = buffer_dir,
    output_dir = output_dir,
    template_path = template_path,
    input_path = input_path
  )
}

make_radius_spatial_fixture <- function() {
  root <- tempfile("generalized-radius-spatial-")
  dir.create(root, recursive = TRUE)

  grid_dir <- file.path(root, "grid")
  buffer_dir <- file.path(root, "buffers")
  output_dir <- file.path(root, "output")
  dir.create(grid_dir)
  dir.create(buffer_dir)
  dir.create(output_dir)

  template_path <- file.path(root, "template.tif")
  input_path <- file.path(root, "xcoord.tif")
  reference_path <- file.path(root, "reference.parquet")
  reference_no_key_path <- file.path(root, "reference_no_key.parquet")

  template <- terra::rast(
    xmin = 0, xmax = 2,
    ymin = 0, ymax = 2,
    resolution = 1,
    crs = "EPSG:3059"
  )
  terra::values(template) <- 1
  terra::writeRaster(template, template_path, overwrite = TRUE)

  input <- template
  xy <- terra::xyFromCell(input, seq_len(terra::ncell(input)))
  terra::values(input) <- xy[, 1]
  terra::writeRaster(input, input_path, overwrite = TRUE)

  bbox <- sf::st_bbox(
    c(xmin = 0, ymin = 0, xmax = 2, ymax = 2),
    crs = sf::st_crs(3059)
  )

  geom <- sf::st_make_grid(
    sf::st_as_sfc(bbox),
    n = c(2, 2),
    what = "polygons"
  )

  grid <- sf::st_sf(
    id = seq_along(geom),
    geometry = geom
  )

  centroids <- suppressWarnings(sf::st_centroid(grid))
  cx <- sf::st_coordinates(centroids)[, 1]
  grid$parent_id <- ifelse(cx < 1, "L", "R")
  grid$rinda300 <- ifelse(cx < 1, 1L, 2L)
  grid$ID1km <- 1L

  # Base tiled rasterization grid. The final underscore-separated component
  # becomes the tile ID used to match buffer files.
  tile_path <- file.path(grid_dir, "grid_T1.parquet")
  write_sf_parquet_silently(grid, tile_path)

  # External rasterization grid. Repeated parent_id values are deliberate:
  # one buffer result should be distributed to multiple output cells.
  write_sf_parquet_silently(grid, reference_path)

  reference_no_key <- grid
  reference_no_key$parent_id <- NULL
  write_sf_parquet_silently(reference_no_key, reference_no_key_path)

  # Direct-to-tile buffer files. Their polygons are the output cells
  # themselves, so a mean extraction should reconstruct the input raster.
  write_sf_parquet_silently(
    grid[, c("id", "geometry")],
    file.path(buffer_dir, "custom_r750_T1.parquet")
  )

  write_sf_parquet_silently(
    grid[, c("id", "geometry")],
    file.path(buffer_dir, "pts100_r3000_T1.parquet")
  )

  write_sf_parquet_silently(
    grid[, c("id", "geometry")],
    file.path(buffer_dir, "pts100_r500_T1.parquet")
  )

  write_sf_parquet_silently(
    grid[, c("id", "geometry")],
    file.path(buffer_dir, "pts100_r1250_T1.parquet")
  )

  # Buffer file with a key that exists in the buffers but not in the base tile.
  alt_key <- grid[, "geometry", drop = FALSE]
  alt_key$alt_id <- seq_len(nrow(alt_key))
  alt_key <- alt_key[, c("alt_id", "geometry")]
  write_sf_parquet_silently(
    alt_key,
    file.path(buffer_dir, "alt_r800_T1.parquet")
  )

  # Duplicate buffer keys should be rejected because one extracted statistic
  # must correspond to one unique buffer key.
  duplicate_key <- grid[, c("id", "geometry")]
  duplicate_key$id[2] <- duplicate_key$id[1]
  write_sf_parquet_silently(
    duplicate_key,
    file.path(buffer_dir, "duplicate_r900_T1.parquet")
  )

  # Two coarse buffers: left and right halves of the test extent.
  parent_levels <- c("L", "R")
  coarse_geoms <- lapply(parent_levels, function(z) {
    sf::st_union(sf::st_geometry(grid[grid$parent_id == z, , drop = FALSE]))[[1]]
  })

  coarse <- sf::st_sf(
    parent_id = parent_levels,
    geometry = sf::st_sfc(coarse_geoms, crs = sf::st_crs(grid))
  )

  write_sf_parquet_silently(
    coarse,
    file.path(buffer_dir, "coarse_r2000_T1.parquet")
  )

  sparse_3000 <- coarse
  sparse_3000$rinda300 <- c(1L, 2L)
  sparse_3000$parent_id <- NULL
  sparse_3000 <- sparse_3000[, c("rinda300", "geometry")]
  write_sf_parquet_silently(
    sparse_3000,
    file.path(buffer_dir, "pts300_r3000_T1.parquet")
  )

  all_geom <- sf::st_union(sf::st_geometry(grid))[[1]]
  sparse_10000 <- sf::st_sf(
    ID1km = 1L,
    geometry = sf::st_sfc(all_geom, crs = sf::st_crs(grid))
  )
  write_sf_parquet_silently(
    sparse_10000,
    file.path(buffer_dir, "pts1000_r10000_T1.parquet")
  )

  list(
    root = root,
    grid_dir = grid_dir,
    buffer_dir = buffer_dir,
    output_dir = output_dir,
    template_path = template_path,
    input_path = input_path,
    reference_path = reference_path,
    reference_no_key_path = reference_no_key_path,
    expected = input
  )
}

run_radius_test <- function(fx, ..., output_subdir = "run") {
  output_dir <- file.path(fx$root, output_subdir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  generalized_radius_function(
    kvadrati_path = fx$grid_dir,
    radii_path = fx$buffer_dir,
    template_path = fx$template_path,
    input_layers = fx$input_path,
    layer_prefixes = "xcoord",
    output_dir = output_dir,
    n_workers = 1,
    fill_missing = FALSE,
    unlink_tiles = TRUE,
    quiet = TRUE,
    ...
  )
}

expect_raster_values_equal <- function(path, expected, tolerance = 1e-6) {
  testthat::expect_true(file.exists(path))

  observed <- terra::rast(path)

  observed_ext <- c(
    xmin = terra::xmin(observed),
    xmax = terra::xmax(observed),
    ymin = terra::ymin(observed),
    ymax = terra::ymax(observed)
  )

  expected_ext <- c(
    xmin = terra::xmin(expected),
    xmax = terra::xmax(expected),
    ymin = terra::ymin(expected),
    ymax = terra::ymax(expected)
  )

  testthat::expect_equal(
    observed_ext,
    expected_ext,
    tolerance = tolerance
  )

  testthat::expect_equal(
    terra::res(observed),
    terra::res(expected),
    tolerance = tolerance
  )

  testthat::expect_equal(
    terra::crs(observed, proj = TRUE),
    terra::crs(expected, proj = TRUE)
  )

  testthat::expect_equal(
    as.numeric(terra::values(observed)),
    as.numeric(terra::values(expected)),
    tolerance = tolerance
  )
}


testthat::test_that("specified mode validates buffer_spec", {
  skip_radius_test_deps()
  fx <- make_radius_argument_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  base_args <- list(
    kvadrati_path = fx$grid_dir,
    radii_path = fx$buffer_dir,
    template_path = fx$template_path,
    input_layers = fx$input_path,
    layer_prefixes = "x",
    output_dir = fx$output_dir,
    n_workers = 1,
    fill_missing = FALSE,
    quiet = TRUE,
    buffer_mode = "specified"
  )

  testthat::expect_error(
    do.call(generalized_radius_function, base_args),
    "`buffer_spec` must be a data.frame"
  )

  missing_column_spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750,
    key_field = "id"
  )

  testthat::expect_error(
    do.call(
      generalized_radius_function,
      c(base_args, list(buffer_spec = missing_column_spec))
    ),
    "missing required column.*grid_path"
  )

  noninteger_radius_spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750.5,
    key_field = "id",
    grid_path = NA_character_
  )

  testthat::expect_error(
    do.call(
      generalized_radius_function,
      c(base_args, list(buffer_spec = noninteger_radius_spec))
    ),
    "must be whole numbers"
  )

  duplicate_radius_spec <- data.frame(
    buffer_grid = c("a", "b"),
    radius_m = c(750, 750),
    key_field = c("id", "id"),
    grid_path = c(NA_character_, NA_character_)
  )

  testthat::expect_error(
    do.call(
      generalized_radius_function,
      c(base_args, list(buffer_spec = duplicate_radius_spec))
    ),
    "Each radius may occur only once"
  )

  nonexistent_grid_spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750,
    key_field = "id",
    grid_path = file.path(fx$root, "does-not-exist.parquet")
  )

  testthat::expect_error(
    do.call(
      generalized_radius_function,
      c(base_args, list(buffer_spec = nonexistent_grid_spec))
    ),
    "External rasterization grid not found"
  )
})


testthat::test_that("requested radii must be defined by the resolved specification", {
  skip_radius_test_deps()
  fx <- make_radius_argument_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750,
    key_field = "id",
    grid_path = NA_character_
  )

  testthat::expect_error(
    generalized_radius_function(
      kvadrati_path = fx$grid_dir,
      radii_path = fx$buffer_dir,
      template_path = fx$template_path,
      input_layers = fx$input_path,
      layer_prefixes = "x",
      output_dir = fx$output_dir,
      n_workers = 1,
      fill_missing = FALSE,
      quiet = TRUE,
      buffer_mode = "specified",
      buffer_spec = spec,
      radii = 1000
    ),
    "Requested radius/radii not defined"
  )
})


testthat::test_that("sparse mode requires a reference grid", {
  skip_radius_test_deps()
  fx <- make_radius_argument_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  testthat::expect_error(
    generalized_radius_function(
      kvadrati_path = fx$grid_dir,
      radii_path = fx$buffer_dir,
      template_path = fx$template_path,
      input_layers = fx$input_path,
      layer_prefixes = "x",
      output_dir = fx$output_dir,
      n_workers = 1,
      fill_missing = FALSE,
      quiet = TRUE,
      buffer_mode = "sparse"
    ),
    "`reference_grid_path` is required"
  )
})


testthat::test_that("sparse mode retains the standard grid-to-radius mapping", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  result <- generalized_radius_function(
    kvadrati_path = fx$grid_dir,
    radii_path = fx$buffer_dir,
    reference_grid_path = fx$reference_path,
    template_path = fx$template_path,
    input_layers = fx$input_path,
    layer_prefixes = "xcoord",
    output_dir = file.path(fx$root, "sparse"),
    n_workers = 1,
    radii = c(500, 3000),
    fill_missing = FALSE,
    unlink_tiles = TRUE,
    quiet = TRUE,
    buffer_mode = "sparse"
  )

  testthat::expect_equal(result$radius, c("r500", "r3000"))
  testthat::expect_equal(result$buffer_grid, c("pts100", "pts300"))

  expect_raster_values_equal(
    file.path(fx$root, "sparse", "xcoord_r500.tif"),
    fx$expected
  )
  expect_raster_values_equal(
    file.path(fx$root, "sparse", "xcoord_r3000.tif"),
    fx$expected
  )
})


testthat::test_that("specified mode processes all specification rows when radii is omitted", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = c("custom", "coarse"),
    radius_m = c(750, 2000),
    key_field = c("id", "parent_id"),
    grid_path = c(NA_character_, fx$reference_path),
    stringsAsFactors = FALSE
  )

  result <- run_radius_test(
    fx,
    buffer_mode = "specified",
    buffer_spec = spec,
    output_subdir = "specified-all"
  )

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_setequal(result$radius, c("r750", "r2000"))
  testthat::expect_setequal(result$buffer_grid, c("custom", "coarse"))
  testthat::expect_equal(result$n_tiles_merged, c(1L, 1L))
  testthat::expect_false(any(result$gap_filled))

  out_dir <- file.path(fx$root, "specified-all")
  expect_raster_values_equal(
    file.path(out_dir, "xcoord_r750.tif"),
    fx$expected
  )
  expect_raster_values_equal(
    file.path(out_dir, "xcoord_r2000.tif"),
    fx$expected
  )
})


testthat::test_that("dense 3000 m buffers join directly to the tile by id", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  result <- run_radius_test(
    fx,
    buffer_mode = "dense",
    radii = 3000,
    output_subdir = "dense-3000"
  )

  testthat::expect_equal(result$radius, "r3000")
  testthat::expect_equal(result$radius_m, 3000)
  testthat::expect_equal(result$buffer_grid, "pts100")
  testthat::expect_equal(result$n_tiles_merged, 1L)

  expect_raster_values_equal(
    file.path(fx$root, "dense-3000", "xcoord_r3000.tif"),
    fx$expected
  )
})


testthat::test_that("buffer key must exist in buffered polygons", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750,
    key_field = "missing_key",
    grid_path = NA_character_
  )

  testthat::expect_error(
    run_radius_test(
      fx,
      buffer_mode = "specified",
      buffer_spec = spec,
      output_subdir = "missing-buffer-key"
    ),
    "Key field 'missing_key' is missing from buffered polygons"
  )
})


testthat::test_that("buffer key must uniquely identify buffered polygons", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = "duplicate",
    radius_m = 900,
    key_field = "id",
    grid_path = NA_character_
  )

  testthat::expect_error(
    run_radius_test(
      fx,
      buffer_mode = "specified",
      buffer_spec = spec,
      output_subdir = "duplicate-buffer-key"
    ),
    "does not uniquely identify buffered polygons"
  )
})


testthat::test_that("tile rasterization grid must contain the specified key", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = "alt",
    radius_m = 800,
    key_field = "alt_id",
    grid_path = NA_character_
  )

  testthat::expect_error(
    run_radius_test(
      fx,
      buffer_mode = "specified",
      buffer_spec = spec,
      output_subdir = "missing-tile-key"
    ),
    "Key field 'alt_id' is missing from the rasterization grid"
  )
})


testthat::test_that("external rasterization grid must contain the specified key", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  spec <- data.frame(
    buffer_grid = "coarse",
    radius_m = 2000,
    key_field = "parent_id",
    grid_path = fx$reference_no_key_path
  )

  testthat::expect_error(
    run_radius_test(
      fx,
      buffer_mode = "specified",
      buffer_spec = spec,
      output_subdir = "missing-reference-key"
    ),
    "Key field 'parent_id' is missing from reference grid"
  )
})


testthat::test_that("tile filename suffix is used as the tile identifier", {
  skip_radius_test_deps()
  fx <- make_radius_spatial_fixture()
  on.exit(unlink(fx$root, recursive = TRUE, force = TRUE), add = TRUE)

  # The fixture uses grid_T1.parquet and custom_r750_T1.parquet. Successful
  # output therefore verifies matching by the final underscore-separated tile
  # identifier rather than by a fixed character position.
  spec <- data.frame(
    buffer_grid = "custom",
    radius_m = 750,
    key_field = "id",
    grid_path = "tile"
  )

  result <- run_radius_test(
    fx,
    buffer_mode = "specified",
    buffer_spec = spec,
    output_subdir = "tile-id"
  )

  testthat::expect_equal(result$n_tiles_merged, 1L)
  testthat::expect_true(
    file.exists(file.path(fx$root, "tile-id", "xcoord_r750.tif"))
  )
})
