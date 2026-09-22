test_that("tile_vector_grid splits by an existing tile field", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs"))

  grid <- sf::st_as_sf(
    data.frame(
      x = 1:6,
      y = 1:6,
      lapa = c("A", "A", "A", "B", "B", "B")
    ),
    coords = c("x", "y"),
    crs = 3059
  )

  src_dir <- tempfile("tile-src-")
  out_dir <- tempfile("tile-out-")
  dir.create(src_dir)
  grid_path <- file.path(src_dir, "tikls100_sauzeme.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(grid, grid_path))

  ans <- suppress_sfarrow_geo_warning(tile_vector_grid(
    grid_path = grid_path,
    out_dir = out_dir,
    tile_field = "lapa",
    quiet = TRUE
  ))

  expect_equal(nrow(ans), 2L)
  expect_setequal(ans$tile_id, c("A", "B"))
  expect_true(all(ans$n_rows == 3L))
  expect_true(all(ans$wrote))
  expect_setequal(
    basename(ans$path),
    c("tikls100_A.parquet", "tikls100_B.parquet")
  )

  a <- sfarrow::st_read_parquet(ans$path[ans$tile_id == "A"])
  expect_equal(nrow(a), 3L)
  expect_true(all(a$lapa == "A"))
})

test_that("tile_vector_grid falls back to deterministic row chunks", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs"))

  grid <- sf::st_as_sf(
    data.frame(x = 1:5, y = 1:5),
    coords = c("x", "y"),
    crs = 3059
  )

  src_dir <- tempfile("tile-src-")
  out_dir <- tempfile("tile-out-")
  dir.create(src_dir)
  grid_path <- file.path(src_dir, "pts100_sauzeme.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(grid, grid_path))

  ans <- suppress_sfarrow_geo_warning(tile_vector_grid(
    grid_path = grid_path,
    out_dir = out_dir,
    tile_field = NULL,
    chunk_size = 2L,
    quiet = TRUE
  ))

  expect_equal(ans$tile_id, c("tile_00001", "tile_00002", "tile_00003"))
  expect_equal(ans$n_rows, c(2L, 2L, 1L))
  expect_true(all(ans$wrote))
})

test_that("tile_vector_grid skips existing tiles", {
  skip_if_spatial_packages_missing(c("sf", "sfarrow", "fs"))

  grid <- sf::st_as_sf(
    data.frame(x = 1:2, y = 1:2, lapa = c("A", "A")),
    coords = c("x", "y"), crs = 3059
  )
  src_dir <- tempfile("tile-src-")
  out_dir <- tempfile("tile-out-")
  dir.create(src_dir)
  grid_path <- file.path(src_dir, "tikls100_sauzeme.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(grid, grid_path))

  first <- suppress_sfarrow_geo_warning(tile_vector_grid(grid_path, out_dir, tile_field = "lapa", quiet = TRUE))
  second <- suppress_sfarrow_geo_warning(tile_vector_grid(grid_path, out_dir, tile_field = "lapa", overwrite = FALSE, quiet = TRUE))

  expect_true(first$wrote)
  expect_false(second$wrote)
})
