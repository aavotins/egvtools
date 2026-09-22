test_that("radius_function validates layer-prefix and extract-function lengths", {
  skip_if_spatial_packages_missing(c(
    "terra", "sf", "raster", "dplyr", "tibble", "purrr", "fs", "sfarrow",
    "fasterize", "exactextractr", "future", "furrr", "whitebox", "glue", "tidyr"
  ))

  expect_error(
    radius_function(
      kvadrati_path = tempdir(),
      radii_path = tempdir(),
      tikls100_path = "x",
      template_path = "x",
      input_layers = c("a", "b"),
      layer_prefixes = "a",
      output_dir = tempfile("radius-") ,
      fill_missing = FALSE,
      quiet = TRUE
    )
  )

  expect_error(
    radius_function(
      kvadrati_path = tempdir(),
      radii_path = tempdir(),
      tikls100_path = "x",
      template_path = "x",
      input_layers = "a",
      layer_prefixes = "a",
      extract_fun = c("mean", "sum"),
      output_dir = tempfile("radius-") ,
      fill_missing = FALSE,
      quiet = TRUE
    ),
    "extract_fun"
  )
})

test_that("radius_function completes a one-tile one-radius mocked extraction", {
  skip_if_spatial_packages_missing(c(
    "terra", "sf", "raster", "dplyr", "tibble", "purrr", "fs", "sfarrow",
    "fasterize", "exactextractr", "future", "furrr", "whitebox", "glue", "tidyr"
  ))
  skip_if_not_installed("testthat", minimum_version = "3.1.0")

  local_mocked_bindings(
    exact_extract = function(x, y, fun, ...) seq_len(nrow(y)),
    .package = "exactextractr"
  )

  old_env <- Sys.getenv(c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"), unset = NA_character_)
  on.exit(restore_envvars(old_env), add = TRUE)
  Sys.setenv(SLURM_CPUS_PER_TASK = "1", SLURM_CPUS_ON_NODE = "1")

  kv_dir <- tempfile("kv-")
  rad_dir <- tempfile("rad-")
  out_dir <- tempfile("radius-out-")
  dir.create(kv_dir)
  dir.create(rad_dir)

  cells <- make_two_zone_sf()
  cells$id <- c("1", "2")
  cells$rinda300 <- c("r1", "r2")
  cells$ID1km <- c("k1", "k2")
  kv_file <- file.path(kv_dir, "tikls100_2434.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(cells, kv_file))

  buffers <- sf::st_sf(
    id = cells$id,
    geometry = sf::st_buffer(sf::st_centroid(sf::st_geometry(cells)), dist = 0.75)
  )
  rad_file <- file.path(rad_dir, "pts100_r500_2434.parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(buffers, rad_file))

  tikls100_file <- tempfile(fileext = ".parquet")
  suppress_sfarrow_geo_warning(sfarrow::st_write_parquet(cells, tikls100_file))

  template <- make_test_raster(values = rep(1, 16))
  template_file <- tempfile(fileext = ".tif")
  write_test_raster(template, template_file)

  cov <- make_test_raster(values = 1:16)
  cov_file <- tempfile(fileext = ".tif")
  write_test_raster(cov, cov_file)

  ans <- radius_function(
    kvadrati_path = kv_dir,
    radii_path = rad_dir,
    tikls100_path = tikls100_file,
    template_path = template_file,
    input_layers = cov_file,
    layer_prefixes = "cov",
    output_dir = out_dir,
    unlink_tiles = TRUE,
    n_workers = 1,
    radii = "r500",
    fill_missing = FALSE,
    radius_mode = "sparse",
    extract_fun = "mean",
    quiet = TRUE
  )

  expect_s3_class(ans, "data.frame")
  expect_equal(nrow(ans), 1L)
  expect_equal(ans$layer, "cov")
  expect_equal(ans$radius, "r500")
  expect_equal(ans$n_tiles_merged, 1L)
  expect_true(file.exists(ans$output_path))

  out <- terra::rast(ans$output_path)
  expect_same_grid(out, template)
  expect_equal(names(out), "cov_r500")
})
