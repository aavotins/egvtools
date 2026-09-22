skip_if_spatial_packages_missing <- function(pkgs = character()) {
  for (pkg in pkgs) testthat::skip_if_not_installed(pkg)
}

make_test_raster <- function(
    nrows = 4L,
    ncols = 4L,
    xmin = 0,
    xmax = ncols,
    ymin = 0,
    ymax = nrows,
    values = seq_len(nrows * ncols),
    crs = "EPSG:3059"
) {
  r <- terra::rast(
    nrows = nrows,
    ncols = ncols,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    crs = crs
  )
  terra::values(r) <- values
  r
}

write_test_raster <- function(r, path) {
  terra::writeRaster(r, path, overwrite = TRUE)
  path
}

make_rect_sf <- function(xmin, xmax, ymin, ymax, crs = 3059, ...) {
  xy <- matrix(
    c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ),
    ncol = 2,
    byrow = TRUE
  )
  geom <- sf::st_sfc(sf::st_polygon(list(xy)), crs = crs)
  sf::st_sf(..., geometry = geom)
}

make_two_zone_sf <- function(crs = 3059) {
  p1 <- matrix(c(0, 0, 2, 0, 2, 4, 0, 4, 0, 0), ncol = 2, byrow = TRUE)
  p2 <- matrix(c(2, 0, 4, 0, 4, 4, 2, 4, 2, 0), ncol = 2, byrow = TRUE)
  sf::st_sf(
    id = c("z1", "z2"),
    tile = c("A", "A"),
    geometry = sf::st_sfc(sf::st_polygon(list(p1)), sf::st_polygon(list(p2)), crs = crs)
  )
}

make_test_points_sf <- function(crs = 3059) {
  sf::st_as_sf(
    data.frame(
      x = c(10, 20, 70, 80),
      y = c(10, 20, 70, 80),
      tks50km = c("A", "A", "B", "B"),
      radius_m = c(5, 10, 15, 20)
    ),
    coords = c("x", "y"),
    crs = crs
  )
}

expect_same_grid <- function(x, y) {
  # Compare the raster grid independently of the CRS representation.
  # GDAL/PROJ may serialize an equivalent CRS differently across platforms
  # (notably Windows versus macOS/Linux), so a single compareGeom(crs = TRUE)
  # assertion can be unnecessarily platform-sensitive.
  testthat::expect_true(
    terra::compareGeom(
      x, y,
      stopOnError = FALSE,
      crs = FALSE,
      ext = TRUE,
      rowcol = TRUE,
      res = TRUE
    ),
    info = "Raster extent, dimensions, or resolution differ"
  )

  # Test CRS equivalence semantically rather than comparing serialized WKT text.
  testthat::expect_true(
    terra::same.crs(x, y),
    info = "Raster coordinate reference systems are not equivalent"
  )
}

restore_envvars <- function(old) {
  for (nm in names(old)) {
    if (is.na(old[[nm]])) {
      Sys.unsetenv(nm)
    } else {
      do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
    }
  }
}


suppress_sfarrow_geo_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl(
        "initial implementation of Parquet/Feather file support",
        conditionMessage(w),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
