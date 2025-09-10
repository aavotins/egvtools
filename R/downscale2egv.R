#' Downscale & Align a Raster to a Template, with optional gap fill & IDW smoothing
#'
#' @description
#' Aligns a source raster to a template grid (CRS, resolution, extent), masks to the
#' template footprint, and optionally:
#' (1) fills NoData gaps using WhiteboxTools' IDW-based `fill_missing_data`, and
#' (2) applies **IDW smoothing** to reduce blockiness from low-resolution inputs.
#'
#' **Mass preservation:** IDW smoothing (and most simple smoothers) are **not mass-preserving**.
#' If conservation of totals/means matters (globally or by zones), post-smoothing rescaling
#' is recommended (see *Details*).
#'
#' If `interpolation_method = "auto"`, integer/factor inputs use `"near"`;
#' continuous inputs use `"bilinear"`.
#'
#' @details
#' **Workflow**
#' 1) **Inputs**: accepts paths or in-memory objects (`SpatRaster` / `sf`).
#' 2) **Fast crop**: transforms the grid bbox to the template CRS, expands by `buffer_m`,
#'    crops the raw raster early to minimise work.
#' 3) **Method auto-choice**: `"near"` for integer/factor rasters, otherwise `"bilinear"`.
#' 4) **Align**: `terra::project(raw_crop, template, method)` then `terra::mask(..., template)`.
#' 5) **Gaps**: optional internal NA counts on the template footprint only.
#' 6) **Gap fill (optional)**: Whitebox `wbt_fill_missing_data`, with window auto-sized from
#'    the widest gap (or user `filter_size_cells`).
#' 7) **Smoothing (optional)**: IDW via `terra::interpIDW()` on points derived from an
#'    aggregated version of the filled raster; masked back to template.
#' 8) **Plotting**: `plot_gaps`, `plot_result` (side-by-side if both).
#' 9) **Write**: atomic write (tmp + rename), LZW, tiling, `BIGTIFF=IF_SAFER`, threaded;
#'    output CRS **forced to template WKT/EPSG** so `crs(out) == crs(template)` is TRUE.
#' 10) **Sink safety & memory**: snapshots/restore sinks (prevents stuck console if interrupted);
#'     `terraOptions(memfrac/tempdir/todisk)` set then restored; optional `gc()`.
#'
#' **IDW smoothing**
#' Aggregate to points to IDW back to template. Controls:
#' - `smooth_radius_km` search radius (km; converted to template units),
#' - `smooth_agg_factor` aggregation factor before point creation,
#' - `smooth_power`, `smooth_epsilon`, `smooth_nmax`.
#' Set `smooth_force=TRUE` to allow smoothing of integer/factor rasters (not recommended).
#'
#' **Mass preservation tips**
#' The smoothed result will generally *not* preserve totals/means.
#' - **Global**: scale smoothed raster to match global sum.
#' - **Zonal**: compute per-zone factors on the original coarse grid and multiply.
#'
#' @param template_path Character path to a template GeoTIFF **or** a `terra::SpatRaster`.
#' @param grid_path Character path to a GeoParquet grid **or** an `sf` object (used only for bbox).
#' @param rawfile_path Character path to the input raster **or** a `SpatRaster`.
#' @param out_path Character. Directory for the output GeoTIFF (created if missing).
#' @param file_name Character. Output file name (e.g., `"result.tif"`).
#' @param layer_name Character. Band name to set on the output.
#' @param interpolation_method `"auto"` (default), or `"bilinear"`, `"near"`, `"cubicspline"`, `"cubic"`.
#' @param buffer_m Numeric. Buffer distance (meters) applied to the grid bbox before cropping. Default `10000`.
#' @param check_na Logical. If `TRUE`, report NA gap count (inside template). Default `FALSE`.
#' @param fill_gaps Logical. If `TRUE`, fill NA gaps via Whitebox. Default `TRUE`.
#' @param idw_weight Numeric. IDW power used by Whitebox gap filling. Default `2`.
#' @param filter_size_cells Integer or `"auto"`. Window size (in **cells**) for gap filling. Default `"auto"`.
#' @param plot_gaps,plot_result Logical flags for diagnostics/plots (side-by-side if both).
#' @param smooth Logical. Enable **IDW** smoothing. Default `FALSE`.
#' @param smooth_radius_km Numeric. IDW smoothing radius (km). Default `10`.
#' @param smooth_agg_factor Integer. Aggregate factor before point creation. Default `10`.
#' @param smooth_power Numeric. IDW power. Default `0.5`.
#' @param smooth_epsilon Numeric. Small smoothing term (if supported by your `terra`). Default `0`.
#' @param smooth_nmax Integer. Max neighbors (if supported by your `terra`). Default `50`.
#' @param smooth_force Logical. Allow smoothing for integer/factor rasters (not recommended). Default `FALSE`.
#' @param gdal_opts Character vector. GDAL creation options; merged with tuned defaults:
#'   `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Character or `NULL`. Optional terra datatype (e.g., `"FLT4S"`).
#' @param NAflag Numeric/integer NoData value to write; `NULL` to omit (default).
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra operations. Default `tempdir()`.
#' @param terra_todisk Logical or `NA`. If `TRUE`, prefer on-disk processing for heavy ops. Default `TRUE`.
#' @param force_gc Logical; call `gc()` at checkpoints. Default `FALSE`.
#' @param return_visible Logical. If `TRUE`, return the summary visibly; otherwise invisibly. Default `FALSE`.
#' @param quiet Suppress progress prints (`cat()`)? Default `FALSE`.
#'
#' @return A **data.frame** with columns:
#' `output`, `gap_count`, `max_gap_distance`, `filter_size_cells`, `smoothed`,
#' `smoothing_radius_used`, `elapsed_sec`.
#'
#' @examples
#' \dontrun{
#' out_dir <- tempdir()
#' df <- downscale2egv(
#'   template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   grid_path     = "./Templates/TemplateGrids/grid_10km.parquet",
#'   rawfile_path  = "./MyCoarse/coarse_indicator.tif",
#'   out_path      = out_dir,
#'   file_name     = "indicator_egv.tif",
#'   layer_name    = "indicator",
#'   fill_gaps     = TRUE,
#'   smooth        = TRUE,
#'   smooth_radius_km = 10,
#'   plot_result   = TRUE
#' )
#' print(df)
#' }
#'
#' @seealso terra::project(), terra::mask(), terra::interpIDW(),
#'          whitebox::wbt_fill_missing_data(), sfarrow::st_read_parquet()
#'
#' @import terra
#' @import sf
#' @import sfarrow
#' @import whitebox
#' @import fs
#' @importFrom utils capture.output
#' @importFrom graphics par
#' @export
downscale2egv <- function(
    template_path,
    grid_path,
    rawfile_path,
    out_path,
    file_name,
    layer_name,
    interpolation_method = "auto",
    buffer_m = 10000,
    check_na = FALSE,
    fill_gaps = TRUE,
    idw_weight = 2,
    filter_size_cells = "auto",
    plot_gaps = FALSE,
    plot_result = FALSE,
    smooth = FALSE,
    smooth_radius_km = 10,
    smooth_agg_factor = 10,
    smooth_power = 0.5,
    smooth_epsilon = 0,
    smooth_nmax = 50,
    smooth_force = FALSE,
    gdal_opts = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
    write_datatype = NULL,
    NAflag = NULL,
    terra_memfrac = 0.7,
    terra_tempdir = tempdir(),
    terra_todisk = TRUE,
    force_gc = FALSE,
    return_visible = FALSE,
    quiet = FALSE
) {
  t0 <- Sys.time()

  # ---- sink safety ----
  orig_out <- sink.number(); orig_msg <- sink.number(type = "message")
  on.exit({ while (sink.number(type="message") > orig_msg) sink(type="message")
    while (sink.number() > orig_out) sink() }, add = TRUE)

  say <- function(...) if (!quiet) cat(..., "\n")
  .maybe_gc <- function() if (isTRUE(force_gc)) gc()


  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("terra", "project/mask/resample/write")
  .need_pkg("sf",    "bbox transform for fast crop")
  .need_pkg("sfarrow", "reading GeoParquet 'grid_path'")
  .need_pkg("whitebox", "IDW gap filling")
  .need_pkg("fs", "filesystem ops")

  # ---- I/O checks
  if (!fs::dir_exists(out_path)) fs::dir_create(out_path, recurse = TRUE)
  if (missing(file_name)  || !nzchar(file_name))  stop("'file_name' must be provided.")
  if (missing(layer_name) || !nzchar(layer_name)) stop("'layer_name' must be provided.")
  out_full <- normalizePath(file.path(out_path, file_name), winslash = "/", mustWork = FALSE)

  # ---- terra options (set + restore)
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac, tempdir = old_opt$tempdir,
                              todisk = old_opt$todisk, progress = old_opt$progress), add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir, progress = FALSE)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
  })

  # ---- helpers
  as_spat_raster <- function(x, role) {
    if (inherits(x, "SpatRaster")) return(x)
    if (is.character(x)) { if (!file.exists(x)) stop(sprintf("'%s' path not found for %s.", x, role)); return(terra::rast(x)) }
    stop(sprintf("%s must be a file path or a terra::SpatRaster.", role))
  }
  as_sf <- function(x, role) {
    if (inherits(x, "sf")) return(x)
    if (is.character(x)) { if (!file.exists(x)) stop(sprintf("'%s' path not found for %s.", x, role)); return(sfarrow::st_read_parquet(x, as_tibble = FALSE, quiet = TRUE)) }
    stop(sprintf("%s must be a file path or an sf object.", role))
  }

  # ---- load inputs
  template_r <- as_spat_raster(template_path, "template raster")
  template_wkt <- terra::crs(template_r)  # ensure WKT/EPSG form
  grid_sf    <- as_sf(grid_path, "grid")
  raw_r      <- as_spat_raster(rawfile_path, "raw raster")

  # ---- Fast bbox transform & early crop
  tpl_crs <- sf::st_crs(template_wkt)
  bb0     <- sf::st_bbox(grid_sf)
  bbox0   <- sf::st_as_sfc(bb0, crs = sf::st_crs(grid_sf))
  bbox_tpl <- try(sf::st_transform(bbox0, tpl_crs), silent = TRUE)
  if (inherits(bbox_tpl, "try-error") || is.null(sf::st_crs(bbox_tpl))) bbox_tpl <- bbox0

  bb <- sf::st_bbox(bbox_tpl)
  if (is.finite(buffer_m) && buffer_m > 0) {
    bb["xmin"] <- bb["xmin"] - buffer_m; bb["xmax"] <- bb["xmax"] + buffer_m
    bb["ymin"] <- bb["ymin"] - buffer_m; bb["ymax"] <- bb["ymax"] + buffer_m
  }

  bbox_tpl_v <- terra::vect(sf::st_as_sfc(bb))              # in template CRS
  bbox_raw_v <- try(terra::project(bbox_tpl_v, terra::crs(raw_r)), silent = TRUE)

  raw_crop <- if (inherits(bbox_raw_v, "try-error")) raw_r else terra::crop(raw_r, terra::ext(bbox_raw_v))
  if (terra::ncell(raw_crop) == 0) stop("Pre-crop produced an empty raster (0 cells). Check CRS/extent overlap or increase 'buffer_m'.")

  # ---- Interpolation method auto-switch
  is_categorical <- any(terra::is.factor(raw_crop)) || all(grepl("^INT", terra::datatype(raw_crop)))
  method <- if (identical(interpolation_method, "auto")) if (is_categorical) "near" else "bilinear" else interpolation_method

  # ---- Align: project -> mask
  say(sprintf("Projecting (%s) and masking to template ...", method))
  projected <- terra::project(raw_crop, template_r, method = method,
                              filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("dwn_proj_", basename(tempfile("")), ".tif")) else "",
                              overwrite = TRUE)
  if (terra::ncell(projected) == 0) stop("Projection produced an empty raster (0 cells).")

  masked <- terra::mask(projected, template_r,
                        filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("dwn_mask_", basename(tempfile("")), ".tif")) else "",
                        overwrite = TRUE)
  if (terra::ncell(masked) == 0) stop("Masking produced an empty raster (0 cells).")
  rm(projected); .maybe_gc()

  # ---- Gaps on template footprint
  gaps_r <- is.na(masked) & !is.na(template_r)
  gap_count <- NA_real_; max_gap_distance <- NA_real_

  if (isTRUE(check_na) || isTRUE(plot_gaps) || (isTRUE(fill_gaps) && identical(filter_size_cells, "auto"))) {
    gap_count <- terra::global(gaps_r, fun = "sum", na.rm = TRUE)[[1]]
    if (isTRUE(fill_gaps) && identical(filter_size_cells, "auto") && !is.na(gap_count) && gap_count > 0) {
      data_mask <- terra::ifel(gaps_r, NA, 1)
      dist_r <- terra::distance(data_mask)
      max_gap_distance <- terra::global(dist_r, fun = "max", na.rm = TRUE)[[1]]
      rm(dist_r); .maybe_gc()
    }
  }
  if (isTRUE(check_na)) say(sprintf("Gap cells (within template): %s", ifelse(is.na(gap_count), "NA", format(gap_count, big.mark=","))))

  # ---- Optional gap fill (Whitebox)
  used_filter <- NA_integer_
  filled_r <- masked

  if (isTRUE(fill_gaps) && !is.na(gap_count) && gap_count > 0) {
    if (!requireNamespace("whitebox", quietly = TRUE)) {
      warning("Package 'whitebox' not installed; skipping gap filling.")
    } else {
      if (identical(filter_size_cells, "auto")) {
        cell_m <- mean(terra::res(template_r))
        fs <- ceiling((max_gap_distance / cell_m) * 2)
        if (fs %% 2 == 0) fs <- fs + 1
        used_filter <- max(3L, as.integer(fs))
      } else {
        used_filter <- max(3L, as.integer(filter_size_cells))
        if (used_filter %% 2 == 0) used_filter <- used_filter + 1
      }

      in_wbt  <- file.path(terra_tempdir, paste0("dwn_wbt_in_",  basename(tempfile("")), ".tif"))
      out_wbt <- file.path(terra_tempdir, paste0("dwn_wbt_out_", basename(tempfile("")), ".tif"))

      ok <- TRUE
      tryCatch({
        terra::writeRaster(masked, in_wbt, overwrite = TRUE,
                           gdal = c("COMPRESS=NONE","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256"))
      }, error = function(e) { warning("Could not write temporary raster for Whitebox; skipping fill. ", conditionMessage(e)); ok <<- FALSE })

      if (ok) {
        wbt_path <- try(whitebox::wbt_exe_path(), silent = TRUE)
        wbt_ok <- !inherits(wbt_path, "try-error") && nzchar(wbt_path)
        if (wbt_ok) {
          say(sprintf("Whitebox fill: filter=%d, weight=%.2f", as.integer(used_filter), idw_weight))
          res_wbt <- try(whitebox::wbt_fill_missing_data(input=in_wbt, output=out_wbt,
                                                         filter=as.integer(used_filter), weight=idw_weight,
                                                         no_edges=FALSE, verbose_mode=FALSE), silent = TRUE)
          if (!inherits(res_wbt, "try-error") && file.exists(out_wbt)) {
            filled_r <- terra::rast(out_wbt)
            names(filled_r) <- names(masked)
            filled_r <- terra::mask(filled_r, template_r)
          } else {
            warning("Whitebox gap filling failed; keeping masked raster.")
          }
        } else {
          warning("WhiteboxTools executable not configured; skipping gap filling.")
        }
      }
      unlink(c(in_wbt, out_wbt), force = TRUE)
    }
  } else if (isTRUE(fill_gaps)) {
    say("No gaps detected; skipping gap filling.")
  }
  rm(masked, gaps_r); .maybe_gc()

  # ---- Optional smoothing (IDW)
  smoothed <- FALSE
  smoothing_radius_used <- NA_real_

  if (isTRUE(smooth)) {
    is_cat_now <- any(terra::is.factor(filled_r)) || all(grepl("^INT", terra::datatype(filled_r)))
    if (is_cat_now && !isTRUE(smooth_force)) {
      warning("Skipping smoothing: raster is integer/factor. Set smooth_force=TRUE to override.")
    } else {
      if (terra::is.lonlat(template_r)) warning("Template is lon/lat; IDW radius will be in degrees after km to deg conversion.")
      km_to_units <- function(km) if (terra::is.lonlat(template_r)) km/111.32 else km*1000
      agg_fac <- max(1L, as.integer(smooth_agg_factor))
      aggr <- try(terra::aggregate(filled_r, agg_fac, fun = mean, na.rm = TRUE), silent = TRUE)
      if (!inherits(aggr, "try-error") && terra::ncell(aggr) > 0) {
        pts <- try(terra::as.points(aggr, values = TRUE, na.rm = TRUE), silent = TRUE)
        if (!inherits(pts, "try-error") && inherits(pts, "SpatVector") && terra::nrow(pts) > 0) {
          fld <- names(pts)[1]
          radius_units <- km_to_units(smooth_radius_km)
          cell_u <- mean(terra::res(template_r))
          rad0 <- max(radius_units, 1.5 * cell_u)

          safe_formals <- function(f) { fml <- try(formals(f), silent = TRUE); if (inherits(fml, "try-error")) character() else names(fml) }
          idw_formals <- safe_formals(terra::interpIDW)
          mk_args <- function(radius_val) {
            a <- list(x = template_r, y = pts, field = fld, radius = radius_val, power = smooth_power)
            if ("smooth" %in% idw_formals) a$smooth <- smooth_epsilon
            if ("nmax"   %in% idw_formals) a$nmax   <- as.integer(smooth_nmax)
            a
          }

          mults <- c(1, 1.5, 3, 5, 10)
          smoothed_r <- NULL
          for (m in mults) {
            rad_try <- rad0 * m
            args_try <- mk_args(rad_try)
            smoothed_r <- try(do.call(terra::interpIDW, args_try), silent = TRUE)
            if (!inherits(smoothed_r, "try-error")) { smoothing_radius_used <- rad_try; break }
          }
          if (!inherits(smoothed_r, "try-error") && !is.null(smoothed_r)) {
            smoothed_r <- terra::mask(smoothed_r, template_r)
            names(smoothed_r) <- layer_name
            filled_r <- smoothed_r
            smoothed <- TRUE
          } else {
            warning("Smoothing: interpIDW failed after retries; skipping smoothing.")
          }
        } else {
          warning("Smoothing: could not create points; skipping smoothing.")
        }
      } else {
        warning("Smoothing: aggregation failed or produced 0 cells; skipping smoothing.")
      }
    }
  }

  # ---- Name, plot, write
  names(filled_r) <- layer_name

  if ((isTRUE(plot_gaps) || isTRUE(plot_result)) && interactive()) {
    gaps_plot <- NULL
    if (isTRUE(plot_gaps)) gaps_plot <- is.na(filled_r) & !is.na(template_r)
    oldpar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(oldpar), add = TRUE)
    if (isTRUE(plot_gaps) && isTRUE(plot_result)) graphics::par(mfrow = c(1, 2))
    if (isTRUE(plot_gaps))   terra::plot(gaps_plot, main = sprintf("Gaps (NA) before fill: %s", ifelse(is.na(gap_count),"NA", format(gap_count, big.mark=","))))
    if (isTRUE(plot_result)) terra::plot(filled_r, main = layer_name)
  }

  if (terra::ncell(filled_r) == 0) stop("Nothing to write: result has 0 cells.")

  # GDAL opts
  def_opts <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")
  gdal_opts <- unique(c(gdal_opts, def_opts))
  out_is_float <- any(grepl("^FLT", terra::datatype(filled_r)))
  if (out_is_float && !any(grepl("^PREDICTOR=", gdal_opts))) gdal_opts <- c(gdal_opts, "PREDICTOR=2")

  # enforce WKT/EPSG CRS for exact string equality with template
  terra::crs(filled_r) <- template_wkt

  # atomic write
  ext  <- tools::file_ext(out_full); ext <- if (nzchar(ext)) ext else "tif"
  stem <- sub(sprintf("\\.%s$", ext), "", out_full)
  tmp_out <- sprintf("%s._tmp.%s", stem, ext)
  if (file.exists(tmp_out)) try(unlink(tmp_out), silent = TRUE)

  write_args <- list(x = filled_r, filename = tmp_out, overwrite = TRUE, gdal = gdal_opts)
  if (!is.null(write_datatype)) write_args$datatype <- write_datatype
  if (!is.null(NAflag))        write_args$NAflag   <- NAflag
  do.call(terra::writeRaster, write_args)

  if (file.exists(out_full)) try(unlink(out_full), silent = TRUE)
  file.rename(tmp_out, out_full)
  say(sprintf("Wrote: %s", out_full))

  # ---- timing & return as data.frame ----
  elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  df <- data.frame(
    output                 = out_full,
    gap_count              = if (is.na(gap_count)) NA_real_ else as.numeric(gap_count),
    max_gap_distance       = if (is.na(max_gap_distance)) NA_real_ else as.numeric(max_gap_distance),
    filter_size_cells      = if (isTRUE(fill_gaps)) as.integer(used_filter) else NA_integer_,
    smoothed               = isTRUE(smoothed),
    smoothing_radius_used  = if (is.na(smoothing_radius_used)) NA_real_ else as.numeric(smoothing_radius_used),
    elapsed_sec            = elapsed_sec,
    stringsAsFactors = FALSE
  )

  if (isTRUE(return_visible)) df else invisible(df)
}
