#' Create constant-background rasters from a directory of GeoTIFFs
#'
#' @description
#' For every raster in `in_dir`, creates a new raster where all **non-NA** cells
#' are set to `background_value` and NA cells are preserved. Outputs go to `out_dir`
#' (created if needed). Files are named with a prefix (by default `nulls_` when
#' `background_value == 0`, otherwise `bg[background_value]_`) followed by the
#' original name, e.g. `bg10_my.tif`. Uses LZW compression and atomic writes.
#' Large rasters can stream to disk with `terra_todisk=TRUE`.
#'
#'
#' @details
#' - Reads each input with `terra::rast()`.
#' - Applies a fast element-wise replacement with `terra::ifel(!is.na(x), background_value, NA)`.
#' - Writes with LZW compression to a temporary `._tmp.tif`, then atomically renames
#'   to the final file.
#' - If you do not set `write_datatype`/`NAflag`, the function picks sensible defaults:
#'     - If `background_value` is an integer and within `[-32768, 32767]`, then `INT2S` + `NAflag=-32768`.
#'     - Otherwise `FLT4S` + float `NAflag` (`-3.4028235e38`).
#' - Only `.tif` / `.tiff` files are processed (optionally refine with `pattern`).
#'
#' @param in_dir Directory containing input rasters.
#' @param out_dir Directory for outputs. If `NULL`, defaults to `file.path(in_dir,"backgrounds")`.
#' @param background_value Numeric value to assign to all non-NA pixels. Default `0`.
#' @param pattern Optional regex to filter filenames (applied after extension filter).
#' @param recursive Logical; search subdirectories in `in_dir`. Default `FALSE`.
#' @param out_prefix Optional filename prefix. If `NULL`, it defaults to:
#'   - `"nulls_"` when `background_value == 0`,
#'   - otherwise `"bg[background_value]_"` (e.g., `"bg10_"`, `"bg0.5_"`).
#' @param NAflag Optional NA flag for writing (passed to GDAL). Default `NULL` (auto).
#' @param gdal_opts Character vector of GDAL creation options (merged with tuned defaults).
#'   Default `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Optional terra datatype for writing (e.g., `"INT2S"`, `"FLT4S"`).
#'   Default `NULL` (auto as described above).
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra operations. Default `tempdir()`.
#' @param terra_todisk Logical or `NA`. If `TRUE`, prefer on-disk processing for the
#'   `ifel` step. Default `FALSE`.
#' @param force_gc Logical; call `gc()` at checkpoints. Default `FALSE`.
#' @param overwrite Overwrite existing outputs? Default `FALSE`.
#' @param quiet Suppress progress prints (`cat()`)? Default `FALSE`.
#'
#' @return Invisibly, a data.frame with columns: `in_file`, `out_file`, `n_cells`, `elapsed_sec`.
#'
#' @examples
#' \dontrun{
#' # Classic "nulls_" backgrounds (fill 0)
#' create_backgrounds(
#'   in_dir = "./Templates/TemplateRasters",
#'   out_dir = "./Outputs/Nulls",
#'   background_value = 0,
#'   overwrite = TRUE,
#'   terra_todisk = TRUE
#' )
#'
#' # Fill non-NA with 10, prefix becomes "bg10_"
#' create_backgrounds(
#'   in_dir = "./Templates/TemplateRasters",
#'   out_dir = "./Outputs/BG10",
#'   background_value = 10,
#'   overwrite = TRUE
#' )
#' }
#'
#' @importFrom fs dir_exists dir_create dir_ls path_ext path_ext_remove path_file file_exists file_delete file_move
#' @importFrom terra rast ncell ifel writeRaster
#' @importFrom utils capture.output
#' @export
create_backgrounds <- function(
    in_dir,
    out_dir = NULL,
    background_value = 0,
    pattern = NULL,
    recursive = FALSE,
    out_prefix = NULL,
    NAflag = NULL,
    gdal_opts = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256"),
    write_datatype = NULL,
    terra_memfrac = 0.7,
    terra_tempdir = tempdir(),
    terra_todisk = FALSE,
    force_gc = FALSE,
    overwrite = FALSE,
    quiet = FALSE
) {
  t_all <- Sys.time()

  # ---- sink safety: snapshot & restore on exit
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)

  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("terra", "read / write / value ops")
  .need_pkg("fs",    "dir creation, file listing")


  # ---- validate dirs
  if (missing(in_dir) || !fs::dir_exists(in_dir)) stop("'in_dir' must be an existing directory.")
  if (is.null(out_dir)) out_dir <- file.path(in_dir, "backgrounds")
  if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)

  # ---- terra options (silenced)
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac, tempdir = old_opt$tempdir, todisk = old_opt$todisk), add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
  })

  say <- function(...) if (!quiet) cat(..., "\n")
  tuned_defaults <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")
  gdal_opts <- unique(c(gdal_opts, tuned_defaults))

  # auto prefix if not specified
  if (is.null(out_prefix)) {
    if (isTRUE(all.equal(background_value, 0))) out_prefix <- "nulls_" else out_prefix <- paste0("bg", background_value, "_")
  }

  # auto datatype / NAflag (unless provided)
  if (is.null(write_datatype)) {
    # choose INT2S if background_value is an integer within [-32768, 32767]; else FLT4S
    is_intlike <- isTRUE(all.equal(background_value, as.integer(background_value)))
    if (is_intlike && background_value >= -32768 && background_value <= 32767) {
      write_datatype <- "INT2S"
    } else {
      write_datatype <- "FLT4S"
    }
  }
  if (is.null(NAflag)) {
    NAflag <- if (write_datatype == "INT2S") -32768 else -3.4028235e38
  }

  # ---- collect files
  files <- fs::dir_ls(in_dir, recurse = recursive, type = "file")
  files <- files[tolower(fs::path_ext(files)) %in% c("tif","tiff")]
  if (!is.null(pattern)) files <- files[grepl(pattern, fs::path_file(files))]
  if (!length(files)) {
    say("No .tif/.tiff files found in:", in_dir)
    return(invisible(data.frame(in_file=character(), out_file=character(), n_cells=integer(), elapsed_sec=numeric())))
  }

  say(sprintf("Found %d raster(s). Writing to: %s", length(files), out_dir))
  results <- vector("list", length(files))

  for (i in seq_along(files)) {
    t0 <- Sys.time()
    in_file <- files[i]
    base_noext <- fs::path_file(fs::path_ext_remove(in_file))
    out_file  <- file.path(out_dir, paste0(out_prefix, base_noext, ".tif"))
    if (fs::file_exists(out_file) && !overwrite) {
      say(sprintf("[%d/%d] Exists (skip): %s", i, length(files), out_file))
      results[[i]] <- data.frame(in_file=in_file, out_file=out_file, n_cells=NA_integer_, elapsed_sec=NA_real_)
      next
    }

    say(sprintf("[%d/%d] Processing: %s", i, length(files), fs::path_file(in_file)))
    r <- terra::rast(in_file)

    # Assign background to all non-NA; keep NA as NA (stream to disk if requested)
    z <- terra::ifel(!is.na(r), background_value, NA,
                     filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, sprintf("cbg_%s.tif", basename(tempfile("")))) else "",
                     overwrite = TRUE)

    # CRS guard: force identical WKT as input
    terra::crs(z) <- terra::crs(r, proj = FALSE)

    # Atomic write with guarded args (preserve .tif)
    ext  <- fs::path_ext(out_file); ext <- if (nzchar(ext)) ext else "tif"
    stem <- fs::path_ext_remove(out_file)
    tmp_out <- sprintf("%s._tmp.%s", stem, ext)
    if (fs::file_exists(tmp_out)) fs::file_delete(tmp_out)

    write_args <- list(x = z, filename = tmp_out, overwrite = TRUE, gdal = gdal_opts)
    if (!is.null(write_datatype)) write_args$datatype <- write_datatype
    if (!is.null(NAflag))        write_args$NAflag   <- NAflag
    do.call(terra::writeRaster, write_args)

    if (fs::file_exists(out_file)) fs::file_delete(out_file)
    fs::file_move(tmp_out, out_file)

    if (force_gc) gc()

    results[[i]] <- data.frame(
      in_file = in_file,
      out_file = out_file,
      n_cells = terra::ncell(r),
      elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
    say(sprintf("  -> Wrote: %s", out_file))
  }

  res_df <- do.call(rbind, results)
  say(sprintf("Done. Total elapsed: %.1f sec", as.numeric(difftime(Sys.time(), t_all, units = "secs"))))
  invisible(res_df)
}
