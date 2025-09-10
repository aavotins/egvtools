#' Convert an input raster to an EGV-aligned raster
#'
#' @description
#' Align a fine-resolution input raster to a (coarser) EGV template, optionally
#' cover missing values and/or fill gaps (IDW via Whitebox), and write the result to disk.
#' Designed for large runs: fast gap counting (inside template footprint only), optional filling,
#' tuned GDAL write options, and controlled `terra` memory/temp behavior.
#'
#' @details
#' **Workflow (high level)**
#' 1. Load `input` and `egv_template`. Optionally warn on CRS/resolution mismatch (`check_alignment`).
#' 2. If `missing_job` contains **CoverInput**, build `input_bg` (numeric + `input_template`, or raster/path)
#'    and `terra::cover()` the input.
#' 3. **Resample/aggregate** the input to the template with `summary_function`
#'    (`"mean"` is linked to `"average"`). If `is_categorical=TRUE` the method is forced to `"near"`.
#' 4. Count **initial gaps** = `is.na(aligned) & !is.na(template)`.
#' 5. If `missing_job` contains **CoverOutput**, apply `terra::cover()` using `output_bg`.
#' 6. If `missing_job` contains **FillOutput** **and** `!is_categorical`,
#'    estimate a **maximum gap width** (via `terra::distance()` on a fillable mask) and set the
#'    Whitebox **filter width** as `ceil(max_gap/pixel_size)*2`.
#' 7. Mask to the template footprint, set the final layer name, optionally plot initial gaps and/or result,
#'    and write with LZW/tiling (adds a `PREDICTOR` appropriate to datatype if missing).
#' 8. Return a single-row **data.frame** with path, method, gap counts, max gap, filter size, and elapsed time.
#'
#' **Missing value workflows** (`missing_job`):
#' - `"CoverInput"`: cover the **input** before resampling (with `input_bg`).
#' - `"CoverOutput"`: cover the **resampled** raster (with `output_bg`).
#' - `"FillOutput"`: if gaps remain after resampling, compute a max-gap distance
#'   (via `terra::distance`) and call `whitebox::wbt_fill_missing_data()` with a filter
#'   width approx. twice the maximum gap (in pixels). Only a minimum clamp (`>= 3`) is applied.
#' - Combinations: `"CoverInput-CoverOutput"`, `"CoverInput-FillOutput"`.
#' - `"none"`: do nothing about gaps.
#'
#' **Resampling**:
#' - Default `summary_function = "average"` (area-preserving). `"mean"` is mapped to `"average"`.
#' - If `is_categorical = TRUE`, method is forced to `"near"`.
#'
#' **I/O and stability**:
#' - Writes are tiled/ threaded with `LZW` compression. Adds `PREDICTOR=2` (floats) or
#'   `PREDICTOR=3` (ints) if not provided.
#' - `terraOptions(memfrac, tempdir, todisk)` are set for the call and **restored** on exit.
#' - Console messages use `cat()` via a `say()` helper and the call is **sink-safe**.
#'
#' @param input `terra::SpatRaster` or character path (fine resolution input).
#' @param egv_template `terra::SpatRaster` or path (coarser EGV template).
#' @param summary_function Character resample method. Default `"average"`; `"mean"` is mapped to `"average"`.
#'   If `is_categorical = TRUE`, `"near"` is enforced.
#' @param missing_job One of `"CoverInput"`, `"CoverOutput"`, `"FillOutput"`,
#'   `"CoverInput-CoverOutput"`, `"CoverInput-FillOutput"`, `"none"`. Default `"none"`.
#' @param input_bg Background for CoverInput: `SpatRaster`, path, or **numeric** (+ `input_template`).
#' @param input_template Template used when `input_bg` is numeric (builds a background raster).
#' @param output_bg Background for CoverOutput: `SpatRaster`, path, or **numeric** (+ `egv_template`).
#' @param idw_weight Numeric power for Whitebox IDW. Default `2`.
#' @param is_categorical Logical; if `TRUE`, use `"near"` and skip WBT fill. Default `FALSE`.
#' @param check_alignment Logical; if `TRUE`, cheaply checks CRS/res; warns if mismatched. Default `TRUE`.
#' @param outlocation Output directory. Default `"./Rastri_100m/RAW/"`.
#' @param outfilename Output filename (e.g., `"layer.tif"`). **Required**.
#' @param layername Layer name to assign before writing. **Required**.
#' @param plot_gaps Logical; plot the **initial** gap map (guarded by `interactive()`). Default `FALSE`.
#' @param plot_final Logical; plot the final raster (guarded by `interactive()`). Default `FALSE`.
#' @param NAflag Optional numeric NA flag for writing. Default `NULL` (terra default).
#' @param gdal_opts Character vector of GDAL creation options (merged with tuned defaults). Default
#'   `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Optional terra datatype (e.g., `"FLT4S"`, `"INT2S"`). Default `NULL`.
#' @param terra_memfrac Memory fraction for `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra/Whitebox operations. Default: `tempdir()`.
#' @param terra_todisk Logical (set `terraOptions(todisk=...)` for this run). Default `TRUE`.
#' @param force_gc Logical; if `TRUE`, call `gc()` at checkpoints. Default `FALSE`.
#' @param return_visible Logical; if `TRUE`, return data.frame is visible; otherwise invisible. Default `FALSE`.
#' @param quiet Logical; suppress console messages. Default `FALSE`.
#'
#' @return A **single-row `data.frame`** with:
#'   - `path` (character): written file path
#'   - `method` (character): resample method used
#'   - `missing_job` (character)
#'   - `n_gaps_initial` (integer; after resample/CoverInput, before output-stage ops)
#'   - `n_gaps_final` (integer; after CoverOutput/FillOutput)
#'   - `max_gap_dist` (numeric; `NA` if not computed)
#'   - `filter_w` (integer; `NA` if not used)
#'   - `elapsed_sec` (numeric)
#'
#' @examples
#' \dontrun{
#' input2egv(
#'   input            = "./TestejuPakotni_early/Forests_StandAge_bg0.tif",
#'   egv_template     = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   summary_function = "average",
#'   missing_job      = "CoverInput-FillOutput",
#'   input_bg         = 0,                     # numeric background + input_template
#'   input_template   = "./Templates/TemplateRasters/LV10m_10km.tif",
#'   outlocation      = "./Rastri_100m/RAW/",
#'   outfilename      = "StandAge_100m.tif",
#'   layername        = "stand_age",
#'   force_gc         = TRUE,
#'   quiet            = FALSE
#' )
#' }
#'
#' @seealso [polygon2input()], [downscale2egv()], [distance2egv()]
#' @import terra
#' @importFrom whitebox wbt_fill_missing_data
#' @export
input2egv <- function(
    input,
    egv_template,
    summary_function   = "average",
    missing_job        = c("CoverInput","CoverOutput","FillOutput",
                           "CoverInput-CoverOutput","CoverInput-FillOutput","none"),
    input_bg           = NULL,
    input_template     = NULL,
    output_bg          = NULL,
    idw_weight         = 2,
    is_categorical     = FALSE,
    check_alignment    = TRUE,
    outlocation        = "./Rastri_100m/RAW/",
    outfilename,
    layername,
    plot_gaps          = FALSE,
    plot_final         = FALSE,
    NAflag             = NULL,
    gdal_opts          = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER",
                           "NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256"),
    write_datatype     = NULL,
    terra_memfrac      = 0.7,
    terra_tempdir      = tempdir(),
    terra_todisk       = TRUE,
    force_gc           = FALSE,
    return_visible     = FALSE,
    quiet              = FALSE
) {
  t0 <- proc.time()[["elapsed"]]

  # ---- sink safety ----
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
  .need_pkg("terra", "resampling / aggregation / writing")
  .need_pkg("fs",    "filesystem ops")
  .need_pkg("whitebox", "gap filling")

  # helpers
  say <- function(...) if (!quiet) cat(..., "\n")
  .maybe_gc <- function() if (isTRUE(force_gc)) gc()
  .ensure_dir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  .as_rast <- function(x) {
    if (inherits(x, "SpatRaster")) return(x)
    if (is.character(x)) return(terra::rast(x))
    stop("Provide a SpatRaster or a filepath for argument: ", deparse(substitute(x)))
  }
  .bg_from <- function(value_or_rast_or_path, template_for_value) {
    if (inherits(value_or_rast_or_path, "SpatRaster")) return(value_or_rast_or_path)
    if (is.character(value_or_rast_or_path)) return(terra::rast(value_or_rast_or_path))
    if (is.numeric(value_or_rast_or_path) && length(value_or_rast_or_path) == 1) {
      if (is.null(template_for_value)) stop("Numeric background provided but corresponding template is missing.")
      tpl <- .as_rast(template_for_value)
      return(terra::mask(terra::init(tpl, value_or_rast_or_path), tpl))
    }
    stop("Unsupported background type. Use SpatRaster, path, or numeric + template.")
  }
  .map_method <- function(m, is_cat) {
    if (is_cat) return("near")
    m <- tolower(m)
    if (m == "mean") return("average")
    m
  }

  # args & setup
  missing_job <- match.arg(missing_job)
  if (missing(outfilename)) stop("'outfilename' is required.")
  if (missing(layername))   stop("'layername' is required.")
  .ensure_dir(outlocation)

  # terra options (silenced)
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac, tempdir = old_opt$tempdir, todisk = old_opt$todisk, progress = old_opt$progress),
          add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir, todisk = isTRUE(terra_todisk), progress = FALSE)
  })

  # load
  input_r    <- .as_rast(input)
  template_r <- .as_rast(egv_template)

  # cheap alignment checks
  if (isTRUE(check_alignment)) {
    if (!terra::same.crs(input_r, template_r)) warning("CRS mismatch between input and template.")
    r_in  <- terra::res(input_r); r_tpl <- terra::res(template_r)
    if (any(is.na(r_in)) || any(is.na(r_tpl))) warning("Unable to read raster resolution for alignment check.")
  }

  # CoverInput
  if (missing_job %in% c("CoverInput", "CoverInput-CoverOutput", "CoverInput-FillOutput")) {
    if (is.null(input_bg)) stop("missing_job requests 'CoverInput' but 'input_bg' is NULL.")
    input_bg_r <- .bg_from(input_bg, template_for_value = input_template)
    input_r <- terra::cover(input_r, input_bg_r)
  }

  # resample
  method <- .map_method(summary_function, is_categorical)
  tmp1   <- tempfile(pattern = "input2egv_align_", tmpdir = terra_tempdir, fileext = ".tif")
  on.exit({ if (file.exists(tmp1)) unlink(tmp1) }, add = TRUE)
  aligned_r <- terra::resample(input_r, template_r, method = method, filename = tmp1, overwrite = TRUE)
  rm(input_r); .maybe_gc()

  # initial gaps (inside template)
  gaps_initial   <- is.na(aligned_r) & !is.na(template_r)
  n_gaps_initial <- as.integer(terra::global(terra::as.int(gaps_initial), fun = "sum", na.rm = TRUE))
  if (is.na(n_gaps_initial)) n_gaps_initial <- 0L

  if (plot_gaps && interactive()) {
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar), add = TRUE)
    if (plot_final) par(mfrow = c(1, 2))
    plot(gaps_initial, main = "Initial gaps (after resample/CoverInput)")
  }

  # CoverOutput
  if (missing_job %in% c("CoverOutput", "CoverInput-CoverOutput")) {
    if (is.null(output_bg)) stop("missing_job requests 'CoverOutput' but 'output_bg' is NULL.")
    output_bg_r <- .bg_from(output_bg, template_for_value = template_r)
    aligned_r   <- terra::cover(aligned_r, output_bg_r)
  }

  # FillOutput
  filter_w <- NA_integer_
  max_gap_dist <- NA_real_
  if (!is_categorical && missing_job %in% c("FillOutput","CoverInput-FillOutput")) {
    if (n_gaps_initial > 0L) {
      fillable <- terra::ifel(!gaps_initial, 1, NA)
      dist_r   <- terra::distance(fillable)
      max_gap_dist <- as.numeric(terra::global(dist_r, fun = "max", na.rm = TRUE))
      rm(fillable, dist_r); .maybe_gc()

      if (!is.na(max_gap_dist) && max_gap_dist > 0) {
        pix_size <- mean(terra::res(template_r))
        filter_w <- max(3L, as.integer(ceiling(max_gap_dist / pix_size) * 2L))  # <- only MIN clamp now
        if (!is.na(filter_w) && filter_w >= 3L) {
          tmp2 <- tempfile(pattern = "input2egv_fill_", tmpdir = terra_tempdir, fileext = ".tif")
          on.exit({ if (file.exists(tmp2)) unlink(tmp2) }, add = TRUE)
          if (!file.exists(tmp1)) terra::writeRaster(aligned_r, filename = tmp1, overwrite = TRUE)
          ok <- TRUE
          tryCatch({
            whitebox::wbt_fill_missing_data(
              input    = tmp1,
              output   = tmp2,
              filter   = filter_w,
              weight   = idw_weight,
              no_edges = FALSE
            )
          }, error = function(e) {
            say(sprintf("Whitebox filling failed: %s", conditionMessage(e)))
            ok <<- FALSE
          })
          if (ok && file.exists(tmp2)) {
            aligned_r <- terra::rast(tmp2)
          } else {
            say("Proceeding without WBT fill due to error or missing output.")
          }
        }
      }
    }
  }
  rm(gaps_initial); .maybe_gc()

  # final gaps (before masking/writing)
  gaps_final   <- is.na(aligned_r) & !is.na(template_r)
  n_gaps_final <- as.integer(terra::global(terra::as.int(gaps_final), fun = "sum", na.rm = TRUE))
  if (is.na(n_gaps_final)) n_gaps_final <- 0L
  rm(gaps_final); .maybe_gc()

  # mask & name
  final_r <- terra::mask(aligned_r, template_r)
  names(final_r) <- layername
  rm(aligned_r); .maybe_gc()

  if (plot_final && interactive()) {
    plot(final_r, main = paste0("Final: ", layername))
  }

  # GDAL opts
  def_opts <- c("NUM_THREADS=ALL_CPUS","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256","COMPRESS=LZW","BIGTIFF=IF_SAFER")
  gdal_opts <- unique(c(gdal_opts, def_opts))
  if (!any(grepl("^PREDICTOR=", gdal_opts))) {
    pred <- if (!is.null(write_datatype) && grepl("^(U?INT)", toupper(write_datatype))) "PREDICTOR=3" else "PREDICTOR=2"
    gdal_opts <- c(gdal_opts, pred)
  }

  # write (atomic)
  out_path <- file.path(outlocation, outfilename)
  ext  <- tools::file_ext(out_path); ext <- if (nzchar(ext)) ext else "tif"
  stem <- sub(sprintf("\\.%s$", ext), "", out_path)
  tmp_out <- sprintf("%s._tmp.%s", stem, ext)
  if (file.exists(tmp_out)) try(unlink(tmp_out), silent = TRUE)

  terra::crs(final_r) <- terra::crs(template_r)  # WKT/EPSG (not proj4)
  write_args <- list(x = final_r, filename = tmp_out, overwrite = TRUE, gdal = gdal_opts)
  if (!is.null(write_datatype)) write_args$datatype <- write_datatype
  if (!is.null(NAflag))        write_args$NAflag   <- NAflag
  do.call(terra::writeRaster, write_args)

  if (file.exists(out_path)) try(unlink(out_path), silent = TRUE)
  file.rename(tmp_out, out_path)

  rm(final_r); .maybe_gc()

  # result as data.frame
  elapsed_sec <- as.numeric(proc.time()[["elapsed"]] - t0)
  res_df <- data.frame(
    path            = out_path,
    method          = method,
    missing_job     = missing_job,
    n_gaps_initial  = n_gaps_initial,
    n_gaps_final    = n_gaps_final,
    max_gap_dist    = max_gap_dist,
    filter_w        = filter_w,
    elapsed_sec     = elapsed_sec,
    stringsAsFactors = FALSE
  )
  if (return_visible) return(res_df) else return(invisible(res_df))
}
