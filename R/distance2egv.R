#' Distance-to-class on the EGV grid
#'
#' @description
#' Computes Euclidean distance (in **map units**) from cells matching a set of
#' class values in an **input raster** to all cells of an **EGV template** grid,
#' then writes a Float32 GeoTIFF aligned to the template. Designed to work with
#' rasters produced by `polygon2input()`.
#'
#' @details
#' **Preparation & CRS:** Prepare inputs with the function `polygon2input()` and the
#' templates from https://zenodo.org/records/14497070 (download via the function
#' `download_raster_templates()`). The function does **not** automatically reproject
#' to the 100 m template. It only warns if CRSs differ. For rare cases, set
#' `project_to_template_input = TRUE` and provide `template_input` (10 m) to reproject
#' the input once at the start. Distance is computed in the (possibly reprojected)
#' inputs map units.
#'
#' **Selecting class values:** `values_as_one` accepts any-length vector combining:
#' - numeric values (e.g., \code{500}, \code{c(610,620,630)}), and/or
#' - range strings in interval notation with inclusive/exclusive bounds:
#'   \code{"[400,600]"}, \code{"(400,600)"}, \code{"[400,600)"}, \code{"(400,600]"}.
#'   Example: \code{c("[600,700)", "500")}. If `values_as_one` is `NULL`, **non-NA**
#'   cells are considered sources.
#'
#' **Distance engine:** Set `use_whitebox = TRUE` (default) to compute the distance with
#' `whitebox::wbt_euclidean_distance()`; otherwise `terra::distance()` is used.
#'
#' **Auto alignment choice:** If input distance grid and the 100 m template are
#' same-CRS, bounding boxes match, and the template resolution is an **integer multiple**
#' of the input resolution (e.g., 10 m to 100 m), the function uses
#' `terra::aggregate(..., fun=mean)` (fast), followed by a light
#' `resample(..., "near")` to lock onto the template grid. Otherwise, it falls back
#' to `terra::resample(..., method="mean")`.
#'
#' **Masking:** A single mask to `template_egv` is applied once **after** alignment
#' and **before** plotting/saving.
#'
#' **Gap filling:** If `fill_gaps = TRUE`, gaps (cells where output is `NA` but the
#' template is not) are filled via `whitebox::wbt_fill_missing_data()` with IDW
#' weight `idw_weight`. If `filter_size` is `NULL`, the function uses
#' `2 * max_gap_width_cells`, with a **minimum of 3** (no max clamp). Whitebox
#' temporaries are written with `COMPRESS=NONE` for speed; the final GeoTIFF uses `gdal_opts`.
#'
#' **Workflow**
#' 1. Load `input` and `template_egv`; optionally reproject `input` once if
#'    `project_to_template_input=TRUE` (using `template_input`).
#' 2. Build a **seeds** raster on the input grid: `1` where value matches `values_as_one`
#'    (or non-`NA` if `values_as_one=NULL`), `NA` elsewhere.
#' 3. Compute **distance** on the input grid via Whitebox (`wbt_euclidean_distance`) or
#'    `terra::distance()` depending on `use_whitebox`.
#' 4. **Align** the distance raster to the template:
#'    - If perfectly nested (same CRS, matching extent, integer resolution ratio),
#'      do `aggregate(mean)` then `resample("near")`.
#'    - Else if same CRS but not nested, `resample(method="mean")`.
#'    - Else, `project(..., method="bilinear")`.
#' 5. Apply a **single final mask** to `template_egv`.
#' 6. Optionally **fill gaps** with Whitebox (if `fill_gaps=TRUE`) on the template grid.
#' 7. Optionally **plot** result and/or gap map (side-by-side if both requested).
#' 8. **Write atomically** to GeoTIFF with LZW compression and tuned GDAL options.
#' 9. Restore `terraOptions()` and sink state on exit (prevents stuck sinks).
#'
#' Console safety: uses `cat()` for progress and snapshots/restores sink state on exit
#' so your console will not remain "sunk" after interrupts.
#'
#' @param input `terra::SpatRaster` or file path. Raster prepared by `polygon2input()`.
#' @param template_egv `terra::SpatRaster` or path. Target EGV template (e.g., 100 m).
#' @param values_as_one `NULL`, numeric vector, and/or range strings in bracket
#'   notation (see Details). If `NULL`, non-`NA` cells are sources.
#' @param project_to_template_input Logical. If `TRUE`, reproject `input` to
#'   `template_input` **once at the start**. Default `FALSE`.
#' @param template_input When `project_to_template_input=TRUE`, a 10 m template
#'   used to define the target CRS/grid for that initial reprojection.
#' @param use_whitebox Logical; use Whitebox for distance (default `TRUE`),
#'   otherwise `terra::distance()`.
#' @param fill_gaps Logical. If `TRUE`, fill remaining gaps on the **template grid**
#'   using Whitebox `wbt_fill_missing_data()`. Default `FALSE`.
#' @param idw_weight IDW power for Whitebox fill. Default `2`.
#' @param filter_size Optional odd integer window for Whitebox. If `NULL`, auto:
#'   `2 * max_gap_width_cells`, with a **minimum of 3**.
#' @param outlocation Output directory. Default `"./Rastri_100m/RAW/"`.
#' @param outfilename Output filename (e.g., `"dist_class_100m.tif"`). **Required**.
#' @param layername Output band name. **Required**.
#' @param check_na Logical. If `TRUE`, report internal NA count on the template footprint. Default `FALSE`.
#' @param plot_result Logical; plot final result.
#' @param plot_gaps Logical; plot gap map (TRUE = gap). If both TRUE, side-by-side.
#' @param NAflag Optional NA flag for writing. Default `NULL` (omitted for Float32).
#' @param gdal_opts GDAL creation options (merged with tuned defaults).
#'   Default `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Terra datatype for writing. Default `NULL` will be coded as `"FLT4S"`.
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra ops. Default `tempdir()`.
#' @param terra_todisk Logical or `NA`. If `TRUE`, prefer on-disk for heavy steps. Default `FALSE`.
#' @param force_gc Logical; call `gc()` at checkpoints. Default `FALSE`.
#' @param quiet Suppress console prints. Default `FALSE`.
#'
#' @return Invisibly, a list with:
#'   - `path` (output file),
#'   - `n_sources` (number of source cells on input grid),
#'   - `n_na_final` (internal NA count on template footprint),
#'   - `min_dist`, `max_dist`,
#'   - `elapsed_sec`.
#'
#' @examples
#' \dontrun{
#' # Distance from classes 1 and [10,20) using Whitebox, with side-by-side plots
#' distance2egv(
#'   input         = "./TestejuPakotni_early/Forests_StandAge_bg0_b.tif",
#'   template_egv  = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   values_as_one = c(1, "[10,20)"),
#'   outlocation   = "./Rastri_100m/RAW/",
#'   outfilename   = "dist_age_100m.tif",
#'   layername     = "dist_age",
#'   use_whitebox  = TRUE,
#'   plot_result   = TRUE,
#'   plot_gaps     = TRUE,
#'   terra_todisk  = TRUE
#' )
#'
#' # Same with terra::distance() and gap filling
#' distance2egv(
#'   input         = "./Masks/roads_mask_10m.tif",
#'   template_egv  = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   values_as_one = NULL,  # non-NA are sources
#'   use_whitebox  = FALSE,
#'   fill_gaps     = TRUE,
#'   idw_weight    = 2,
#'   outlocation   = "./Rastri_100m/RAW/",
#'   outfilename   = "dist_roads_100m.tif",
#'   layername     = "dist_roads",
#'   terra_todisk  = TRUE
#' )
#' }
#'
#' @seealso
#'   [polygon2input()] for preparing the input raster,
#'   [input2egv()] for aggregation/resampling to the EGV grid,
#'   [downscale2egv()], [landscape_function()], [radius_function()],
#'   and the template downloaders: [download_raster_templates()], [download_vector_templates()].
#'
#' @importFrom utils capture.output
#' @importFrom graphics par
#' @importFrom terra rast crs same.crs res ext nrow ncol mask project aggregate resample ifel distance global writeRaster as.int ncell app
#' @importFrom whitebox wbt_euclidean_distance wbt_fill_missing_data
#' @export
distance2egv <- function(
    input,
    template_egv,
    values_as_one = NULL,
    project_to_template_input = FALSE,
    template_input = NULL,
    use_whitebox  = TRUE,
    fill_gaps     = FALSE,
    idw_weight    = 2,
    filter_size   = NULL,
    outlocation   = "./Rastri_100m/RAW/",
    outfilename,
    layername,
    check_na      = FALSE,
    plot_result   = FALSE,
    plot_gaps     = FALSE,
    NAflag        = NULL,
    gdal_opts     = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256"),
    write_datatype = NULL,   # default to FLT4S
    terra_memfrac  = 0.7,
    terra_tempdir  = tempdir(),
    terra_todisk   = FALSE,
    force_gc       = FALSE,
    quiet          = FALSE
) {
  t0 <- Sys.time()

  # ---- sink safety ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)

  say <- function(...) if (!quiet) cat(..., "\n")
  .maybe_gc <- function() if (isTRUE(force_gc)) gc()
  .as_rast <- function(x) if (inherits(x, "SpatRaster")) x else terra::rast(x)

  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("terra", "distance / resample / write")
  .need_pkg("whitebox", "Euclidean distance and/or gap fill via Whitebox")


  # robust range parser: "[a,b]", "(a,b]", "[a,)", "(-inf,10)", "[1,+inf]" etc.
  in_range_str <- function(x, spec) {
    s <- trimws(as.character(spec)[1])
    s <- trimws(s)
    if (nchar(s) < 3) stop("Bad range spec: ", spec)

    lb <- substr(s, 1, 1)
    rb <- substr(s, nchar(s), nchar(s))
    if (!(lb %in% c("(", "[")) || !(rb %in% c(")", "]"))) stop("Bad range spec: ", spec)

    inner <- substr(s, 2, nchar(s) - 1)
    # split on first comma only
    comma_pos <- regexpr(",", inner, fixed = TRUE)[1]
    if (comma_pos < 1) stop("Bad range spec: ", spec)

    a_txt <- trimws(substr(inner, 1, comma_pos - 1))
    b_txt <- trimws(substr(inner, comma_pos + 1, nchar(inner)))

    parse_bound <- function(z) {
      if (!nzchar(z)) return(NA_real_)  # unbounded
      zl <- tolower(z)
      if (zl %in% c("inf","+inf","infinity","+infinity","*")) return( Inf)
      if (zl %in% c("-inf","-infinity"))                   return(-Inf)
      as.numeric(z)
    }

    a <- parse_bound(a_txt)
    b <- parse_bound(b_txt)

    # Lower test
    lower_ok <- if (is.na(a) || (is.infinite(a) && a < 0)) {
      TRUE
    } else if (lb == "(") {
      x > a
    } else {
      x >= a
    }

    # Upper test
    upper_ok <- if (is.na(b) || (is.infinite(b) && b > 0)) {
      TRUE
    } else if (rb == ")") {
      x < b
    } else {
      x <= b
    }

    lower_ok & upper_ok
  }


  # ---- args & I/O ----
  if (missing(outfilename)) stop("'outfilename' is required.")
  if (missing(layername))   stop("'layername' is required.")
  if (!dir.exists(outlocation)) dir.create(outlocation, recursive = TRUE, showWarnings = FALSE)

  # ---- terra tuning (silenced) + restore on exit ----
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac, tempdir = old_opt$tempdir, todisk = old_opt$todisk, progress = old_opt$progress),
          add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
    terra::terraOptions(progress = FALSE)
  })

  # ---- load rasters ----
  r_in <- .as_rast(input)
  tmpl <- .as_rast(template_egv)

  # optional early reprojection (once), if requested
  if (isTRUE(project_to_template_input)) {
    if (is.null(template_input)) stop("Set 'template_input' (10 m) when 'project_to_template_input=TRUE'.")
    tpl_in <- .as_rast(template_input)
    say("Reprojecting input to the 10 m template CRS/grid (once) ...")
    r_in <- terra::project(r_in, tpl_in, method = "near",
                           filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_proj_", basename(tempfile("")), ".tif")) else "",
                           overwrite = TRUE)
  } else if (!terra::same.crs(r_in, tmpl)) {
    warning("CRS of input and template_egv differ. Distance will be computed on the input CRS/grid.\n",
            "Consider 'project_to_template_input=TRUE' with 'template_input' to reproject once at the start.")
  }

  # ---- build source (seed) raster on the input grid ----
  if (is.null(values_as_one)) {
    seeds <- terra::ifel(!is.na(r_in), 1, NA,
                         filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_seeds_", basename(tempfile("")), ".tif")) else "",
                         overwrite = TRUE)
  } else {
    nums_all <- suppressWarnings(as.numeric(values_as_one))
    is_num   <- !is.na(nums_all)
    nums     <- values_as_one[is_num]
    rngs     <- values_as_one[!is_num]
    seeds <- terra::app(
      r_in,
      fun = function(x) {
        keep <- rep(FALSE, length(x))
        if (length(nums)) {
          nu <- suppressWarnings(as.numeric(nums))
          keep <- keep | (x %in% nu)
        }
        if (length(rngs)) {
          for (sp in rngs) {
            spn <- suppressWarnings(as.numeric(sp))
            if (!is.na(spn)) {
              keep <- keep | (x == spn)
            } else {
              z <- in_range_str(x, sp)
              keep <- keep | (isTRUE(z) | (!is.na(z) & z))
            }
          }
        }
        y <- rep(NA_real_, length(x)); y[keep] <- 1; y
      },
      filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_seeds_", basename(tempfile("")), ".tif")) else "",
      overwrite = TRUE
    )
  }

  # tally sources
  n_sources <- as.integer(terra::global(terra::as.int(!is.na(seeds)), fun = "sum", na.rm = TRUE))
  if (is.na(n_sources)) n_sources <- 0L
  if (n_sources == 0L) say("Warning: no source cells found (no class matches).")

  .maybe_gc()

  # ---- distance on the input grid ----
  say(if (use_whitebox) "Computing distance with WhiteboxTools ..." else "Computing distance with terra::distance() ...")
  if (isTRUE(use_whitebox)) {
    wb_in  <- file.path(terra_tempdir, paste0("d2e_wb_in_",  basename(tempfile("")), ".tif"))
    wb_out <- file.path(terra_tempdir, paste0("d2e_wb_out_", basename(tempfile("")), ".tif"))
    bin <- terra::ifel(!is.na(seeds), 1, 0)
    terra::writeRaster(bin, wb_in, overwrite = TRUE,
                       gdal = c("COMPRESS=NONE","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256"))
    rm(bin); .maybe_gc()

    ok <- TRUE
    tryCatch({
      whitebox::wbt_euclidean_distance(input = wb_in, output = wb_out)
    }, error = function(e) {
      say(sprintf("Whitebox euclidean_distance failed: %s ; falling back to terra::distance().", conditionMessage(e)))
      ok <<- FALSE
    })
    d_in <- if (ok && file.exists(wb_out)) terra::rast(wb_out) else {
      terra::distance(seeds,
                      filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_dist_", basename(tempfile("")), ".tif")) else "",
                      overwrite = TRUE)
    }
  } else {
    d_in <- terra::distance(seeds,
                            filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_dist_", basename(tempfile("")), ".tif")) else "",
                            overwrite = TRUE)
  }
  rm(seeds); .maybe_gc()

  # ---- align to template (aggregate with near if perfectly nested; else resample mean) ----
  same_crs <- terra::same.crs(d_in, tmpl)
  bbox_match <- all(terra::ext(d_in) == terra::ext(tmpl))
  res_in  <- terra::res(d_in);  res_t  <- terra::res(tmpl)
  mult_ok <- all(abs(res_t / res_in - round(res_t / res_in)) < 1e-7)
  aligned <- NULL

  if (same_crs && bbox_match && mult_ok) {
    fact <- as.integer(round(res_t / res_in))
    fact <- ifelse(is.finite(fact), fact, 1L)
    say(sprintf("Align: aggregate (mean) by %dx%dx, then resample (near).", fact[1], fact[2]))
    a1 <- terra::aggregate(d_in, fact = fact, fun = "mean",
                           filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_aggr_", basename(tempfile("")), ".tif")) else "",
                           overwrite = TRUE)
    aligned <- terra::resample(a1, tmpl, method = "near",
                               filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_align_", basename(tempfile("")), ".tif")) else "",
                               overwrite = TRUE)
    rm(a1)
  } else if (same_crs) {
    say("Align: resample(mean) to template grid.")
    aligned <- terra::resample(d_in, tmpl, method = "mean",
                               filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_align_", basename(tempfile("")), ".tif")) else "",
                               overwrite = TRUE)
  } else {
    say("CRS differs: projecting distance raster to template CRS (bilinear) for output alignment.")
    aligned <- terra::project(d_in, tmpl, method = "bilinear",
                              filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_proj_", basename(tempfile("")), ".tif")) else "",
                              overwrite = TRUE)
  }
  rm(d_in); .maybe_gc()

  # ---- final mask to template ----
  aligned <- terra::mask(aligned, tmpl,
                         filename = if (isTRUE(terra_todisk)) file.path(terra_tempdir, paste0("d2e_mask_", basename(tempfile("")), ".tif")) else "",
                         overwrite = TRUE)
  names(aligned) <- layername

  # ---- optional gap fill on template grid (Whitebox) ----
  if (isTRUE(fill_gaps)) {
    if (is.null(filter_size)) {
      gaps    <- is.na(aligned) & !is.na(tmpl)
      fillable <- terra::ifel(!gaps, 1, NA)
      gdist    <- terra::distance(fillable)
      pix      <- mean(terra::res(tmpl))
      max_gap  <- suppressWarnings(terra::global(gdist, fun = "max", na.rm = TRUE)[[1]])
      rm(fillable, gdist); .maybe_gc()
      if (is.finite(max_gap) && !is.na(max_gap)) {
        fs <- as.integer(ceiling(max_gap / pix) * 2L)
        if (fs < 3L) fs <- 3L
        if (fs %% 2L == 0L) fs <- fs + 1L
        filter_size <- fs
      } else {
        filter_size <- 3L
      }
    }
    in_wbt  <- file.path(terra_tempdir, paste0("d2e_wbt_in_", basename(tempfile("")), ".tif"))
    out_wbt <- file.path(terra_tempdir, paste0("d2e_wbt_out_", basename(tempfile("")), ".tif"))
    terra::writeRaster(aligned, in_wbt, overwrite = TRUE,
                       gdal = c("COMPRESS=NONE","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256"))
    ok <- TRUE
    say(sprintf("Whitebox fill: filter=%d, weight=%.2f", as.integer(filter_size), idw_weight))
    tryCatch({
      whitebox::wbt_fill_missing_data(
        input    = in_wbt,
        output   = out_wbt,
        filter   = as.integer(filter_size),
        weight   = idw_weight,
        no_edges = FALSE
      )
    }, error = function(e) {
      say(sprintf("Whitebox failed: %s (continuing without fill)", conditionMessage(e)))
      ok <<- FALSE
    })
    if (ok && file.exists(out_wbt)) {
      aligned <- terra::rast(out_wbt)
      names(aligned) <- layername
    }
  }

  # ---- NA count (inside template footprint only) ----
  n_na_final <- NA_integer_
  if (isTRUE(check_na) || isTRUE(plot_gaps)) {
    gaps <- is.na(aligned) & !is.na(tmpl)
    n_na_final <- as.integer(terra::global(terra::as.int(gaps), fun = "sum", na.rm = TRUE))
    if (is.na(n_na_final)) n_na_final <- 0L
  }

  # ---- stats & plots ----
  min_dist <- suppressWarnings(terra::global(aligned, fun = "min", na.rm = TRUE)[[1]])
  max_dist <- suppressWarnings(terra::global(aligned, fun = "max", na.rm = TRUE)[[1]])
  if ((isTRUE(plot_result) || isTRUE(plot_gaps)) && interactive()) {
    oldpar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(oldpar), add = TRUE)
    if (isTRUE(plot_result) && isTRUE(plot_gaps)) graphics::par(mfrow = c(1, 2))
    if (isTRUE(plot_result)) terra::plot(aligned, main = paste0("distance2egv: ", layername))
    if (isTRUE(plot_gaps)) {
      gap_r <- terra::ifel(is.na(aligned) & !is.na(tmpl), 1, NA)
      terra::plot(gap_r, main = "Gaps (1 = NA inside template)", col = c(NA, "red"), legend = FALSE)
    }
  }

  # ---- write (atomic, LZW, Float32 default) ----
  def_opts <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")
  gdal_opts <- unique(c(gdal_opts, def_opts))
  if (is.null(write_datatype)) write_datatype <- "FLT4S"
  if (!any(grepl("^PREDICTOR=", gdal_opts))) gdal_opts <- c(gdal_opts, "PREDICTOR=2")

  out_path <- file.path(outlocation, outfilename)
  ext  <- tools::file_ext(out_path); ext <- if (nzchar(ext)) ext else "tif"
  stem <- sub(sprintf("\\.%s$", ext), "", out_path)
  tmp_out <- sprintf("%s._tmp.%s", stem, ext)
  if (file.exists(tmp_out)) try(unlink(tmp_out), silent = TRUE)

  write_args <- list(x = aligned, filename = tmp_out, overwrite = TRUE, gdal = gdal_opts, datatype = write_datatype)
  if (!is.null(NAflag)) write_args$NAflag <- NAflag
  do.call(terra::writeRaster, write_args)

  if (file.exists(out_path)) try(unlink(out_path), silent = TRUE)
  file.rename(tmp_out, out_path)

  elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  invisible(list(
    path        = out_path,
    n_sources   = n_sources,
    n_na_final  = n_na_final,
    min_dist    = min_dist,
    max_dist    = max_dist,
    elapsed_sec = elapsed_sec
  ))
}
