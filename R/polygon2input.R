#' Rasterize polygons to a template grid, optionally restrict & cover gaps, then write GeoTIFF
#'
#' @description
#' Rasterizes polygon/multipolygon **sf** data to a raster aligned to a **template** GeoTIFF.
#' Rasterization targets a `raster::RasterLayer` built from the template (so grids normally match).
#' Projection is optional (`project_mode`). Missing values are counted **only** over valid template
#' cells. You may optionally restrict the result with a raster mask (`restrict_to`) using numeric
#' values or **bracketed range strings** (e.g., `"(0,5]"`, `"[10,)"`). Remaining `NA` cells can be
#' filled by covering with a background raster (`background_raster`) or a constant (`background_value`).
#' For large rasters, heavy steps (projection/mask/cover) can stream to disk via `terra_todisk=TRUE`.
#'
#' @details
#' **Workflow**
#' 1. (Optional) `prepare=TRUE` ensures valid polygons and transforms to the template CRS
#'    (`sf::st_make_valid()` + `sf::st_transform()`). Set `prepare=FALSE` if you guarantee this.
#' 2. Rasterize with `fasterize::fasterize()` onto a `raster::RasterLayer` built from the template.
#' 3. Depending on `project_mode`, **project** to the template grid:
#'    - `"auto"` (default): checks CRS/resolution/extent/rows/cols; projects **only if needed**.
#'    - `"never"`: skips projection.
#'    - `"always"`: forces projection (`"near"` for `value_type="categorical"`, `"bilinear"` for `"continuous"`).
#' 4. Apply a **template mask** (keeps template extent and `NA` footprint).
#' 5. (Optional) **Restrict** by `restrict_to` (path or in-memory `SpatRaster`):
#'    - If `restrict_values` is `NULL`, keep **non-NA** cells of `restrict_to`.
#'    - If supplied, keep cells matching **numbers** (scalar/vector) and/or **ranges** using bracket
#'      syntax: `"(a,b)"` open; `"[a,b]"` closed; mix ends like `"(a,b]"`. Use `-inf`/`+inf` for unbounded,
#'      e.g., `"[10,)"`, `"(-inf,0)"`. Multiple entries are OR-ed.
#' 6. (Optional) **Cover** remaining `NA` cells:
#'    - with `background_raster` (aligned to template; projected if needed), or
#'    - with a **constant** using `background_value` (fast and memory-light).
#' 7. A **final mask** to the template is applied **before plotting and saving**.
#' 8. NA counts are computed with `terra::global()` on raster masks (memory-safe).
#'
#' **Datatype & NAflag auto-chooser**
#' If you do **not** provide `write_datatype`/`NAflag`, the function picks:
#' - `value_type="categorical"` to `write_datatype="INT2S"`, `NAflag=-32768`
#' - `value_type="continuous"` to `write_datatype="FLT4S"`, `NAflag` omitted
#' You can always override via arguments.
#'
#' **Performance & stability**
#' - For very large rasters (e.g., ~1.3B cells), set `terra_todisk=TRUE` and consider a fast SSD for
#'   `terra_tempdir`. This streams big operations to disk and avoids "vector memory limit" errors.
#' - Projection is often unnecessary here because rasterization targets the template grid; `"auto"`
#'   will detect alignment and skip it.
#' - Writes are **atomic**: a temporary file is written and then moved into place.
#'
#' @param vector_data sf polygons/multipolygons to rasterize.
#' @param template_path Path to the template raster (.tif).
#' @param out_path Output directory. Default `"./"`.
#' @param file_name Output filename, e.g., `"my_input.tif"`.
#' @param value_field Character or `NULL`. Attribute to burn; if `NULL`, uses `constant_value`.
#' @param constant_value Numeric scalar burned when `value_field=NULL`. Default `1`.
#' @param fun Aggregation in `fasterize` for overlaps: one of `"max"`, `"sum"`, `"first"`, `"last"`, `"min"`, `"count"`.
#'   Default `"max"`.
#' @param value_type Guides resampling **if projection happens**:
#'   `"categorical"` with method `"near"`, `"continuous"` with `"bilinear"`. Default `"categorical"`.
#' @param project_mode `"auto"` (default; project only if misaligned), `"never"`, or `"always"`.
#' @param prepare Logical. If `TRUE`, run `st_make_valid()` and transform to template CRS. Default `TRUE`.
#' @param restrict_to Optional raster mask (path or `terra` `SpatRaster`).
#' @param restrict_values Optional values to keep from `restrict_to`: numbers (scalar/vector) and/or
#'   **range strings** using bracket inclusivity, e.g. `"(0,5]"`, `"[10,)"`, `"(-inf,0)"`. Multiple entries are OR-ed.
#'   If `NULL`, keep non-NA cells of `restrict_to`.
#' @param background_raster Optional path/`SpatRaster` used to cover remaining NAs.
#' @param background_value Numeric constant to fill where result is NA **if** no `background_raster`.
#'   Default `NULL` (no covering).
#' @param check_na Logical. If `TRUE`, report NA counts before/after. Default `TRUE`.
#' @param plot_result Logical. Plot final raster (after all processing). Default `FALSE`.
#' @param plot_gaps Logical. Plot NA gaps (only within template footprint). If both plot flags are `TRUE`,
#'   plots are side-by-side. Default `FALSE`.
#' @param NAflag Optional NA flag for writing (passed to GDAL). Default `NULL` (auto if needed).
#' @param gdal_opts Character vector of GDAL creation options (merged with tuned defaults).
#'   Default `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Optional terra datatype for writing (e.g., `"FLT4S"`, `"INT2S"`). Default `NULL` (auto).
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra operations. Default `tempdir()`.
#' @param terra_todisk Logical or `NA`. If `TRUE`, prefer on-disk processing for heavy ops (project/mask/cover).
#'   If `FALSE`, prefer memory. If `NA`, leave session default unchanged. Default `FALSE`.
#' @param force_gc Logical; call `gc()` at checkpoints. Default `FALSE`.
#' @param overwrite Overwrite output? Default `FALSE`.
#' @param quiet Suppress progress prints (`cat()`)? Default `FALSE`.
#'
#' @return Invisibly, a list with `out_file`, `n_cells`, `n_na_initial`, `n_na_final`, `elapsed_sec`, and `crs`.
#'
#' @section Range string syntax for restrict_values:
#' Use `(a,b)` for open interval, `[a,b]` for closed; mix ends like `(a,b]`.
#' Use `-inf`/`+inf` (or `inf`) for unbounded, e.g. `"[10,)"`, `"(-inf,0)"`.
#' Supply multiple strings to OR them, e.g. `c("(0,5]","[10,15)")`.
#'
#' @examples
#' \dontrun{
#' # Basic: burn constant "1", continuous output, no covering
#' polygon2input(
#' vector_data = my_polys_sf,
#' template_path = "./Templates/TemplateRasters/template10m.tif",
#' out_path = "./Outputs",
#' file_name = "mask_const1.tif",
#' value_field = NULL,
#' constant_value = 1,
#' value_type = "continuous",
#' project_mode = "auto",
#' prepare = FALSE,
#' check_na = TRUE,
#' plot_result = TRUE,
#' overwrite = TRUE
#' )
#'
#' # Restrict to classes 1 and 2, plus (10,20], then fill remaining NAs with 0
#' polygon2input(
#' vector_data = my_polys_sf,
#' template_path = "./Templates/TemplateRasters/template10m.tif",
#' out_path = "./Outputs",
#' file_name = "mask_restricted_bg0.tif",
#' value_field = "attr",
#' value_type = "categorical",
#' restrict_to = "./mask_classes.tif",
#' restrict_values = c(1, 2, "(10,20]"),
#' background_value = 0,
#' terra_todisk = TRUE, # stream big ops to disk
#' terra_tempdir = tempdir(), # or a fast SSD scratch
#' check_na = TRUE,
#' plot_result = TRUE,
#' plot_gaps = TRUE,
#' overwrite = TRUE
#' )
#' }
#'
#' @seealso [tiled_buffers()], [tile_vector_grid()]
#' @importFrom fs dir_exists dir_create file_exists file_delete file_move path_ext path_ext_remove
#' @importFrom sf st_make_valid st_crs st_transform
#' @importFrom terra rast crs ncell mask project cover writeRaster plot global ifel nrow ncol res ext app same.crs
#' @importFrom fasterize fasterize
#' @importFrom raster raster
#' @importFrom utils capture.output
#' @importFrom graphics par
#' @export
polygon2input <- function(
    vector_data,
    template_path,
    out_path = "./",
    file_name,
    value_field = NULL,
    constant_value = 1,
    fun = "max",
    value_type = c("categorical", "continuous"),
    project_mode = c("auto","never","always"),
    prepare = TRUE,
    restrict_to = NULL,
    restrict_values = NULL,
    background_raster = NULL,
    background_value = NULL,
    check_na = TRUE,
    plot_result = FALSE,
    plot_gaps = FALSE,
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
  t0 <- Sys.time()
  value_type   <- match.arg(value_type)
  project_mode <- match.arg(project_mode)

  # sink safety
  orig_out <- sink.number(); orig_msg <- sink.number(type = "message")
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
  .need_pkg("terra",       "raster operations")
  .need_pkg("sf",          "vector handling")
  .need_pkg("raster",      "building the fasterize target RasterLayer")
  .need_pkg("fasterize",   "fast polygon rasterization")
  .need_pkg("fs",          "filesystem ops (paths, exists)")


  # basic I/O
  if (missing(vector_data) || !inherits(vector_data, "sf")) stop("'vector_data' must be an sf object.")
  if (missing(template_path) || !fs::file_exists(template_path)) stop("'template_path' must be an existing .tif file.")
  if (missing(file_name) || !nzchar(file_name)) stop("'file_name' must be provided.")
  if (!fs::dir_exists(out_path)) fs::dir_create(out_path, recurse = TRUE)
  out_file <- file.path(out_path, file_name)
  if (fs::file_exists(out_file) && !overwrite) {
    if (!quiet) cat("Output exists. Set overwrite=TRUE to replace:", out_file, "\n")
    return(invisible(list(out_file = out_file, n_cells = NA_integer_,
                          n_na_initial = NA_integer_, n_na_final = NA_integer_,
                          elapsed_sec = NA_real_, crs = NA_character_)))
  }

  # terra options (silenced)
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac, tempdir = old_opt$tempdir, todisk = old_opt$todisk), add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
  })

  # helpers
  mktemp <- function(tag) file.path(terra_tempdir, sprintf("p2i_%s_%s.tif", tag, basename(tempfile(""))))
  aligned_tol <- function(a, b, eps = 1e-7) {
    if (terra::nrow(a) != terra::nrow(b) || terra::ncol(a) != terra::ncol(b)) return(FALSE)
    r1 <- terra::res(a); r2 <- terra::res(b)
    e1 <- terra::ext(a); e2 <- terra::ext(b)
    res_ok <- all(abs(r1 - r2) <= pmax(abs(r1), abs(r2), 1) * eps)
    ext_ok <- all(abs(c(e1$xmin, e1$xmax, e1$ymin, e1$ymax) - c(e2$xmin, e2$xmax, e2$ymin, e2$ymax)) <=
                    pmax(abs(c(e1$xmin, e1$xmax, e1$ymin, e1$ymax)),
                         abs(c(e2$xmin, e2$xmax, e2$ymin, e2$ymax)), 1) * eps)
    crs_ok <- identical(terra::crs(a, proj = TRUE), terra::crs(b, proj = TRUE))
    res_ok && ext_ok && crs_ok
  }
  in_range_str <- function(x, spec) {
    s <- trimws(spec)
    m <- regexec("^([\\(\\[])(.*),(.*)([\\)\\]])$", s); p <- regmatches(s, m)[[1]]
    if (length(p) != 5) stop("Bad range spec: ", spec)
    lb <- p[2]; a_txt <- trimws(p[3]); b_txt <- trimws(p[4]); rb <- p[5]
    to_num <- function(z) {
      zlow <- tolower(z)
      if (zlow %in% c("inf","+inf","infinity","+infinity")) return(Inf)
      if (zlow %in% c("-inf","-infinity")) return(-Inf)
      as.numeric(z)
    }
    a <- to_num(a_txt); b <- to_num(b_txt)
    gt <- if (lb == "(") x >  a else x >= a
    lt <- if (rb == ")") x <  b else x <= b
    gt & lt
  }

  # template & (optional) prepare
  tmpl <- terra::rast(template_path)
  tmpl_crs <- terra::crs(tmpl)             # <-- WKT/EPSG (fix)
  vec <- vector_data
  if (isTRUE(prepare)) {
    vec <- sf::st_make_valid(vec)
    if (!identical(sf::st_crs(vec)$wkt, tmpl_crs)) vec <- sf::st_transform(vec, tmpl_crs)
  }

  # rasterize onto RasterLayer built from the template
  tmpl_r <- raster::raster(tmpl)
  fld <- value_field
  if (is.null(fld)) {
    fld <- ".__const__"; vec[[fld]] <- constant_value
    on.exit({ if (fld %in% names(vec)) vec[[fld]] <- NULL }, add = TRUE)
  }
  if (!quiet) cat("Rasterizing polygons (fasterize) ...\n")
  r <- terra::rast(fasterize::fasterize(sf = vec, raster = tmpl_r, field = fld, fun = fun))
  terra::crs(r) <- tmpl_crs                  # <-- enforce WKT/EPSG on the raster
  if (isTRUE(terra_todisk)) r <- terra::writeRaster(r, filename = mktemp("rasterized"), overwrite = TRUE)

  # project if needed
  need_project <- switch(project_mode, "always" = TRUE, "never" = FALSE, "auto" = !aligned_tol(r, tmpl))
  if (need_project) {
    proj_method <- if (value_type == "continuous") "bilinear" else "near"
    if (!quiet) cat("Projecting to template grid with method =", proj_method, "...\n")
    r <- if (isTRUE(terra_todisk)) {
      terra::project(r, tmpl, method = proj_method, filename = mktemp("proj"), overwrite = TRUE)
    } else {
      terra::project(r, tmpl, method = proj_method)
    }
  } else if (!quiet) cat("Projection skipped (grids already aligned).\n")

  # keep template extent / NA mask
  r <- if (isTRUE(terra_todisk)) {
    terra::mask(r, tmpl, filename = mktemp("masktmpl"), overwrite = TRUE)
  } else {
    terra::mask(r, tmpl)
  }

  if (force_gc) gc()

  # NA count before restrict/cover (inside template)
  n_na_initial <- NA_integer_
  if (isTRUE(check_na)) {
    gaps_pre <- is.na(r) & !is.na(tmpl)
    n_na_initial <- terra::global(gaps_pre, fun = "sum", na.rm = TRUE)[[1]]
    if (!quiet) cat("Missing cells before restrict/cover:", n_na_initial, "\n")
  }

  # restriction
  if (!is.null(restrict_to)) {
    if (!quiet) cat("Applying restriction by values ...\n")
    rr <- if (inherits(restrict_to, "SpatRaster")) restrict_to else terra::rast(restrict_to)
    if (!aligned_tol(rr, tmpl)) rr <- terra::project(rr, tmpl, method = "near")
    mask_r <- if (is.null(restrict_values)) {
      terra::app(rr, fun = function(x){ y <- rep(NA_real_, length(x)); y[!is.na(x)] <- 1; y },
                 filename = if (isTRUE(terra_todisk)) mktemp("mask_keep") else "", overwrite = TRUE)
    } else {
      nums <- if (is.numeric(restrict_values)) restrict_values else numeric(0)
      rngs <- if (is.character(restrict_values)) restrict_values else character(0)
      terra::app(
        rr,
        fun = function(x) {
          keep <- rep(FALSE, length(x))
          if (length(nums)) keep <- keep | (x %in% nums)
          if (length(rngs)) {
            for (sp in rngs) {
              spn <- suppressWarnings(as.numeric(sp))
              if (!is.na(spn)) {
                keep <- keep | (x == spn)
              } else {
                z <- in_range_str(x, sp); keep <- keep | (isTRUE(z) | (!is.na(z) & z))
              }
            }
          }
          y <- rep(NA_real_, length(x)); y[keep] <- 1; y
        },
        filename = if (isTRUE(terra_todisk)) mktemp("mask_vals") else "", overwrite = TRUE
      )
    }
    r <- if (isTRUE(terra_todisk)) {
      terra::mask(r, mask_r, filename = mktemp("maskapply"), overwrite = TRUE)
    } else {
      terra::mask(r, mask_r)
    }
  }

  # background cover
  if (!is.null(background_raster) || !is.null(background_value)) {
    if (!quiet) cat("Covering remaining NAs with background ...\n")
    if (!is.null(background_raster)) {
      bg <- if (inherits(background_raster, "SpatRaster")) background_raster else terra::rast(background_raster)
      if (!aligned_tol(bg, tmpl)) {
        proj_method <- if (value_type == "continuous") "bilinear" else "near"
        bg <- terra::project(bg, tmpl, method = proj_method)
      }
      r <- if (isTRUE(terra_todisk)) {
        terra::cover(r, bg, filename = mktemp("cover"), overwrite = TRUE)
      } else {
        terra::cover(r, bg)
      }
    } else {
      r <- terra::ifel(is.na(r), background_value, r,
                       filename = if (isTRUE(terra_todisk)) mktemp("coverconst") else "", overwrite = TRUE)
    }
  }

  if (force_gc) gc()

  # final mask to template
  r <- if (isTRUE(terra_todisk)) {
    terra::mask(r, tmpl, filename = mktemp("maskfinal"), overwrite = TRUE)
  } else {
    terra::mask(r, tmpl)
  }

  # NA count after processing
  n_na_final <- NA_integer_
  if (isTRUE(check_na)) {
    gaps_post <- is.na(r) & !is.na(tmpl)
    n_na_final <- terra::global(gaps_post, fun = "sum", na.rm = TRUE)[[1]]
    if (!quiet) cat("Missing cells after all processing:", n_na_final, "\n")
  }

  # plotting
  if (isTRUE(plot_result) || isTRUE(plot_gaps)) {
    oldpar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(oldpar), add = TRUE)
    if (isTRUE(plot_result) && isTRUE(plot_gaps)) graphics::par(mfrow = c(1, 2))
    if (isTRUE(plot_result)) terra::plot(r, main = "polygon2input: result")
    if (isTRUE(plot_gaps))  {
      gaps_plot <- terra::ifel(is.na(r) & !is.na(tmpl), 1, NA,
                               filename = if (isTRUE(terra_todisk)) mktemp("gapplot") else "", overwrite = TRUE)
      terra::plot(gaps_plot, main = "NA gaps", col = c(NA, "red"), legend = FALSE)
    }
  }

  # atomic write (+ auto dtype/NAflag)
  tuned_defaults <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")
  gdal_opts <- unique(c(gdal_opts, tuned_defaults))
  if (is.null(write_datatype)) write_datatype <- if (value_type == "categorical") "INT2S" else "FLT4S"
  if (is.null(NAflag) && write_datatype == "INT2S") NAflag <- -32768

  # enforce WKT/EPSG before writing
  terra::crs(r) <- tmpl_crs

  ext <- fs::path_ext(out_file); if (!nzchar(ext)) ext <- "tif"
  stem <- fs::path_ext_remove(out_file)
  tmp_out <- sprintf("%s._tmp.%s", stem, ext)
  if (fs::file_exists(tmp_out)) fs::file_delete(tmp_out)

  write_args <- list(x = r, filename = tmp_out, overwrite = TRUE, gdal = gdal_opts)
  if (!is.null(write_datatype)) write_args$datatype <- write_datatype
  if (!is.null(NAflag))        write_args$NAflag   <- NAflag
  do.call(terra::writeRaster, write_args)

  if (fs::file_exists(out_file)) fs::file_delete(out_file)
  fs::file_move(tmp_out, out_file)
  if (!quiet) cat("Wrote:", out_file, "\n")

  invisible(list(
    out_file = out_file,
    n_cells = terra::ncell(tmpl),
    n_na_initial = n_na_initial,
    n_na_final = n_na_final,
    elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    crs = tmpl_crs
  ))
}
