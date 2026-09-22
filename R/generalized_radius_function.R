#' Extract and Rasterize Buffered-Radius Summary Statistics
#'
#' @description
#' Extracts one summary statistic per buffered polygon from one or more raster
#' layers, transfers those values to a rasterization grid using a user-defined
#' key field, rasterizes the results to a common template, and merges the
#' per-tile outputs into final rasters.
#'
#' The relationship between buffered point grids, radii, join keys, and
#' rasterization grids is controlled by `buffer_mode` and, for fully custom
#' workflows, `buffer_spec`. This avoids hard-coding particular grid names or
#' key fields in the extraction logic.
#'
#' @param kvadrati_path Character. Directory containing tiled GeoParquet files
#'   for the base rasterization grid. The tile identifier is taken from the
#'   final underscore-separated component of each file name, before the file
#'   extension. For example, `grid_2434.parquet` has tile ID `2434`.
#' @param radii_path Character. Directory containing tiled buffered-polygon
#'   GeoParquet files. Fixed-radius files are expected to follow the naming
#'   convention `<buffer_grid>_r<radius>_<tileid>.parquet`, matching the output
#'   convention of [tiled_buffers()].
#' @param reference_grid_path Character or `NULL`. Optional path to a
#'   GeoParquet rasterization grid used by the built-in `"sparse"` mode for
#'   coarser buffer grids. The grid must contain the relevant `key_field`
#'   values. It is not required for `"dense"` mode or when a custom
#'   `buffer_spec` supplies all required grid paths.
#' @param template_path Character. Path to the template raster defining final
#'   extent, resolution, alignment, mask, and CRS.
#' @param input_layers Named or unnamed character vector of raster paths from
#'   which summary statistics are extracted. Names are replaced by
#'   `layer_prefixes`.
#' @param layer_prefixes Character vector with one output prefix for each
#'   element of `input_layers`. Values must be unique.
#' @param output_dir Character. Output directory root. Default
#'   `"./Extracted_Layers"`.
#' @param unlink_tiles Logical. If `TRUE`, delete intermediate per-tile raster
#'   directories after successful mosaicking. Default `TRUE`.
#' @param n_workers Integer. Number of parallel workers. On SLURM systems this
#'   is capped at the detected CPU allocation. Default `1`.
#' @param radii Character or numeric vector selecting radii to process. Numeric
#'   values are interpreted as distances and converted to IDs such as `500` to
#'   `"r500"`. Character values may be supplied either as `"500"` or
#'   `"r500"`. In `"sparse"` and `"dense"` modes the default processes
#'   `r500`, `r1250`, `r3000`, and `r10000`. In `"specified"` mode, if
#'   `radii` is omitted, all rows of `buffer_spec` are processed.
#' @param fill_missing Logical. If `TRUE`, fill NA gaps within the template
#'   footprint using `whitebox::wbt_fill_missing_data()`. Default `TRUE`.
#' @param radius_mode Deprecated compatibility argument corresponding to the
#'   former `radius_mode`. Prefer `buffer_mode`. If `buffer_mode` is `NULL`,
#'   this value is used. Default `"sparse"`.
#' @param IDW_weight Numeric. Inverse-distance weighting power passed to
#'   `whitebox::wbt_fill_missing_data()`. Default `2`.
#' @param extract_fun A function or a single character string such as
#'   `"mean"`, passed to `exactextractr::exact_extract()`. It must resolve to
#'   exactly one numeric value per polygon and raster layer.
#' @param future_max_size Numeric. Maximum size, in bytes, of globals exported
#'   to each future worker. Default `8 * 1024^3`.
#' @param gdal_opts Character vector of GDAL creation options used for raster
#'   writes. By default LZW compression, tiling, BigTIFF-if-needed, and 256-cell
#'   blocks are used. Threading is reduced to one GDAL thread per worker when
#'   parallel processing is active.
#' @param write_datatype Character or `NULL`. Optional terra datatype such as
#'   `"FLT4S"` or `"INT2S"`. Default `NULL` uses terra's default.
#' @param NAflag Numeric or `NULL`. Optional NA flag used when writing rasters.
#'   Default `NULL`.
#' @param terra_memfrac Numeric. Fraction passed to
#'   `terra::terraOptions(memfrac = ...)`. Default `0.7`.
#' @param terra_tempdir Character. Temporary directory used by terra and for
#'   intermediate mosaics. Default `tempdir()`.
#' @param terra_todisk Logical. If `TRUE`, prefer on-disk terra operations.
#'   Default `TRUE`.
#' @param quiet Logical. If `TRUE`, suppress progress messages. Default
#'   `FALSE`.
#' @param buffer_mode Character or `NULL`. Preferred replacement for
#'   `radius_mode`. One of `"sparse"`, `"dense"`, or `"specified"`.
#'
#'   * `"sparse"` reproduces the package's standard sparse setup:
#'     `pts100` at 500 and 1250 m, `pts300` at 3000 m, and `pts1000` at
#'     10000 m. The 500 and 1250 m results join directly to each tile by `id`;
#'     the 3000 and 10000 m results use `reference_grid_path` with
#'     `rinda300` and `ID1km`, respectively.
#'   * `"dense"` uses `pts100` for all four standard radii and joins every
#'     result directly to the current tile by `id`.
#'   * `"specified"` uses `buffer_spec` and contains no built-in assumptions
#'     about buffer-grid names, radii, or key fields.
#'
#'   If `NULL`, the value of the compatibility argument `radius_mode` is used.
#' @param buffer_spec A data.frame used when `buffer_mode = "specified"`.
#'   It must contain exactly one processing rule per radius and the columns:
#'
#'   * `buffer_grid`: character prefix identifying buffered files. For example,
#'     `"pts250"` matches files such as `pts250_r2000_2434.parquet`.
#'   * `radius_m`: positive whole-number buffer radius used to construct the
#'     radius ID `r<radius_m>` and locate buffered files.
#'   * `key_field`: character name of the field that links each buffered
#'     polygon to the grid on which its extracted value is rasterized.
#'   * `grid_path`: character path to a GeoParquet rasterization grid, or
#'     `NA`/`"tile"` to use the current tile from `kvadrati_path`.
#'
#'   Each radius must occur only once because each layer-radius combination
#'   produces one output raster.
#' @param tikls100_path Deprecated compatibility alias for
#'   `reference_grid_path`. If supplied, `reference_grid_path` must be `NULL`.
#' @param crop_buffer Numeric. Additional map-unit buffer around the bounding
#'   box of the largest available buffered polygons when cropping input rasters
#'   and the template. Default `1000`.
#'
#' @details
#' ## Buffer specification
#'
#' Internally all three modes are converted to the same four-column
#' specification: `buffer_grid`, `radius_m`, `key_field`, and `grid_path`.
#' Consequently, the extraction and rasterization code does not contain
#' radius-specific rules.
#'
#' If `grid_path` is `NA`, empty, or `"tile"`, extracted values are joined to
#' the current tiled grid from `kvadrati_path`. Otherwise the referenced
#' GeoParquet file is loaded once per worker, filtered to the key values present
#' in the current buffer polygons, and used as the rasterization grid.
#'
#' The `key_field` must exist in both the buffered polygons and the selected
#' rasterization grid. It must uniquely identify buffered polygons within a
#' buffer file. It may be repeated in the rasterization grid; this is what
#' permits a statistic calculated for a coarser buffer center to be assigned to
#' multiple finer output cells.
#'
#' ## Workflow
#'
#' 1. Resolve `buffer_mode` to a validated buffer specification.
#' 2. Discover base-grid tiles and buffered files using their file-name
#'    conventions.
#' 3. For each tile, read the available buffered polygons and identify the
#'    largest available radius.
#' 4. Crop the template and raster stack to the largest buffered extent plus
#'    `crop_buffer`.
#' 5. Run `exactextractr::exact_extract()` once per radius over the cropped
#'    raster stack.
#' 6. Join each extracted statistic to the rule's rasterization grid using
#'    `key_field`, rasterize with `fasterize::fasterize()`, and write per-tile
#'    GeoTIFFs.
#' 7. VRT-mosaic the per-tile rasters, project and mask them to `template_path`,
#'    optionally fill internal NA gaps, and write the final GeoTIFF.
#'
#' ## Parallel and memory behavior
#'
#' Multisession futures are used rather than forking. Large reusable objects,
#' including the template, input raster stack, and external rasterization grids,
#' are cached once per worker through the internal egvtools worker cache. Cache
#' keys include the normalized source paths (and raster layer names where
#' relevant), preventing consecutive calls with different inputs from reusing
#' stale file-backed raster objects. BLAS, OpenMP, and GDAL thread counts are
#' constrained inside workers to reduce oversubscription on HPC systems.
#'
#' ## Gap filling
#'
#' After mosaicking and alignment, gaps are defined as NA cells inside the
#' non-NA template footprint. When `fill_missing = TRUE`, the maximum distance
#' into the gap mask is used to derive an odd Whitebox filter width with a
#' minimum of three cells. Gap counts before and after filling are returned.
#'
#' @return Invisibly returns a data.frame with one row per successfully written
#'   layer-radius combination and columns:
#'
#'   * `layer`: output layer prefix;
#'   * `radius`: radius ID such as `"r3000"`;
#'   * `radius_m`: numeric radius;
#'   * `buffer_grid`: buffered-file grid prefix;
#'   * `output_path`: final raster path;
#'   * `n_tiles_merged`: number of tile rasters included in the mosaic;
#'   * `gaps_before`: number of NA cells inside the template before filling;
#'   * `gaps_after`: number remaining after optional filling;
#'   * `filter_size_used`: Whitebox filter width, or `NA` if not used; and
#'   * `gap_filled`: logical indicating whether gap filling was successfully
#'     applied.
#'
#' @seealso [tiled_buffers()]
#'
#' @examples
#' \dontrun{
#' # Standard sparse workflow.
#' generalized_radius_function(
#'   kvadrati_path = "./Templates/TemplateGrids/lapas/",
#'   radii_path = "./Templates/TemplateGridPoints/lapas/",
#'   reference_grid_path =
#'     "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
#'   template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   input_layers = c(
#'     Soils_txtSand = "./Rastri_100m/RAW/Soils_txtSand_cell.tif",
#'     Soils_txtSilt = "./Rastri_100m/RAW/Soils_txtSilt_cell.tif"
#'   ),
#'   layer_prefixes = c("Soils_txtSand", "Soils_txtSilt"),
#'   buffer_mode = "sparse",
#'   n_workers = 4
#' )
#'
#' # Fully custom relationship between buffered files and rasterization grids.
#' spec <- data.frame(
#'   buffer_grid = c("pts100", "pts100", "pts300", "pts1000"),
#'   radius_m = c(500, 1250, 3000, 10000),
#'   key_field = c("id", "id", "rinda300", "ID1km"),
#'   grid_path = c(
#'     NA,
#'     NA,
#'     "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
#'     "./Templates/TemplateGrids/tikls100_sauzeme.parquet"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#'
#' generalized_radius_function(
#'   kvadrati_path = "./Templates/TemplateGrids/lapas/",
#'   radii_path = "./Templates/TemplateGridPoints/lapas/",
#'   template_path = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   input_layers = c(NDVI = "./Rastri_100m/NDVI_cell.tif"),
#'   layer_prefixes = "NDVI",
#'   buffer_mode = "specified",
#'   buffer_spec = spec,
#'   n_workers = 4
#' )
#'
#' # A custom grid can join directly back to the current tile.
#' custom_spec <- data.frame(
#'   buffer_grid = "pts200",
#'   radius_m = 750,
#'   key_field = "id",
#'   grid_path = NA_character_,
#'   stringsAsFactors = FALSE
#' )
#' }
#'
#' @importFrom terra rast crop mask project writeRaster vrt global distance ifel crs res as.int terraOptions
#' @importFrom sf st_drop_geometry st_as_sfc st_bbox st_buffer st_geometry st_is_empty
#' @importFrom dplyr left_join filter
#' @importFrom furrr future_map2 furrr_options
#' @importFrom future plan sequential multisession
#' @importFrom sfarrow st_read_parquet
#' @importFrom exactextractr exact_extract
#' @importFrom fasterize fasterize
#' @importFrom whitebox wbt_fill_missing_data
#' @importFrom raster raster
#' @importFrom fs dir_exists dir_create is_file
#' @importFrom glue glue
#' @export
generalized_radius_function <- function(
    kvadrati_path,
    radii_path,
    reference_grid_path = NULL,
    template_path,
    input_layers,
    layer_prefixes,
    output_dir = "./Extracted_Layers",
    unlink_tiles = TRUE,
    n_workers = 1,
    radii = c("r500", "r1250", "r3000", "r10000"),
    fill_missing = TRUE,
    radius_mode = "sparse",
    IDW_weight = 2,
    extract_fun = "mean",
    future_max_size = 8 * 1024^3,
    gdal_opts = c(
      "COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER",
      "NUM_THREADS=ALL_CPUS", "BLOCKXSIZE=256", "BLOCKYSIZE=256"
    ),
    write_datatype = NULL,
    NAflag = NULL,
    terra_memfrac = 0.7,
    terra_tempdir = tempdir(),
    terra_todisk = TRUE,
    quiet = FALSE,
    buffer_mode = NULL,
    buffer_spec = NULL,
    tikls100_path = NULL,
    crop_buffer = 1000
) {
  radii_were_supplied <- !missing(radii)
  radius_mode_was_supplied <- !missing(radius_mode)

  # ---- sink safety ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)

  say <- function(...) if (!quiet) cat(..., "\n")

  # ---- dependency checks ----
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' is required for %s. Please install it.", p, why),
        call. = FALSE
      )
    }
  }

  .need_pkg("terra", "raster IO, cropping, mosaicking, and masking")
  .need_pkg("sf", "vector geometry handling")
  .need_pkg("raster", "RasterLayer conversion for fasterize")
  .need_pkg("dplyr", "attribute joins")
  .need_pkg("fs", "filesystem operations")
  .need_pkg("sfarrow", "GeoParquet input")
  .need_pkg("fasterize", "polygon rasterization")
  .need_pkg("exactextractr", "zonal raster extraction")
  .need_pkg("future", "parallel execution")
  .need_pkg("furrr", "parallel mapping")
  .need_pkg("glue", "output paths and messages")
  if (isTRUE(fill_missing)) {
    .need_pkg("whitebox", "gap filling")
  }

  # ---- compatibility aliases ----
  if (!is.null(tikls100_path)) {
    if (!is.null(reference_grid_path)) {
      stop(
        "Supply only one of `reference_grid_path` and deprecated `tikls100_path`.",
        call. = FALSE
      )
    }
    reference_grid_path <- tikls100_path
    warning(
      "`tikls100_path` is deprecated; use `reference_grid_path` instead.",
      call. = FALSE
    )
  }

  if (is.null(buffer_mode)) {
    buffer_mode <- radius_mode
    if (radius_mode_was_supplied) {
      warning(
        "`radius_mode` is deprecated; use `buffer_mode` instead.",
        call. = FALSE
      )
    }
  } else if (radius_mode_was_supplied && !identical(buffer_mode, radius_mode)) {
    stop(
      "`buffer_mode` and deprecated `radius_mode` specify different modes. Use only `buffer_mode`.",
      call. = FALSE
    )
  }

  buffer_mode <- match.arg(buffer_mode, c("sparse", "dense", "specified"))

  # ---- basic argument checks ----
  if (!is.character(kvadrati_path) || length(kvadrati_path) != 1L ||
      !fs::dir_exists(kvadrati_path)) {
    stop("`kvadrati_path` must be an existing directory.", call. = FALSE)
  }
  if (!is.character(radii_path) || length(radii_path) != 1L ||
      !fs::dir_exists(radii_path)) {
    stop("`radii_path` must be an existing directory.", call. = FALSE)
  }
  if (!is.character(template_path) || length(template_path) != 1L ||
      !fs::is_file(template_path)) {
    stop("`template_path` must be an existing file.", call. = FALSE)
  }
  if (!is.character(input_layers) || !length(input_layers)) {
    stop("`input_layers` must be a non-empty character vector of raster paths.", call. = FALSE)
  }
  if (any(!file.exists(input_layers))) {
    missing_layers <- input_layers[!file.exists(input_layers)]
    stop(
      "Input raster file(s) not found: ",
      paste(missing_layers, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.character(layer_prefixes) ||
      length(layer_prefixes) != length(input_layers) ||
      anyNA(layer_prefixes) || any(!nzchar(layer_prefixes))) {
    stop(
      "`layer_prefixes` must be a non-empty character vector with one value per `input_layers` element.",
      call. = FALSE
    )
  }
  if (anyDuplicated(layer_prefixes)) {
    stop("`layer_prefixes` must be unique.", call. = FALSE)
  }
  if (length(extract_fun) != 1L) {
    stop(
      "`extract_fun` must be a single function or string such as \"mean\".",
      call. = FALSE
    )
  }
  if (length(n_workers) != 1L || !is.numeric(n_workers) ||
      !is.finite(n_workers) || n_workers < 1 || n_workers != floor(n_workers)) {
    stop("`n_workers` must be a single whole number >= 1.", call. = FALSE)
  }
  n_workers <- as.integer(n_workers)

  if (length(crop_buffer) != 1L || !is.numeric(crop_buffer) ||
      !is.finite(crop_buffer) || crop_buffer < 0) {
    stop("`crop_buffer` must be a single finite number >= 0.", call. = FALSE)
  }
  if (length(IDW_weight) != 1L || !is.numeric(IDW_weight) ||
      !is.finite(IDW_weight) || IDW_weight <= 0) {
    stop("`IDW_weight` must be a single positive finite number.", call. = FALSE)
  }
  if (length(future_max_size) != 1L || !is.numeric(future_max_size) ||
      !is.finite(future_max_size) || future_max_size <= 0) {
    stop("`future_max_size` must be a single positive finite number.", call. = FALSE)
  }
  if (length(terra_memfrac) != 1L || !is.numeric(terra_memfrac) ||
      !is.finite(terra_memfrac) || terra_memfrac <= 0 || terra_memfrac > 1) {
    stop("`terra_memfrac` must be in (0, 1].", call. = FALSE)
  }

  names(input_layers) <- layer_prefixes

  if (!fs::dir_exists(output_dir)) {
    fs::dir_create(output_dir, recurse = TRUE)
  }

  # ---- normalize paths used repeatedly ----
  kvadrati_path <- normalizePath(kvadrati_path, winslash = "/", mustWork = TRUE)
  radii_path <- normalizePath(radii_path, winslash = "/", mustWork = TRUE)
  template_path <- normalizePath(template_path, winslash = "/", mustWork = TRUE)
  input_layers <- vapply(
    input_layers,
    normalizePath,
    character(1),
    winslash = "/",
    mustWork = TRUE
  )
  names(input_layers) <- layer_prefixes
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (!dir.exists(terra_tempdir)) {
    dir.create(terra_tempdir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(terra_tempdir)) {
    stop("Could not create `terra_tempdir`: ", terra_tempdir, call. = FALSE)
  }
  terra_tempdir <- normalizePath(terra_tempdir, winslash = "/", mustWork = TRUE)

  if (!is.null(reference_grid_path)) {
    if (!is.character(reference_grid_path) || length(reference_grid_path) != 1L ||
        !fs::is_file(reference_grid_path)) {
      stop("`reference_grid_path` must be an existing file or NULL.", call. = FALSE)
    }
    reference_grid_path <- normalizePath(
      reference_grid_path,
      winslash = "/",
      mustWork = TRUE
    )
  }

  # ---- common radius helpers ----
  .radius_id <- function(x) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0)) {
      stop("`radius_m` values must be positive finite numbers.", call. = FALSE)
    }
    tol <- sqrt(.Machine$double.eps)
    if (any(abs(x - round(x)) > tol)) {
      stop(
        "`radius_m` values must be whole numbers because tiled buffer file names use integer radii.",
        call. = FALSE
      )
    }
    paste0("r", as.integer(round(x)))
  }

  .normalize_radii_filter <- function(x) {
    if (is.numeric(x)) {
      return(.radius_id(x))
    }
    if (!is.character(x) || anyNA(x) || any(!nzchar(x))) {
      stop("`radii` must be numeric or a non-empty character vector.", call. = FALSE)
    }
    out <- ifelse(startsWith(tolower(x), "r"), x, paste0("r", x))
    paste0("r", sub("^[rR]", "", out))
  }

  # ---- resolve every mode to one common specification ----
  if (identical(buffer_mode, "sparse")) {
    if (is.null(reference_grid_path)) {
      stop(
        "`reference_grid_path` is required for `buffer_mode = \"sparse\"`. ",
        "Alternatively use `buffer_mode = \"specified\"` and provide `buffer_spec`.",
        call. = FALSE
      )
    }
    spec <- data.frame(
      buffer_grid = c("pts100", "pts100", "pts300", "pts1000"),
      radius_m = c(500, 1250, 3000, 10000),
      key_field = c("id", "id", "rinda300", "ID1km"),
      grid_path = c(NA_character_, NA_character_,
                    reference_grid_path, reference_grid_path),
      stringsAsFactors = FALSE
    )
  } else if (identical(buffer_mode, "dense")) {
    spec <- data.frame(
      buffer_grid = rep("pts100", 4L),
      radius_m = c(500, 1250, 3000, 10000),
      key_field = rep("id", 4L),
      grid_path = rep(NA_character_, 4L),
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(buffer_spec) || !is.data.frame(buffer_spec)) {
      stop(
        "`buffer_spec` must be a data.frame when `buffer_mode = \"specified\"`.",
        call. = FALSE
      )
    }
    required_cols <- c("buffer_grid", "radius_m", "key_field", "grid_path")
    missing_cols <- setdiff(required_cols, names(buffer_spec))
    if (length(missing_cols)) {
      stop(
        "`buffer_spec` is missing required column(s): ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    spec <- buffer_spec[, required_cols, drop = FALSE]
  }

  if (!nrow(spec)) {
    stop("The resolved buffer specification contains no rows.", call. = FALSE)
  }

  spec$buffer_grid <- as.character(spec$buffer_grid)
  spec$key_field <- as.character(spec$key_field)
  spec$grid_path <- as.character(spec$grid_path)
  spec$radius_m <- suppressWarnings(as.numeric(spec$radius_m))

  if (anyNA(spec$buffer_grid) || any(!nzchar(spec$buffer_grid))) {
    stop("`buffer_spec$buffer_grid` must contain non-empty values.", call. = FALSE)
  }
  if (anyNA(spec$key_field) || any(!nzchar(spec$key_field))) {
    stop("`buffer_spec$key_field` must contain non-empty field names.", call. = FALSE)
  }

  spec$radius_id <- .radius_id(spec$radius_m)

  if (anyDuplicated(spec$radius_id)) {
    dup <- unique(spec$radius_id[duplicated(spec$radius_id)])
    stop(
      "Each radius may occur only once in `buffer_spec`. Duplicated radius ID(s): ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }

  # NA, empty, and literal "tile" all mean the per-tile base grid.
  use_tile_grid <- is.na(spec$grid_path) |
    !nzchar(spec$grid_path) |
    tolower(spec$grid_path) == "tile"
  spec$grid_path[use_tile_grid] <- NA_character_

  external_idx <- which(!is.na(spec$grid_path))
  if (length(external_idx)) {
    for (i in external_idx) {
      p <- spec$grid_path[i]
      if (!fs::is_file(p)) {
        stop(
          "External rasterization grid not found for ", spec$radius_id[i],
          ": ", p,
          call. = FALSE
        )
      }
      spec$grid_path[i] <- normalizePath(p, winslash = "/", mustWork = TRUE)
    }
  }

  # Preserve legacy default radii for sparse/dense, but in specified mode an
  # omitted `radii` means all radii in the specification.
  if (identical(buffer_mode, "specified") && !radii_were_supplied) {
    requested_radius_ids <- spec$radius_id
  } else {
    requested_radius_ids <- .normalize_radii_filter(radii)
  }

  requested_radius_ids <- unique(requested_radius_ids)

  missing_requested <- setdiff(requested_radius_ids, spec$radius_id)
  if (length(missing_requested)) {
    stop(
      "Requested radius/radii not defined by the resolved buffer specification: ",
      paste(missing_requested, collapse = ", "),
      call. = FALSE
    )
  }

  spec <- spec[spec$radius_id %in% requested_radius_ids, , drop = FALSE]
  spec <- spec[match(requested_radius_ids, spec$radius_id), , drop = FALSE]
  rownames(spec) <- NULL

  if (!nrow(spec)) {
    stop("No radii remain after applying `radii`.", call. = FALSE)
  }

  # ---- terra options: set and restore ----
  old_opt <- NULL
  utils::capture.output({
    old_opt <- terra::terraOptions()
  })
  on.exit(
    terra::terraOptions(
      memfrac = old_opt$memfrac,
      tempdir = old_opt$tempdir,
      todisk = old_opt$todisk,
      progress = old_opt$progress
    ),
    add = TRUE
  )
  utils::capture.output({
    terra::terraOptions(
      memfrac = terra_memfrac,
      tempdir = terra_tempdir,
      progress = FALSE
    )
    if (!is.na(terra_todisk)) {
      terra::terraOptions(todisk = isTRUE(terra_todisk))
    }
  })

  # ---- parallel plan ----
  options(future.fork.enable = FALSE)

  slurm_vals <- suppressWarnings(as.integer(c(
    Sys.getenv("SLURM_CPUS_PER_TASK"),
    Sys.getenv("SLURM_CPUS_ON_NODE")
  )))
  slurm_vals <- slurm_vals[is.finite(slurm_vals) & slurm_vals > 0L]
  if (length(slurm_vals)) {
    n_workers <- min(n_workers, max(slurm_vals))
  }

  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = future_max_size)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)

  is_hpc <- function() {
    any(nzchar(Sys.getenv(c(
      "SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID",
      "APPTAINER_NAME", "SINGULARITY_NAME"
    ))))
  }

  show_progress <- !quiet && !is_hpc()
  scratch <- Sys.getenv("SLURM_TMPDIR", unset = tempdir())
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    BLIS_NUM_THREADS = "1",
    GOTO_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1",
    GDAL_NUM_THREADS = "1",
    MALLOC_ARENA_MAX = "2",
    GDAL_CACHEMAX = "256",
    CPL_VSIL_CURL_CACHE_SIZE = "0",
    TMPDIR = scratch,
    TEMP = scratch,
    R_TEMPORARY_DIR = scratch,
    CPL_TMPDIR = scratch
  )

  if (n_workers <= 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }

  # ---- discover tiled base grids ----
  tile_files <- list.files(
    kvadrati_path,
    pattern = "\\.parquet$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(tile_files)) {
    stop("No GeoParquet tile files found in `kvadrati_path`.", call. = FALSE)
  }

  tile_stems <- tools::file_path_sans_ext(basename(tile_files))
  tile_ids <- sub("^.*_", "", tile_stems)
  if (any(!nzchar(tile_ids))) {
    stop("Could not derive a tile ID from one or more files in `kvadrati_path`.", call. = FALSE)
  }
  if (anyDuplicated(tile_ids)) {
    dup <- unique(tile_ids[duplicated(tile_ids)])
    stop(
      "Multiple base-grid files resolve to the same tile ID: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }

  tile_table <- data.frame(
    tile_id = tile_ids,
    tile_path = tile_files,
    stringsAsFactors = FALSE
  )

  # ---- discover buffered files according to the resolved spec ----
  buffer_files <- list.files(
    radii_path,
    pattern = "\\.parquet$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(buffer_files)) {
    stop("No GeoParquet buffer files found in `radii_path`.", call. = FALSE)
  }

  buffer_stems <- tools::file_path_sans_ext(basename(buffer_files))
  index_parts <- vector("list", nrow(spec))

  for (i in seq_len(nrow(spec))) {
    prefix <- paste0(spec$buffer_grid[i], "_", spec$radius_id[i], "_")
    hits <- startsWith(tolower(buffer_stems), tolower(prefix))
    if (!any(hits)) next

    matched_stems <- buffer_stems[hits]
    matched_files <- buffer_files[hits]
    matched_tile_ids <- substring(matched_stems, nchar(prefix) + 1L)

    if (any(!nzchar(matched_tile_ids))) {
      stop(
        "Could not derive tile IDs from buffered files for ",
        spec$buffer_grid[i], " / ", spec$radius_id[i], ".",
        call. = FALSE
      )
    }

    index_parts[[i]] <- data.frame(
      tile_id = matched_tile_ids,
      radius_id = spec$radius_id[i],
      buffer_grid = spec$buffer_grid[i],
      buffer_path = matched_files,
      stringsAsFactors = FALSE
    )
  }

  buffer_index <- do.call(rbind, index_parts[!vapply(index_parts, is.null, logical(1))])
  if (is.null(buffer_index) || !nrow(buffer_index)) {
    stop(
      "No buffered files in `radii_path` matched the resolved buffer specification.",
      call. = FALSE
    )
  }

  index_key <- paste(buffer_index$tile_id, buffer_index$radius_id, sep = "\r")
  if (anyDuplicated(index_key)) {
    dup <- unique(index_key[duplicated(index_key)])
    stop(
      "More than one buffered file was found for the same tile/radius combination: ",
      paste(gsub("\r", " / ", dup, fixed = TRUE), collapse = ", "),
      call. = FALSE
    )
  }

  # Build one compact job object per base-grid tile.
  tile_jobs <- lapply(seq_len(nrow(tile_table)), function(i) {
    id <- tile_table$tile_id[i]
    list(
      tile_path = tile_table$tile_path[i],
      buffers = buffer_index[buffer_index$tile_id == id, , drop = FALSE]
    )
  })
  names(tile_jobs) <- tile_table$tile_id

  # ---- template ----
  template_r_full <- terra::rast(template_path)

  # ---- GDAL options ----
  tuned_defaults <- c(
    "COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER",
    "BLOCKXSIZE=256", "BLOCKYSIZE=256",
    if (n_workers > 1L) "NUM_THREADS=1" else "NUM_THREADS=ALL_CPUS"
  )
  gdal_opts <- gdal_opts[!grepl("^NUM_THREADS=", gdal_opts)]
  gdal_opts <- unique(c(gdal_opts, tuned_defaults))

  # ---- atomic raster write ----
  .write_r <- function(r, path, gdal_opts, write_datatype, NAflag) {
    ext <- tools::file_ext(path)
    if (!nzchar(ext)) ext <- "tif"
    stem <- sub(sprintf("\\.%s$", ext), "", path)
    tmp <- sprintf("%s._tmp.%s", stem, ext)

    if (file.exists(tmp)) {
      try(unlink(tmp), silent = TRUE)
    }

    args <- list(
      x = r,
      filename = tmp,
      overwrite = TRUE,
      gdal = gdal_opts
    )
    if (!is.null(write_datatype)) args$datatype <- write_datatype
    if (!is.null(NAflag)) args$NAflag <- NAflag

    do.call(terra::writeRaster, args)

    if (file.exists(path)) {
      try(unlink(path), silent = TRUE)
    }
    if (!file.rename(tmp, path)) {
      stop("Failed to atomically move temporary raster to: ", path, call. = FALSE)
    }
    invisible(path)
  }

  # ---- per-worker cache wrapper ----
  .get_cached <- function(name, loader) {
    x <- .egv_cache_get(name)
    if (!is.null(x)) return(x)
    x <- loader()
    .egv_cache_set(name, x)
    x
  }

  .get_reference_grid <- function(path) {
    cache_key <- paste0(".radius_join_grid::", path)
    .get_cached(cache_key, function() sfarrow::st_read_parquet(path))
  }

  # ---- parser for exact_extract output ----
  .parse_extract <- function(res, layer_names, fun_name_or_fn, n_features) {
    if (is.matrix(res)) {
      res <- as.data.frame(res)
    }

    # Single numeric vector: valid only for a single raster layer.
    if (is.numeric(res) && is.null(dim(res))) {
      if (length(layer_names) != 1L) {
        stop(
          "exact_extract() returned one numeric vector for multiple raster layers.",
          call. = FALSE
        )
      }
      if (length(res) != n_features) {
        stop(
          "exact_extract() returned ", length(res),
          " values for ", n_features, " polygons.",
          call. = FALSE
        )
      }
      return(stats::setNames(list(as.numeric(res)), layer_names[1]))
    }

    # Data frame: normally one column per raster layer.
    if (is.data.frame(res)) {
      if (nrow(res) != n_features) {
        stop(
          "exact_extract() returned a data.frame with ", nrow(res),
          " rows for ", n_features, " polygons.",
          call. = FALSE
        )
      }

      cn <- colnames(res)
      fun_str <- if (is.character(fun_name_or_fn)) {
        tolower(fun_name_or_fn)
      } else {
        NA_character_
      }

      stripped <- if (!is.na(fun_str)) {
        sub(paste0("^", fun_str, "\\."), "", cn)
      } else {
        cn
      }

      if (!all(stripped %in% layer_names)) {
        stripped <- sub("^.*\\.", "", cn)
      }

      final_names <- if (all(stripped %in% layer_names)) stripped else cn
      keep_idx <- final_names %in% layer_names

      if (!any(keep_idx)) {
        stop(
          "Could not match exact_extract() output columns to raster layer names.",
          call. = FALSE
        )
      }

      res2 <- res[, keep_idx, drop = FALSE]
      final_names <- final_names[keep_idx]

      if (anyDuplicated(final_names)) {
        stop(
          "exact_extract() produced duplicated layer names after parsing.",
          call. = FALSE
        )
      }

      out <- stats::setNames(vector("list", length(final_names)), final_names)
      for (i in seq_along(final_names)) {
        out[[i]] <- as.numeric(res2[[i]])
      }
      return(out)
    }

    # List: one element per feature, each element representing the layers.
    if (is.list(res)) {
      if (length(res) != n_features) {
        stop(
          "exact_extract() returned a list of length ", length(res),
          " for ", n_features, " polygons.",
          call. = FALSE
        )
      }

      norm_one <- function(el) {
        if (is.numeric(el)) {
          v <- as.numeric(el)
          if (length(v) == 1L && length(layer_names) == 1L) {
            names(v) <- layer_names[1]
            return(v)
          }
          if (length(v) == length(layer_names)) {
            names(v) <- layer_names
            return(v)
          }
          stop(
            "A list element from exact_extract() has length ", length(v),
            " but ", length(layer_names), " raster layer(s) were supplied.",
            call. = FALSE
          )
        }

        if (is.data.frame(el)) {
          if (nrow(el) != 1L) {
            stop(
              "A data.frame list element from exact_extract() has ", nrow(el),
              " rows; expected one.",
              call. = FALSE
            )
          }
          cn <- colnames(el)
          cn2 <- sub("^.*\\.", "", cn)

          if (all(cn2 %in% layer_names)) {
            v <- as.numeric(el[1, , drop = TRUE])
            names(v) <- cn2
            return(v)
          }
          if (length(cn) == length(layer_names)) {
            v <- as.numeric(el[1, , drop = TRUE])
            names(v) <- layer_names
            return(v)
          }
          stop(
            "Cannot map exact_extract() data.frame columns to raster layer names.",
            call. = FALSE
          )
        }

        stop(
          "Unsupported list element type from exact_extract(): ", class(el)[1],
          call. = FALSE
        )
      }

      mat <- lapply(res, norm_one)
      M <- do.call(
        rbind,
        lapply(mat, function(v) stats::setNames(v[layer_names], layer_names))
      )

      out <- stats::setNames(vector("list", length(layer_names)), layer_names)
      for (i in seq_along(layer_names)) {
        out[[i]] <- as.numeric(M[, i])
      }
      return(out)
    }

    stop(
      "Unsupported result type from exact_extract(): ", class(res)[1],
      call. = FALSE
    )
  }

  # ---- per-tile worker ----
  process_tile <- function(
    solis,
    job,
    spec,
    template_path,
    input_layers,
    layer_prefixes,
    output_dir,
    extract_fun,
    gdal_opts,
    write_datatype,
    NAflag,
    crop_buffer
  ) {
    say("Processing tile: ", solis)

    # Cache keys include the normalized source paths (and raster names) so that
    # consecutive calls to generalized_radius_function() cannot accidentally
    # reuse SpatRaster objects created for different input files.
    template_cache_key <- paste0(
      ".radius_template_cache::",
      template_path
    )
    stack_cache_key <- paste0(
      ".radius_stack_cache::",
      paste0(layer_prefixes, "=", input_layers, collapse = "\r")
    )

    template_100 <- .get_cached(
      template_cache_key,
      function() terra::rast(template_path)
    )
    stack_all <- .get_cached(
      stack_cache_key,
      function() {
        rs <- terra::rast(input_layers)
        names(rs) <- layer_prefixes
        rs
      }
    )

    # Base grid for this tile.
    tile_grid <- sfarrow::st_read_parquet(job$tile_path)
    if (!inherits(tile_grid, "sf") || is.null(sf::st_geometry(tile_grid))) {
      stop("Base-grid tile is not a valid sf object: ", job$tile_path, call. = FALSE)
    }

    # Read each available buffer rule for the tile.
    r_polys <- stats::setNames(vector("list", nrow(spec)), spec$radius_id)

    for (i in seq_len(nrow(spec))) {
      rad <- spec$radius_id[i]
      hit <- job$buffers$buffer_path[job$buffers$radius_id == rad]
      if (!length(hit)) next

      vec <- tryCatch(
        sfarrow::st_read_parquet(hit[1]),
        error = function(e) NULL
      )
      if (!inherits(vec, "sf") || is.null(sf::st_geometry(vec))) next

      vec <- vec[!sf::st_is_empty(vec), , drop = FALSE]
      if (!nrow(vec)) next

      key <- spec$key_field[i]
      if (!key %in% names(vec)) {
        stop(
          "Key field '", key, "' is missing from buffered polygons: ", hit[1],
          call. = FALSE
        )
      }
      if (anyNA(vec[[key]])) {
        stop(
          "Key field '", key, "' contains NA values in buffered polygons: ",
          hit[1],
          call. = FALSE
        )
      }
      if (anyDuplicated(vec[[key]])) {
        stop(
          "Key field '", key, "' does not uniquely identify buffered polygons for ",
          rad, " in tile ", solis, ".",
          call. = FALSE
        )
      }

      r_polys[[rad]] <- vec
    }

    available <- names(r_polys)[!vapply(r_polys, is.null, logical(1))]
    if (!length(available)) {
      stop("No valid buffered polygons found for tile ", solis, ".", call. = FALSE)
    }

    available_rows <- match(available, spec$radius_id)
    largest_rad <- available[which.max(spec$radius_m[available_rows])]
    largest <- r_polys[[largest_rad]]

    telpa2 <- sf::st_buffer(
      sf::st_as_sfc(sf::st_bbox(largest)),
      dist = crop_buffer
    )
    template_crop <- terra::crop(template_100, telpa2)
    templateRL <- raster::raster(template_crop)
    stack_crop <- terra::crop(stack_all, telpa2)

    for (i in seq_len(nrow(spec))) {
      rad <- spec$radius_id[i]
      vec <- r_polys[[rad]]
      if (is.null(vec)) next

      key <- spec$key_field[i]
      grid_path <- spec$grid_path[i]

      res <- suppressWarnings(
        exactextractr::exact_extract(stack_crop, vec, fun = extract_fun)
      )
      vals_by_layer <- .parse_extract(
        res,
        layer_names = names(stack_crop),
        fun_name_or_fn = extract_fun,
        n_features = nrow(vec)
      )

      # Select the grid that receives the extracted buffer statistic.
      if (is.na(grid_path)) {
        joined_base <- tile_grid
      } else {
        ref_grid <- .get_reference_grid(grid_path)
        if (!inherits(ref_grid, "sf") || is.null(sf::st_geometry(ref_grid))) {
          stop("Reference grid is not a valid sf object: ", grid_path, call. = FALSE)
        }
        if (!key %in% names(ref_grid)) {
          stop(
            "Key field '", key, "' is missing from reference grid: ", grid_path,
            call. = FALSE
          )
        }
        ids <- unique(vec[[key]])
        joined_base <- ref_grid[ref_grid[[key]] %in% ids, , drop = FALSE]
      }

      if (!key %in% names(joined_base)) {
        stop(
          "Key field '", key, "' is missing from the rasterization grid for ",
          rad, " in tile ", solis, ".",
          call. = FALSE
        )
      }
      if (!nrow(joined_base)) {
        stop(
          "No rasterization-grid features matched key field '", key,
          "' for ", rad, " in tile ", solis, ".",
          call. = FALSE
        )
      }

      xdf <- sf::st_drop_geometry(vec)[key]

      for (prefix in names(vals_by_layer)) {
        v <- vals_by_layer[[prefix]]
        if (length(v) != nrow(vec)) {
          stop(
            "Row mismatch between extracted values and buffered polygons at ",
            rad, " for layer ", prefix, ".",
            call. = FALSE
          )
        }

        xdf$vertibas <- as.numeric(v)
        joined <- dplyr::left_join(joined_base, xdf, by = key)

        rr <- fasterize::fasterize(joined, templateRL, field = "vertibas")
        rr <- terra::rast(rr)
        rr <- terra::mask(rr, template_crop)

        out_path <- glue::glue(
          "{output_dir}/{prefix}_{rad}/{prefix}_{rad}_{solis}.tif"
        )
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
        .write_r(rr, out_path, gdal_opts, write_datatype, NAflag)
      }
    }

    solis
  }

  # ---- run tiles in parallel ----
  opts <- furrr::furrr_options(
    seed = TRUE,
    globals = c(
      "process_tile",
      ".get_cached",
      ".get_reference_grid",
      ".parse_extract",
      ".write_r",
      "say",
      ".egvtools_cache",
      ".egv_cache_get",
      ".egv_cache_set",
      ".egv_cache_drop",
      ".egv_cache_clear"
    ),
    packages = c(
      "terra", "sf", "sfarrow", "fasterize", "exactextractr",
      "dplyr", "glue", "raster"
    )
  )

  furrr::future_map2(
    .x = names(tile_jobs),
    .y = tile_jobs,
    .f = function(x, y) {
      process_tile(
        solis = x,
        job = y,
        spec = spec,
        template_path = template_path,
        input_layers = input_layers,
        layer_prefixes = layer_prefixes,
        output_dir = output_dir,
        extract_fun = extract_fun,
        gdal_opts = gdal_opts,
        write_datatype = write_datatype,
        NAflag = NAflag,
        crop_buffer = crop_buffer
      )
    },
    .progress = show_progress,
    .options = opts
  )

  # ---- mosaic, align, gap analysis, optional fill, final write ----
  template_r <- template_r_full
  pix_size <- mean(terra::res(template_r))
  results_ls <- list()

  for (prefix in layer_prefixes) {
    for (i in seq_len(nrow(spec))) {
      rad <- spec$radius_id[i]
      tif_dir <- glue::glue("{output_dir}/{prefix}_{rad}")
      out_files <- list.files(
        tif_dir,
        pattern = "\\.tif$",
        full.names = TRUE,
        ignore.case = TRUE
      )
      if (!length(out_files)) next

      vrt <- terra::vrt(out_files)
      names(vrt) <- glue::glue("{prefix}_{rad}")

      temp_raw <- file.path(
        terra_tempdir,
        glue::glue("mosaic_{prefix}_{rad}.tif")
      )
      .write_r(vrt, temp_raw, gdal_opts, write_datatype, NAflag)

      proj_r <- terra::project(terra::rast(temp_raw), template_r)
      proj_r <- terra::mask(proj_r, template_r)

      robi <- is.na(proj_r) & !is.na(template_r)
      gaps_before <- as.integer(
        terra::global(terra::as.int(robi), fun = "sum", na.rm = TRUE)[[1]]
      )
      if (is.na(gaps_before)) gaps_before <- 0L

      gap_filled <- FALSE
      filter_used <- NA_integer_

      if (isTRUE(fill_missing) && gaps_before > 0L) {
        say(glue::glue(
          "[{prefix}_{rad}] Filling {format(gaps_before, big.mark = ',')} NA cells"
        ))

        fillable <- terra::ifel(!robi, 1, NA)
        dist_r <- terra::distance(fillable)
        max_att <- suppressWarnings(
          terra::global(dist_r, fun = "max", na.rm = TRUE)[[1]]
        )
        rm(fillable, dist_r)

        if (is.finite(max_att)) {
          filter_used <- max(
            3L,
            as.integer(ceiling(as.numeric(max_att) / pix_size) * 2L)
          )
          if (filter_used %% 2L == 0L) {
            filter_used <- filter_used + 1L
          }

          temp_filled <- file.path(
            terra_tempdir,
            glue::glue("filled_{prefix}_{rad}.tif")
          )

          ok <- TRUE
          tryCatch(
            {
              whitebox::wbt_fill_missing_data(
                input = temp_raw,
                output = temp_filled,
                filter = filter_used,
                weight = IDW_weight,
                no_edges = FALSE
              )
            },
            error = function(e) {
              say(glue::glue(
                "[{prefix}_{rad}] Whitebox failed: {conditionMessage(e)}"
              ))
              ok <<- FALSE
            }
          )

          if (ok && file.exists(temp_filled)) {
            proj_r <- terra::project(terra::rast(temp_filled), template_r)
            proj_r <- terra::mask(proj_r, template_r)
            unlink(temp_filled)
            gap_filled <- TRUE
          }
        }
      }

      robi2 <- is.na(proj_r) & !is.na(template_r)
      gaps_after <- as.integer(
        terra::global(terra::as.int(robi2), fun = "sum", na.rm = TRUE)[[1]]
      )
      if (is.na(gaps_after)) gaps_after <- 0L

      names(proj_r) <- glue::glue("{prefix}_{rad}")
      terra::crs(proj_r) <- terra::crs(template_r, proj = FALSE)

      out_final <- glue::glue("{output_dir}/{prefix}_{rad}.tif")
      .write_r(proj_r, out_final, gdal_opts, write_datatype, NAflag)

      if (file.exists(temp_raw)) unlink(temp_raw)
      if (unlink_tiles) {
        try(unlink(tif_dir, recursive = TRUE), silent = TRUE)
      }

      results_ls[[length(results_ls) + 1L]] <- data.frame(
        layer = prefix,
        radius = rad,
        radius_m = spec$radius_m[i],
        buffer_grid = spec$buffer_grid[i],
        output_path = out_final,
        n_tiles_merged = length(out_files),
        gaps_before = gaps_before,
        gaps_after = gaps_after,
        filter_size_used = if (isTRUE(fill_missing)) filter_used else NA_integer_,
        gap_filled = isTRUE(gap_filled),
        stringsAsFactors = FALSE
      )
    }
  }

  res_df <- if (length(results_ls)) {
    do.call(rbind, results_ls)
  } else {
    data.frame(
      layer = character(),
      radius = character(),
      radius_m = numeric(),
      buffer_grid = character(),
      output_path = character(),
      n_tiles_merged = integer(),
      gaps_before = integer(),
      gaps_after = integer(),
      filter_size_used = integer(),
      gap_filled = logical(),
      stringsAsFactors = FALSE
    )
  }

  invisible(res_df)
}
