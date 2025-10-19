#' Extract and Rasterize Summary Statistics from Buffered Radii Using exactextractr
#'
#' Extracts summary statistics from raster layers using buffered polygon zones of multiple radii
#' and rasterizes them onto a common template grid. Designed for parallel execution.
#'
#' @details
#' **Workflow**
#' 1. Discover inputs (1 ha tiles and buffered radii files) according to `radius_mode`.
#' 2. **Per tile (parallel):** read tile zones and available radii; crop to bbox + 1 km;
#'    **run `exactextractr::exact_extract()` once per radius** over the (cropped) raster stack,
#'    producing one value per polygon per layer; join to the appropriate 100 m grid (or its
#'    parent for r3000/r10000), and rasterize to the cropped template via `fasterize`,
#'    writing **per-tile GeoTIFFs**.
#' 3. For each `layer at radius`: VRT-merge tiles then **project to the full template grid**
#'    then mask to the template footprint.
#' 4. **Gap analysis & optional fill:** count NA gaps inside the template. If `fill_missing = TRUE`
#'    and gaps exist, estimate max gap width on the template (`terra::distance()` on a fillable
#'    mask), set Whitebox `filter = 2 * ceil(max_gap / pixel_size)` with **min clamp = 3**,
#'    run `whitebox::wbt_fill_missing_data()`, then mask again.
#' 5. **CRS guard & write:** set `terra::crs(out) <- terra::crs(template, proj = FALSE)` and
#'    write with LZW tiling via **atomic writes**.
#' 6. **Return:** a data.frame with
#'    `layer, radius, output_path, n_tiles_merged, gaps_before, gaps_after, filter_size_used, gap_filled`.
#'
#' **Output parsing guard**
#' Handles single vectors, data.frames (`"mean.Layer"` or layer names), and lists
#' (one element per feature: numeric vectors or 1-row data.frames). Clear error otherwise.
#'
#' **Per-worker cache**
#' Large `sf` objects like a ~6.5M-feature `tikls100` are **loaded once per worker**
#' to avoid repeated I/O and pointer invalidation issues in parallel sessions.
#'
#' **CRS guard**
#' Outputs copy the CRS string from the template using
#' `terra::crs(template, proj = FALSE)` to keep a Proj.4-style string for maximal compatibility.
#'
#' @param kvadrati_path Path to directory containing 1 ha grid-cell GeoParquet tiles (e.g., `.../lapas/`).
#' @param radii_path Path to directory with buffered polygon radii tiles (`r500`, `r1250`, `r3000`, `r10000`).
#' @param tikls100_path Path to GeoParquet for the 100 m grid (must contain `rinda300` and `ID1km`).
#' @param template_path Path to the template raster (final alignment target; e.g., `LV100m_10km.tif`).
#' @param input_layers Named character vector of covariate raster paths (values are file paths,
#'   names are the **layer prefixes** used in outputs).
#' @param layer_prefixes Character vector of names (same length/order as `input_layers`) to use as
#'   output layer prefixes.
#' @param output_dir Output directory root. Default `"./Extracted_Layers"`.
#' @param unlink_tiles Logical; delete per-tile rasters after merging. Default `TRUE`.
#' @param n_workers Number of parallel workers. Default `1`.
#' @param radii Character vector of radii to process (subset of: `"r500"`, `"r1250"`, `"r3000"`, `"r10000"`).
#' @param fill_missing If `TRUE`, run Whitebox IDW gap filling on the projected mosaic.
#' @param radius_mode `"sparse"` or `"dense"`. Controls which buffered files are used:
#'   *`"sparse"`* uses `pts100` for r500/r1250, `pts300` for r3000, and `pts1000/pts1km` for r10000;
#'   *`"dense"`* uses only `pts100` for all radii.
#' @param IDW_weight Numeric power for Whitebox IDW. Default `2`.
#' @param extract_fun Function or single string (e.g., `"mean"`) passed to `exactextractr::exact_extract()`
#'   that returns **one scalar per polygon**. If it returns multiple values per polygon, the function stops.
#' @param future_max_size Max globals per worker (bytes). Default `8 * 1024^3`.
#' @param gdal_opts GDAL creation options for writes. Default
#'   `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`.
#' @param write_datatype Optional terra datatype (e.g., `"FLT4S"`, `"INT2S"`). Default `NULL` (terra default).
#' @param NAflag Optional NA flag for writing. Default `NULL` (terra default).
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for heavy ops. Default `tempdir()`.
#' @param terra_todisk If `TRUE`, prefer on-disk operations. Default `TRUE`.
#' @param quiet Suppress progress prints. Default `FALSE`.
#'
#' @return Invisibly, a data.frame with columns:
#'   `layer, radius, output_path, n_tiles_merged, gaps_before, gaps_after, filter_size_used, gap_filled`.
#'
#' @examples
#' \dontrun{
#' radius_function(
#'   kvadrati_path  = "./Templates/TemplateGrids/lapas/",
#'   radii_path     = "./Templates/TemplateGridPoints/lapas/",
#'   tikls100_path  = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
#'   template_path  = "./Templates/TemplateRasters/LV100m_10km.tif",
#'   input_layers   = c(Soils_txtSand  = "./Rastri_100m/RAW/Soils_txtSand_cell.tif",
#'                      Soils_txtSilt  = "./Rastri_100m/RAW/Soils_txtSilt_cell.tif",
#'                      Soils_txtClay  = "./Rastri_100m/RAW/Soils_txtClay_cell.tif"),
#'   layer_prefixes = c("Soils_txtSand","Soils_txtSilt","Soils_txtClay"),
#'   output_dir     = "./Extracted_Layers",
#'   n_workers      = 4,
#'   radii          = c("r500","r1250","r3000","r10000"),
#'   radius_mode    = "sparse",
#'   extract_fun    = "mean",
#'   fill_missing   = TRUE,
#'   IDW_weight     = 2
#' )
#' }
#'
#' @importFrom stats setNames
#' @importFrom rlang .data :=
#' @importFrom terra rast crop mask project writeRaster vrt global distance ifel crs res as.int
#' @importFrom sf st_drop_geometry st_as_sfc st_bbox st_buffer st_geometry
#' @importFrom dplyr mutate select left_join case_when filter rename_with starts_with
#' @importFrom purrr map set_names
#' @importFrom glue glue
#' @importFrom tibble tibble
#' @importFrom tidyr separate pivot_wider
#' @importFrom furrr future_map2 furrr_options
#' @importFrom future plan sequential multisession
#' @importFrom sfarrow st_read_parquet
#' @importFrom exactextractr exact_extract
#' @importFrom whitebox wbt_fill_missing_data
#' @importFrom raster raster
#' @importFrom fs dir_exists dir_create
#' @export
radius_function <- function(
    kvadrati_path,
    radii_path,
    tikls100_path,
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
    gdal_opts = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER",
                  "NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256"),
    write_datatype = NULL,
    NAflag = NULL,
    terra_memfrac = 0.7,
    terra_tempdir = tempdir(),
    terra_todisk = TRUE,
    quiet = FALSE
) {
  # ---- sink safety ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)
  say <- function(...) if (!quiet) cat(..., "\n")


  # ---- HPC/container safety ----
  options(future.fork.enable = FALSE)  # never fork

  # ---- deps ----
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("terra",        "raster IO & masks")
  .need_pkg("sf",           "zones handling")
  .need_pkg("raster",       "RasterLayer for fasterize grid")
  .need_pkg("dplyr",        "joins & wrangling")
  .need_pkg("tibble",       "small data frames")
  .need_pkg("purrr",        "lists & mapping")
  .need_pkg("fs",           "filesystem ops")
  .need_pkg("sfarrow",      "reading GeoParquet")
  .need_pkg("fasterize",    "polygon rasterization")
  .need_pkg("exactextractr","fast zonal statistics")
  .need_pkg("future",       "parallel execution")
  .need_pkg("furrr",        "parallel mapping")
  .need_pkg("whitebox",     "gap filling")
  .need_pkg("glue",         "paths & messages")
  .need_pkg("tidyr",        "file name parsing")

  # ---- terra options (set & restore) ----
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(memfrac = old_opt$memfrac,
                              tempdir = old_opt$tempdir,
                              todisk  = old_opt$todisk,
                              progress = old_opt$progress), add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir, progress = FALSE)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
  })


  # ---- basic checks ----
  if (!fs::dir_exists(output_dir)) fs::dir_create(output_dir, recurse = TRUE)
  stopifnot(length(input_layers) == length(layer_prefixes))
  names(input_layers) <- layer_prefixes
  if (length(extract_fun) != 1) stop("`extract_fun` must be a single function or string (e.g., \"mean\").")

  # ---- parallel plan (auto; PSOCK on HPC, multisession elsewhere) ----
  # Cap workers to SLURM allocation
  slurm_ct <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  slurm_on <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_ON_NODE")))
  cap <- max(c(slurm_ct, slurm_on), na.rm = TRUE)
  if (is.finite(cap) && cap > 0L) n_workers <- min(n_workers, cap)

  old_plan <- future::plan()
  on.exit(try(future::plan(old_plan), silent = TRUE), add = TRUE)

  old_max <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = future_max_size)
  on.exit(options(future.globals.maxSize = old_max), add = TRUE)

  is_hpc <- function() {
    any(nzchar(Sys.getenv(c(
      "SLURM_JOB_ID","PBS_JOBID","LSB_JOBID","APPTAINER_NAME","SINGULARITY_NAME"
    ))))
  }

  show_progress <- (!quiet) && !is_hpc()

  scratch <- Sys.getenv("SLURM_TMPDIR", unset = tempdir())
  Sys.setenv(
    OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1", BLIS_NUM_THREADS="1", GOTO_NUM_THREADS="1",
    RCPP_PARALLEL_NUM_THREADS="1", GDAL_NUM_THREADS="1",
    MALLOC_ARENA_MAX="2", GDAL_CACHEMAX="256", CPL_VSIL_CURL_CACHE_SIZE="0",
    TMPDIR=scratch, TEMP=scratch, R_TEMPORARY_DIR=scratch,CPL_TMPDIR = scratch
  )

  if (n_workers <= 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }

  # ---- list tiles & radii files ----
  kvadratiem <- tibble::tibble(fails_c = list.files(kvadrati_path, full.names = TRUE)) |>
    dplyr::mutate(cels_c = .data$fails_c, numurs = substr(basename(.data$fails_c), 10, 13))

  kv_rad <- tibble::tibble(fails_r = list.files(radii_path, full.names = TRUE)) |>
    dplyr::mutate(cels_radiuss = .data$fails_r, filename = basename(.data$fails_r)) |>
    tidyr::separate(.data$filename, into = c("sakums","veids","lapa","beigas"),remove = FALSE)

  if (identical(radius_mode, "sparse")) {
    kv_rad <- dplyr::filter(
      kv_rad,
      (.data$sakums == "pts100"  & .data$veids %in% c("r500","r1250")) |
        (.data$sakums == "pts300"  & .data$veids == "r3000") |
        (.data$sakums %in% c("pts1000","pts1km") & .data$veids == "r10000")
    )
  } else if (identical(radius_mode, "dense")) {
    kv_rad <- dplyr::filter(kv_rad, .data$sakums == "pts100" & .data$veids %in% c("r500","r1250","r3000","r10000"))
  } else stop("Invalid `radius_mode`. Use 'sparse' or 'dense'.")

  kv_rad <- tidyr::pivot_wider(
    kv_rad,
    id_cols   = lapa,
    names_from  = veids,
    values_from = c(fails_r, cels_radiuss),
    names_glue  = "{.value}_{veids}"
  ) |>
    dplyr::rename_with(~ gsub("fails_r_", "fails_", .x), dplyr::starts_with("fails_r_")) |>
    dplyr::rename_with(~ gsub("cels_radiuss_", "cels_", .x), dplyr::starts_with("cels_radiuss_"))

  kvadrati <- dplyr::left_join(kvadratiem, kv_rad, by = dplyr::join_by(numurs == lapa))
  kvadrati_list <- split(kvadrati, kvadrati$numurs)

  # ---- template (final alignment) ----
  template_r_full <- terra::rast(template_path)

  # ---- GDAL options (avoid oversubscription when parallel) ----
  tuned_defaults <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","BLOCKXSIZE=256","BLOCKYSIZE=256")
  tuned_defaults <- c(tuned_defaults, if (n_workers > 1L) "NUM_THREADS=1" else "NUM_THREADS=ALL_CPUS")

  # ensure only one NUM_THREADS setting
  gdal_opts <- c(setdiff(gdal_opts, gdal_opts[grepl("^NUM_THREADS=", gdal_opts)]), tuned_defaults)
  gdal_opts <- unique(gdal_opts)

  # ---- atomic write helper ----
  .write_r <- function(r, path, gdal_opts, write_datatype, NAflag) {
    ext <- tools::file_ext(path); ext <- if (nzchar(ext)) ext else "tif"
    stem <- sub(sprintf("\\.%s$", ext), "", path)
    tmp  <- sprintf("%s._tmp.%s", stem, ext)
    if (file.exists(tmp)) try(unlink(tmp), silent = TRUE)
    args <- list(x = r, filename = tmp, overwrite = TRUE, gdal = gdal_opts)
    if (!is.null(write_datatype)) args$datatype <- write_datatype
    if (!is.null(NAflag))        args$NAflag   <- NAflag
    do.call(terra::writeRaster, args)
    if (file.exists(path)) try(unlink(path), silent = TRUE)
    file.rename(tmp, path)
  }

  # ---- per-worker cache for heavy opens ----
  .get_cached <- function(name, loader) {
    x <- .egv_cache_get(name)
    if (!is.null(x)) return(x)
    x <- loader()
    .egv_cache_set(name, x)
    x
  }


  # ---- parser for exact_extract output (robust to list/data.frame/vector) ----
  .parse_extract <- function(res, layer_names, fun_name_or_fn, n_features) {
    # returns named list: layer -> numeric vector length n_features
    # case 1: single numeric vector (single layer)
    if (is.numeric(res)) {
      return(setNames(list(as.numeric(res)), layer_names[1]))
    }
    # case 2: data.frame with columns per layer (maybe prefixed like "mean.Layer")
    if (is.data.frame(res)) {
      cn <- colnames(res)
      fun_str <- if (is.character(fun_name_or_fn)) tolower(fun_name_or_fn) else NA_character_
      stripped <- if (!is.na(fun_str)) sub(paste0("^", fun_str, "\\."), "", cn) else cn
      if (!all(stripped %in% layer_names))
        stripped <- sub("^.*\\.", "", cn)   # last token after dot
      final_names <- if (all(stripped %in% layer_names)) stripped else cn
      keep_idx <- final_names %in% layer_names
      if (!any(keep_idx)) stop("Could not match extracted column names to layer names.")
      res2 <- res[, keep_idx, drop = FALSE]
      final_names <- final_names[keep_idx]
      out <- setNames(vector("list", length(final_names)), final_names)
      for (i in seq_along(final_names)) out[[i]] <- as.numeric(res2[[i]])
      return(out)
    }
    # case 3: list (one element per feature)
    if (is.list(res)) {
      # elements may be numeric vectors length = n_layers OR 1-row data.frames
      if (length(res) != n_features) {
        stop("exact_extract() returned a list of unexpected length: ", length(res), " (features: ", n_features, ")")
      }
      # normalize each element to a named numeric vector of length n_layers
      norm_one <- function(el) {
        if (is.numeric(el)) {
          v <- as.numeric(el)
          if (length(v) == 1L) {
            # single-layer case in list form
            names(v) <- layer_names[1]
            return(v)
          } else if (length(v) == length(layer_names)) {
            names(v) <- layer_names
            return(v)
          } else {
            stop("List element has length ", length(v), " but layers = ", length(layer_names))
          }
        } else if (is.data.frame(el)) {
          if (nrow(el) != 1L) stop("List element is a data.frame with ", nrow(el), " rows (expected 1).")
          cn <- colnames(el)
          cn2 <- sub("^.*\\.", "", cn)  # handle "mean.Layer"
          if (all(cn2 %in% layer_names)) {
            v <- as.numeric(el[1, , drop = TRUE]); names(v) <- cn2; return(v)
          } else if (length(cn) == length(layer_names)) {
            v <- as.numeric(el[1, , drop = TRUE]); names(v) <- layer_names; return(v)
          } else {
            stop("Cannot map data.frame columns to layer names in list element.")
          }
        } else {
          stop("Unsupported list element type from exact_extract(): ", class(el)[1])
        }
      }
      mat <- lapply(res, norm_one)                      # list of named numeric vectors
      # bind by row (features) then columns ordered by layer_names
      M <- do.call(rbind, lapply(mat, function(v) setNames(v[layer_names], layer_names)))
      out <- setNames(vector("list", length(layer_names)), layer_names)
      for (i in seq_along(layer_names)) out[[i]] <- as.numeric(M[, i])
      return(out)
    }
    stop("Unsupported result type from exact_extract(): ", class(res)[1])
  }

  # ---- per-tile worker (extract once per radius) ----
  process_tile <- function(solis, kv_row,
                           tikls100_path, template_path,
                           input_layers, layer_prefixes,
                           output_dir, radii, extract_fun,
                           gdal_opts, write_datatype, NAflag,
                           terra_tempdir) {

    say("Processing tile: ", solis)

    # cache heavy opens ONCE per worker
    template_100 <- .get_cached(".tmpl_cache", function() terra::rast(template_path))
    tikls100     <- .get_cached(".tikls100_cache", function() sfarrow::st_read_parquet(tikls100_path))
    stack_all    <- .get_cached(".stack_cache", function() {
      rs <- terra::rast(input_layers); names(rs) <- layer_prefixes; rs
    })

    # 1ha polygons for this tile
    sunas <- sfarrow::st_read_parquet(kv_row$cels_c)

    # available radius polygons
    r_polys <- purrr::map(purrr::set_names(radii), function(rad) {
      path <- kv_row[[paste0("cels_", rad)]]
      if (!is.null(path) && file.exists(path)) {
        vec <- tryCatch(sfarrow::st_read_parquet(path), error = function(e) NULL)
        if (!inherits(vec, "sf") || is.null(sf::st_geometry(vec))) NULL else vec
      } else NULL
    })

    # bbox of the largest available radius + 1000 m buffer
    largest <- NULL
    for (rad in rev(radii)) if (!is.null(r_polys[[rad]])) { largest <- r_polys[[rad]]; break }
    if (is.null(largest)) stop("No valid radius polygon for tile ", solis)

    telpa2        <- sf::st_buffer(sf::st_as_sfc(sf::st_bbox(largest)), dist = 1000)
    template_crop <- terra::crop(template_100, telpa2)
    templateRL    <- raster::raster(template_crop)  # grid for fasterize
    stack_crop    <- terra::crop(stack_all, telpa2)

    # parent grids for 3 km and 10 km joins
    sunas300 <- if (!is.null(r_polys$r3000)) {
      vals <- r_polys$r3000$rinda300; dplyr::filter(tikls100, .data$rinda300 %in% vals)
    } else NULL
    sunas1000 <- if (!is.null(r_polys$r10000)) {
      vals <- r_polys$r10000$ID1km;  dplyr::filter(tikls100, .data$ID1km %in% vals)
    } else NULL

    # ---- extract ONCE per radius over the (cropped) stack ----
    for (rad in radii) {
      vec <- r_polys[[rad]]
      if (is.null(vec)) next

      vec <- vec[!sf::st_is_empty(vec), , drop = FALSE]
      if (nrow(vec) == 0) next

      res <- suppressWarnings(
        exactextractr::exact_extract(stack_crop, vec, fun = extract_fun)
      )
      vals_by_layer <- .parse_extract(
        res, layer_names = names(stack_crop),
        fun_name_or_fn = extract_fun, n_features = nrow(vec)
      )

      # choose join base and key
      join_key <- dplyr::case_when(
        rad == "r500"   ~ "id",
        rad == "r1250"  ~ "id",
        rad == "r3000"  ~ "rinda300",
        rad == "r10000" ~ "ID1km"
      )
      joined_base <- switch(rad, r500 = sunas, r1250 = sunas, r3000 = sunas300, r10000 = sunas1000)
      if (is.null(joined_base)) next

      # For each layer: build xdf and rasterize
      for (prefix in names(vals_by_layer)) {
        v <- vals_by_layer[[prefix]]
        if (length(v) != nrow(vec)) stop("Row mismatch for extracted values vs polygons at ", rad, " layer ", prefix)

        xdf <- sf::st_drop_geometry(vec)
        xdf <- dplyr::select(xdf, dplyr::all_of(join_key))
        xdf$vertibas <- as.numeric(v)

        joined <- dplyr::left_join(joined_base, xdf, by = join_key)

        rr <- fasterize::fasterize(joined, templateRL, field = "vertibas")
        rr <- terra::rast(rr)
        rr <- terra::mask(rr, template_crop)

        out_path <- glue::glue("{output_dir}/{prefix}_{rad}/{prefix}_{rad}_{solis}.tif")
        dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
        .write_r(rr, out_path, gdal_opts, write_datatype, NAflag)
      }
    }

    return(solis)
  }

  # ---- run tiles in parallel ----
  opts <- furrr::furrr_options(
    seed    = TRUE,
    globals = c(
      "process_tile",
      ".get_cached",
      ".parse_extract",
      ".write_r",
      "say",
      ".egvtools_cache",
      ".egv_cache_get",
      ".egv_cache_set",
      ".egv_cache_drop",
      ".egv_cache_clear"
    ),
    packages = c("terra","sf","sfarrow","fasterize","exactextractr",
                 "dplyr","tibble","tidyr","glue","raster")
  )


  furrr::future_map2(
    .x = names(kvadrati_list),
    .y = kvadrati_list,
    .f = function(x, y) process_tile(x, y,
                                     tikls100_path, template_path,
                                     input_layers, layer_prefixes,
                                     output_dir, radii, extract_fun,
                                     gdal_opts, write_datatype, NAflag,
                                     terra_tempdir),
    .progress = show_progress,
    .options = opts
  )

  # ---- merge tiles then project/mask then gap analysis then optional fill then write ----
  template_r <- template_r_full
  pix_size   <- mean(terra::res(template_r))
  results_ls <- list()

  for (prefix in layer_prefixes) {
    for (rad in radii) {
      tif_dir   <- glue::glue("{output_dir}/{prefix}_{rad}")
      out_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)
      if (!length(out_files)) next

      vrt <- terra::vrt(out_files)
      names(vrt) <- glue::glue("{prefix}_{rad}")

      temp_raw <- file.path(terra_tempdir, glue::glue("mosaic_{prefix}_{rad}.tif"))
      .write_r(vrt, temp_raw, gdal_opts, write_datatype, NAflag)

      proj_r <- terra::project(terra::rast(temp_raw), template_r)
      proj_r <- terra::mask(proj_r, template_r)

      robi <- is.na(proj_r) & !is.na(template_r)
      gaps_before <- as.integer(terra::global(terra::as.int(robi), fun = "sum", na.rm = TRUE)[[1]])
      if (is.na(gaps_before)) gaps_before <- 0L

      gap_filled <- FALSE
      filter_used <- NA_integer_

      if (isTRUE(fill_missing) && gaps_before > 0L) {
        say(glue::glue("[{prefix}_{rad}] Filling {format(gaps_before, big.mark=',')} NA cells"))
        fillable <- terra::ifel(!robi, 1, NA)
        dist_r   <- terra::distance(fillable)
        max_att  <- suppressWarnings(terra::global(dist_r, fun = "max", na.rm = TRUE)[[1]])
        rm(fillable, dist_r)

        filter_used <- max(3L, as.integer(ceiling(as.numeric(max_att) / pix_size) * 2L))  # min clamp only
        if (filter_used %% 2 == 0) filter_used <- filter_used + 1L

        temp_filled <- file.path(terra_tempdir, glue::glue("filled_{prefix}_{rad}.tif"))
        ok <- TRUE
        tryCatch({
          whitebox::wbt_fill_missing_data(
            input  = temp_raw,
            output = temp_filled,
            filter = filter_used,
            weight = IDW_weight,
            no_edges = FALSE
          )
        }, error = function(e) {
          say(glue::glue("[{prefix}_{rad}] Whitebox failed: {conditionMessage(e)}"))
          ok <<- FALSE
        })

        if (ok && file.exists(temp_filled)) {
          proj_r <- terra::project(terra::rast(temp_filled), template_r)
          proj_r <- terra::mask(proj_r, template_r)
          terra::readAll(proj_r)
          unlink(temp_filled)
          gap_filled <- TRUE
        }
      }

      robi2 <- is.na(proj_r) & !is.na(template_r)
      gaps_after <- as.integer(terra::global(terra::as.int(robi2), fun = "sum", na.rm = TRUE)[[1]])
      if (is.na(gaps_after)) gaps_after <- 0L

      names(proj_r) <- glue::glue("{prefix}_{rad}")

      terra::crs(proj_r) <- terra::crs(template_r, proj = FALSE)

      out_final <- glue::glue("{output_dir}/{prefix}_{rad}.tif")
      .write_r(proj_r, out_final, gdal_opts, write_datatype, NAflag)

      if (file.exists(temp_raw)) unlink(temp_raw)
      if (unlink_tiles) try(unlink(tif_dir, recursive = TRUE), silent = TRUE)

      results_ls[[length(results_ls)+1L]] <- data.frame(
        layer            = prefix,
        radius           = rad,
        output_path      = out_final,
        n_tiles_merged   = length(out_files),
        gaps_before      = gaps_before,
        gaps_after       = gaps_after,
        filter_size_used = if (isTRUE(fill_missing)) filter_used else NA_integer_,
        gap_filled       = isTRUE(gap_filled),
        stringsAsFactors = FALSE
      )
    }
  }

  res_df <- if (length(results_ls)) do.call(rbind, results_ls) else
    data.frame(layer=character(), radius=character(), output_path=character(),
               n_tiles_merged=integer(), gaps_before=integer(), gaps_after=integer(),
               filter_size_used=integer(), gap_filled=logical())

  invisible(res_df)
}
