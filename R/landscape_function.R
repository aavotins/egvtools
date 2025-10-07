#' Compute landscape-level metrics per zone (tiled), merge, analyze gaps, optionally fill with IDW, and write rasters
#'
#' @description
#' Computes a \pkg{landscapemetrics} metric (default `"lsm_l_shdi"`), optionally with
#' extra `lm_args`, that yields one value per zone and **per input layer**. Runs tile-by-tile
#' (by `tile_field`), writes per-tile rasters, merges to final per-layer GeoTIFF(s),
#' then performs **gap analysis** (NA count within the template footprint and optional
#' maximum gap width) and **optional IDW gap filling** via WhiteboxTools.
#' Returns a compact **data.frame** with per-layer stats and timing.
#'
#' @details
#' **Workflow**
#' 1. **Inputs & CRS**
#'    - Accept paths or in-memory (`SpatRaster`, `sf`). Zones and landscape are aligned to the **template** CRS/grid.
#' 2. **Tiling**
#'    - Split `zones` by `tile_field`; process each tile independently.
#' 3. **Per-tile (all layers)**
#'    - Quick extent crop (buffered by `buffer_m`), `landscapemetrics::sample_lsm()` per layer,
#'      join back by `id_field`, rasterize (engine: `"fasterize"` or `"terra"`), mask to template crop,
#'      and write tile raster. Skips empty crops/joins cleanly.
#' 4. **Merge per layer**
#'    - VRT over all tile rasters to final GeoTIFF with tuned GDAL options (LZW, tiling, BIGTIFF=IF_SAFER, threads).
#'    - CRS is **forced to the template WKT** so `crs(output) == crs(template)` string-matches.
#' 5. **Gap analysis & optional fill (per layer)**
#'    - `gap_count`: number of NA cells **inside** the template footprint.
#'    - If `report_gap_size=TRUE` or `filter_size_cells="auto"`, compute **max gap width** (distance to nearest non-NA).
#'    - If `fill_gaps=TRUE` and gaps exist, fill with `whitebox::wbt_fill_missing_data()` using
#'      `idw_weight` and `filter_size_cells` (or auto from the max gap width; odd, >=3).
#'    - Plot result and/or gaps if requested (side-by-side if both).
#' 6. **Return**
#'    - One row per output layer with paths, tile counters, gap stats, fill parameters, and elapsed seconds.
#'
#' @param landscape SpatRaster or path(s) to .tif (class labels). Multi-layer supported.
#' @param zones sf polygons or path to \pkg{sfarrow} GeoParquet. Default `"./Templates/TemplateGrids/tikls500_sauzeme.parquet"`.
#' @param id_field Zone id field (default `"rinda500"`).
#' @param tile_field Tiling field (default `"tks50km"`).
#' @param template SpatRaster or path defining target grid/CRS/cellsize. Default `"./Templates/TemplateRasters/LV500m_10km.tif"`.
#' @param out_dir Output directory (created if missing).
#' @param out_filename Character vector: final filename(s) (no dir); one per input layer, matching order.
#' @param out_layername Character vector: final layer name(s); one per input layer, matching order.
#' @param what Landscapemetrics metric (default `"lsm_l_shdi"`).
#' @param lm_args Named list of extra args for [landscapemetrics::sample_lsm()] compatible with `what`.
#' @param buffer_m Numeric buffer (m) for tile bbox before crop (default `1000`).
#' @param rasterize_engine `"fasterize"` (default) or `"terra"`.
#' @param n_workers Integer; `1` = sequential.
#' @param os_type `"auto","windows","mac","linux","slurm"` (default `"auto"`).
#' @param future_max_size Max globals per worker (bytes); default `4 * 1024^3` (~4 GiB).
#' @param gdal_opts GDAL creation options (merged with tuned defaults:
#'   `c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")`).
#' @param write_datatype terra datatype string; if `NULL`, defaults to `"FLT4S"`.
#' @param NAflag Numeric NA value to write; if `NA`/`NULL`, a sensible default is chosen
#'   (`FLT* to -9999`, `INT2S to -32768`, `INT4S to -2147483648`).
#' @param keep_tiles Logical; keep temp tiles dir (debug). Default `FALSE`.
#' @param skip_existing Logical; skip tiles/finals that already exist. Default `TRUE`.
#' @param terra_memfrac `terraOptions(memfrac=...)`. Default `0.7`.
#' @param terra_tempdir Temp dir for terra ops. Default `tempdir()`.
#' @param terra_todisk Logical or `NA`. If `TRUE`, prefer on-disk. Default `TRUE`.
#' @param force_gc Logical; call `gc()` at key checkpoints. Default `FALSE`.
#' @param quiet Suppress progress prints (`cat()`)? Default `FALSE`.
#' @param report_gaps Logical; print/compute `gap_count`. Default `TRUE`.
#' @param report_gap_size Logical; compute/print maximum gap width (expensive). Default `TRUE`.
#' @param fill_gaps Logical; run Whitebox IDW fill if gaps exist. Default `FALSE`.
#' @param idw_weight Numeric; IDW power for Whitebox fill. Default `2`.
#' @param filter_size_cells Integer or `"auto"`; Whitebox window size (cells). Default `"auto"`.
#' @param plot_result Logical; plot final per-layer raster(s). Default `FALSE`.
#' @param plot_gaps Logical; plot gap map(s) (1 = NA inside template). Default `FALSE`.
#'
#' @return A **data.frame** (one row per layer) with:
#'   - `layer_name`, `output_path`,
#'   - `tiles_written`, `tiles_skipped_existing`, `tiles_skipped_empty_crop`, `tiles_skipped_empty_join`,
#'   - `merge_skipped`, `gap_count`, `max_gap_distance`, `filter_size_cells_used`, `gap_filled`,
#'   - `n_tiles`, `n_zones`, `n_layers`, `elapsed_sec`.
#'   Attributes: `"tile_dir"` (path or `NULL`), `"run_params"` (list).
#'
#' @examples
#' \dontrun{
#' res_tbl <- landscape_function(
#'   landscape      = "./Rastri_10m/Ainava_KopejaiDaudzveidibai.tif",
#'   zones          = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
#'   id_field       = "rinda500",
#'   tile_field     = "tks50km",
#'   template       = "./Templates/TemplateRasters/LV500m_10km.tif",
#'   out_dir        = "./out/lsm/",
#'   out_filename   = "Landscape_diversity.tif",
#'   out_layername  = "Landscape_diversity",
#'   what           = "lsm_l_shdi",
#'   rasterize_engine = "fasterize",
#'   n_workers      = 8,
#'   fill_gaps      = TRUE,
#'   plot_gaps      = TRUE,
#'   plot_result    = TRUE
#' )
#' print(res_tbl)
#'
#' #' ## --- Total edge length per zone (ignore outside map boundary) ---
#' ## Using a binary "water vs other" landscape, compute lsm_l_te per 500 m zone.
#' rez_edges <- landscape_function(
#'   landscape        = "./Rastri_10m/Ainava_vienk_mask.tif",
#'   zones            = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
#'   id_field         = "rinda500",
#'   tile_field       = "tks50km",
#'   template         = "./Templates/TemplateRasters/LV500m_10km.tif",
#'   out_dir          = "./",
#'   out_filename     = "edges_water.tif",
#'   out_layername    = "edges_water",
#'   what             = "lsm_l_te",
#'   lm_args          = list(count_boundary = FALSE),
#'   rasterize_engine = "terra",
#'   n_workers        = 8,
#'   future_max_size  = 2 * 1024^3,
#'   report_gaps      = TRUE,
#'   plot_result      = TRUE,
#'   plot_gaps        = TRUE
#' )
#' rez_edges
#' }
#'
#' @seealso landscapemetrics::sample_lsm, terra::writeRaster, fasterize::fasterize, whitebox::wbt_fill_missing_data
#' @import terra
#' @import sf
#' @import landscapemetrics
#' @importFrom sfarrow st_read_parquet
#' @importFrom stringr str_ends str_detect str_replace_all fixed
#' @importFrom dplyr left_join mutate filter
#' @importFrom tibble tibble
#' @importFrom purrr map imap
#' @importFrom fs path
#' @importFrom tools file_ext
#' @importFrom utils capture.output modifyList
#' @importFrom future plan sequential multisession multicore cluster
#' @importFrom furrr furrr_options future_imap
#' @importFrom raster stack projection raster ncell
#' @importFrom whitebox wbt_fill_missing_data
#' @importFrom methods as
#' @importFrom sp proj4string
#' @importFrom graphics par box
#' @importFrom stats setNames
#' @importFrom rlang .data :=
#' @export
landscape_function <- function(
    landscape,
    zones = "./Templates/TemplateGrids/tikls500_sauzeme.parquet",
    id_field   = "rinda500",
    tile_field = "tks50km",
    template   = "./Templates/TemplateRasters/LV500m_10km.tif",
    out_dir,
    out_filename,
    out_layername,
    what = "lsm_l_shdi",
    lm_args = NULL,
    buffer_m = 1000,
    rasterize_engine = c("fasterize","terra"),
    n_workers = 1,
    os_type = "auto",
    future_max_size = 4 * 1024^3,
    gdal_opts = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"),
    write_datatype = NULL,
    NAflag = NA_real_,
    keep_tiles = FALSE,
    skip_existing = TRUE,
    terra_memfrac = 0.7,
    terra_tempdir = tempdir(),
    terra_todisk = TRUE,
    force_gc = FALSE,
    quiet = FALSE,
    report_gaps = TRUE,
    report_gap_size = TRUE,
    fill_gaps = FALSE,
    idw_weight = 2,
    filter_size_cells = "auto",
    plot_result = FALSE,
    plot_gaps = FALSE
){
  t_start <- Sys.time()
  rasterize_engine <- match.arg(rasterize_engine)

  # ---- sink safety
  orig_out <- sink.number(); orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)

  say <- function(...) if (!quiet) cat(..., "\n")
  `%||%` <- function(x, y) if (is.null(x)) y else x
  .maybe_gc <- function() if (isTRUE(force_gc)) gc()

  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("terra", "raster IO & masks")
  .need_pkg("sf",    "zones handling")
  .need_pkg("raster","multi-layer stacks for landscapemetrics")
  .need_pkg("landscapemetrics", "metric computation")
  .need_pkg("stringr","tile/layer filename handling")
  .need_pkg("dplyr", "joins")
  .need_pkg("tibble","small data frames")
  .need_pkg("purrr", "lists & mapping")
  .need_pkg("fs",    "filesystem ops")
  .need_pkg("sfarrow", "reading GeoParquet 'zones'")
  .need_pkg("fasterize", "polygon rasterization")
  .need_pkg("future", "parallel execution")
  .need_pkg("furrr",  "parallel mapping")
  .need_pkg("whitebox",    "gap filling")


  # terra options (set + restore)
  old_opt <- NULL
  utils::capture.output({ old_opt <- terra::terraOptions() })
  on.exit(terra::terraOptions(
    memfrac = old_opt$memfrac, tempdir = old_opt$tempdir,
    todisk = old_opt$todisk, progress = old_opt$progress
  ), add = TRUE)
  utils::capture.output({
    terra::terraOptions(memfrac = terra_memfrac, tempdir = terra_tempdir, progress = FALSE)
    if (!is.na(terra_todisk)) terra::terraOptions(todisk = isTRUE(terra_todisk))
  })

  # helpers
  load_landscape <- function(x) if (inherits(x,"SpatRaster")) x else terra::rast(x)
  load_zones <- function(z) if (inherits(z,"sf")) z else {
    zlow <- tolower(z)
    if (stringr::str_ends(zlow, "\\.parquet$")) sfarrow::st_read_parquet(z, as_tibble = FALSE, quiet = TRUE) else sf::st_read(z, quiet = TRUE)
  }
  sanitize_fs <- function(x) x |>
    stringr::str_replace_all("[^A-Za-z0-9._-]+", "_") |>
    stringr::str_replace("^_+", "") |>
    stringr::str_replace("_+$", "")
  on_slurm <- function() any(nzchar(Sys.getenv(c("SLURM_JOB_ID","SLURM_CPUS_ON_NODE","SLURM_JOB_NODELIST"))))
  safe_plan <- function(n_workers=1, os_type="auto"){
    if (n_workers==1) return(future::plan(future::sequential))
    osd <- if (identical(tolower(os_type),"auto")) { if (on_slurm()) "slurm" else tolower(Sys.info()[["sysname"]]) } else tolower(os_type)
    if (osd=="slurm") { requireNamespace("parallel", quietly=TRUE); cl <- parallel::makeCluster(n_workers); return(future::plan(future::cluster, workers=cl)) }
    if (osd %in% c("windows","mac","darwin")) return(future::plan(future::multisession, workers=n_workers))
    if (osd=="linux") { tryCatch(future::plan(future::multicore, workers=n_workers),
                                 error=function(e) future::plan(future::multisession, workers=n_workers)); return(invisible()) }
    stop("Unknown os_type: ", osd)
  }
  ensure_raster_path <- function(r, tmp_dir, label){
    if (inherits(r,"SpatRaster")){
      src <- try(terra::sources(r), silent=TRUE)
      if (!inherits(src,"try-error") && length(src)>0 && all(file.exists(src))) return(src)
      p <- file.path(tmp_dir, paste0("tmp_", label, ".tif")); terra::writeRaster(r, p, overwrite=TRUE); return(p)
    } else { stopifnot(is.character(r), length(r)>=1, all(file.exists(r))); return(r) }
  }
  valid_spat <- function(x) inherits(x,"SpatRaster") && !inherits(x,"try-error") && !is.null(x) && terra::ncell(x) > 0
  safe_datatype <- function(dt) if (is.null(dt) || is.na(dt)) "FLT4S" else toupper(dt)
  infer_default_naflag <- function(datatype) {
    if (is.null(datatype) || is.na(datatype)) return(NA_real_)
    dt <- toupper(datatype)
    if (startsWith(dt, "FLT")) return(-9999)
    if (dt == "INT2S")        return(-32768)
    if (dt == "INT4S")        return(-2147483648)
    NA_real_
  }

  # inputs & checks
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  r_in <- load_landscape(landscape)
  z_in <- load_zones(zones)
  tmpl <- if (inherits(template,"SpatRaster")) template else terra::rast(template)

  if (!id_field   %in% names(z_in)) stop("id_field '", id_field,   "' not found in zones.")
  if (!tile_field %in% names(z_in)) stop("tile_field '", tile_field, "' not found in zones.")

  # CRS harmonization (WKT equality)
  tmpl_wkt <- terra::crs(tmpl, proj = TRUE)
  if (!identical(sf::st_crs(z_in)$wkt, tmpl_wkt)) z_in <- sf::st_transform(z_in, tmpl_wkt)
  if (!identical(terra::crs(r_in, proj = TRUE), tmpl_wkt)) r_in <- terra::project(r_in, tmpl, method = "near")

  nlyr <- terra::nlyr(r_in)
  if (nlyr==1L) { stopifnot(length(out_filename)==1L, length(out_layername)==1L) } else {
    stopifnot(length(out_filename)==nlyr, length(out_layername)==nlyr)
  }

  # tiles & futures
  tile_root <- file.path(out_dir, paste0("tiles_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(tile_root, recursive = TRUE, showWarnings = FALSE)

  old_plan <- future::plan()
  old_max  <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = future_max_size)
  on.exit({
    options(future.globals.maxSize = old_max)
    try(future::plan(old_plan), silent = TRUE)
    if (!keep_tiles) try(unlink(tile_root, recursive = TRUE, force = TRUE), silent = TRUE)
  }, add = TRUE)
  safe_plan(n_workers = n_workers, os_type = os_type)

  landscape_path <- ensure_raster_path(r_in, tile_root, "landscape")
  template_path  <- ensure_raster_path(tmpl,  tile_root, "template")
  tiles_sf <- split(z_in, f = z_in[[tile_field]], drop = TRUE)

  layer_counts <- data.frame(
    layer_name = out_layername,
    tiles_written = 0L,
    tiles_skipped_existing = 0L,
    tiles_skipped_empty_crop = 0L,
    tiles_skipped_empty_join = 0L,
    stringsAsFactors = FALSE
  )

  # core tile worker
  process_one_tile <- function(tile_sf_in, tile_id){
    tryCatch({
      if (inherits(tile_sf_in, "list")) tile_sf_in <- tile_sf_in[[1]]
      if (!inherits(tile_sf_in, "sf"))  tile_sf_in <- sf::st_as_sf(tile_sf_in)
      safe_id  <- sanitize_fs(as.character(tile_id))
      if (!quiet) cat("Tile", safe_id, ": start (n=", nrow(tile_sf_in), ")\n", sep = "")

      lc <- setNames(as.list(rep(0L, nlyr)), paste0("w_", out_layername))
      se <- setNames(as.list(rep(0L, nlyr)), paste0("se_", out_layername))
      sc <- setNames(as.list(rep(0L, nlyr)), paste0("scrop_", out_layername))
      sj <- setNames(as.list(rep(0L, nlyr)), paste0("sjoin_", out_layername))

      expected_files <- file.path(tile_root, paste0("tile_", safe_id, "_", sanitize_fs(out_layername), ".tif"))
      if (isTRUE(skip_existing) && length(expected_files) > 0 && all(file.exists(expected_files))) {
        if (!quiet) cat("Tile", safe_id, ": all per-layer tiles exist; skip compute\n")
        for (i in seq_len(nlyr)) se[[paste0("se_", out_layername[i])]] <- se[[paste0("se_", out_layername[i])]] + 1L
        return(list(lc=lc,se=se,sc=sc,sj=sj))
      }

      tile_sf <- sf::st_make_valid(tile_sf_in)
      if (nrow(tile_sf)==0 || all(sf::st_is_empty(tile_sf))) return(list(lc=lc,se=se,sc=sc,sj=sj))
      tile_sf[[id_field]] <- as.character(tile_sf[[id_field]])

      r_in  <- terra::rast(landscape_path)
      tmpl  <- terra::rast(template_path)

      sv   <- terra::vect(tile_sf); e <- terra::ext(sv)
      ebuf <- terra::ext(e[1]-buffer_m, e[2]+buffer_m, e[3]-buffer_m, e[4]+buffer_m)

      tmpl_crop <- try(terra::crop(tmpl, ebuf), silent = TRUE)
      if (inherits(tmpl_crop,"try-error") || is.null(tmpl_crop) || terra::ncell(tmpl_crop)==0L) {
        for (i in seq_len(nlyr)) sc[[paste0("scrop_", out_layername[i])]] <- sc[[paste0("scrop_", out_layername[i])]] + 1L
        return(list(lc=lc,se=se,sc=sc,sj=sj))
      }
      r_crop <- try(terra::crop(r_in, ebuf), silent = TRUE)
      if (inherits(r_crop,"try-error") || is.null(r_crop) || terra::ncell(r_crop)==0L) {
        for (i in seq_len(nlyr)) sc[[paste0("scrop_", out_layername[i])]] <- sc[[paste0("scrop_", out_layername[i])]] + 1L
        return(list(lc=lc,se=se,sc=sc,sj=sj))
      }

      r_crop_path <- file.path(tile_root, paste0("landscape_", safe_id, ".tif"))
      terra::writeRaster(r_crop, r_crop_path, overwrite = TRUE)
      r_crop_r <- raster::stack(r_crop_path)

      crs_proj4 <- terra::crs(r_crop, proj = FALSE)
      if (!is.null(crs_proj4) && !is.na(crs_proj4) && nzchar(crs_proj4)) raster::projection(r_crop_r) <- crs_proj4
      tile_sp <- methods::as(tile_sf, "Spatial")
      if (!is.null(crs_proj4) && !is.na(crs_proj4) && nzchar(crs_proj4)) suppressWarnings(sp::proj4string(tile_sp) <- crs_proj4)

      base_args <- list(landscape = r_crop_r, y = tile_sp,
                        plot_id = as.character(tile_sf[[id_field]]), what = what)
      call_args <- if (is.null(lm_args)) base_args else utils::modifyList(base_args, lm_args, keep.null = TRUE)
      lsm_tb <- try(suppressWarnings(do.call(landscapemetrics::sample_lsm, call_args)), silent = TRUE)
      if (inherits(lsm_tb, "try-error") || is.null(lsm_tb) || !all(c("value","plot_id") %in% names(lsm_tb))) {
        for (i in seq_len(nlyr)) sj[[paste0("sjoin_", out_layername[i])]] <- sj[[paste0("sjoin_", out_layername[i])]] + 1L
        return(list(lc=lc,se=se,sc=sc,sj=sj))
      }
      lsm_tb$value[is.na(lsm_tb$value)] <- 0
      layer_vec <- if ("layer" %in% names(lsm_tb)) sort(unique(lsm_tb$layer)) else 1L

      dt_out <- safe_datatype(write_datatype)
      na_out <- if (is.null(NAflag) || is.na(NAflag)) infer_default_naflag(dt_out) else NAflag

      for (li in layer_vec) {
        sub_tb <- if (length(layer_vec)==1L && !("layer" %in% names(lsm_tb))) lsm_tb else dplyr::filter(lsm_tb, .data$layer==!!li)
        lyr_index <- if (length(layer_vec)==1L && !("layer" %in% names(lsm_tb))) 1L else li
        lyr_name  <- out_layername[lyr_index]
        safe_lyr  <- sanitize_fs(lyr_name)
        lyr_file  <- file.path(tile_root, paste0("tile_", safe_id, "_", safe_lyr, ".tif"))

        join_df <- tibble::tibble(!!id_field := as.character(sub_tb$plot_id), value = sub_tb$value)
        tile_sf_val <- dplyr::left_join(
          dplyr::mutate(tile_sf, !!id_field := as.character(.data[[id_field]])),
          join_df, by = id_field
        )
        if (!"value" %in% names(tile_sf_val) || all(is.na(tile_sf_val$value))) {
          sj[[paste0("sjoin_", lyr_name)]] <- sj[[paste0("sjoin_", lyr_name)]] + 1L
          next
        }

        if (isTRUE(skip_existing) && file.exists(lyr_file)) {
          se[[paste0("se_", lyr_name)]] <- se[[paste0("se_", lyr_name)]] + 1L
          next
        }

        rr <- NULL
        if (identical(rasterize_engine, "fasterize") && requireNamespace("fasterize", quietly = TRUE)) {
          e    <- as.vector(terra::ext(tmpl_crop)); rres <- terra::res(tmpl_crop)
          crs_p4 <- terra::crs(tmpl_crop, proj = FALSE)
          rtpl <- try({
            if (!is.null(crs_p4) && !is.na(crs_p4) && !identical(crs_p4, "")) {
              raster::raster(xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], res=rres, crs=crs_p4)
            } else raster::raster(xmn=e[1], xmx=e[2], ymn=e[3], ymx=e[4], res=rres)
          }, silent = TRUE)
          if (!inherits(rtpl,"try-error") && !is.null(rtpl) && raster::ncell(rtpl)>0) {
            tf <- try(sf::st_cast(tile_sf_val, "MULTIPOLYGON", warn = FALSE), silent = TRUE)
            if (!inherits(tf,"try-error")) tile_sf_val <- tf
            rr <- try({ rr_r <- fasterize::fasterize(tile_sf_val, rtpl, field = "value", fun = "last"); terra::rast(rr_r) }, silent = TRUE)
          }
          if (!valid_spat(rr)) {
            svv <- terra::vect(tile_sf_val)
            rr  <- try(terra::rasterize(x = svv, y = tmpl_crop, field = "value", fun = "mean", background = NA_real_), silent = TRUE)
          }
        } else {
          svv <- terra::vect(tile_sf_val)
          rr  <- try(terra::rasterize(x = svv, y = tmpl_crop, field = "value", fun = "mean", background = NA_real_), silent = TRUE)
        }

        if (!valid_spat(rr)) {
          sc[[paste0("scrop_", lyr_name)]] <- sc[[paste0("scrop_", lyr_name)]] + 1L
          next
        }
        names(rr) <- lyr_name
        rr <- terra::mask(rr, tmpl_crop)

        terra::writeRaster(rr, filename = lyr_file, overwrite = TRUE,
                           gdal = gdal_opts, datatype = dt_out, NAflag = na_out)
        lc[[paste0("w_", lyr_name)]] <- lc[[paste0("w_", lyr_name)]] + 1L
      }

      list(lc=lc,se=se,sc=sc,sj=sj)

    }, error = function(e){
      if (!quiet) cat("Tile", as.character(tile_id), ": hard error ->", conditionMessage(e), "\n")
      list(
        lc = setNames(as.list(rep(0L, nlyr)), paste0("w_", out_layername)),
        se = setNames(as.list(rep(0L, nlyr)), paste0("se_", out_layername)),
        sc = setNames(as.list(rep(0L, nlyr)), paste0("scrop_", out_layername)),
        sj = setNames(as.list(rep(0L, nlyr)), paste0("sjoin_", out_layername))
      )
    })
  }

  furrr_opts <- furrr::furrr_options(seed = TRUE, globals = FALSE)
  tile_stats <- furrr::future_imap(tiles_sf, ~process_one_tile(.x, .y),
                                   .options = furrr_opts, .progress = !quiet)

  # aggregate tile counters
  for (ts in tile_stats) {
    if (is.null(ts)) next
    for (i in seq_len(nlyr)) {
      ln <- out_layername[i]
      layer_counts$tiles_written[i]            <- layer_counts$tiles_written[i]            + (ts$lc[[paste0("w_", ln)]]    %||% 0L)
      layer_counts$tiles_skipped_existing[i]   <- layer_counts$tiles_skipped_existing[i]   + (ts$se[[paste0("se_", ln)]]   %||% 0L)
      layer_counts$tiles_skipped_empty_crop[i] <- layer_counts$tiles_skipped_empty_crop[i] + (ts$sc[[paste0("scrop_", ln)]]%||% 0L)
      layer_counts$tiles_skipped_empty_join[i] <- layer_counts$tiles_skipped_empty_join[i] + (ts$sj[[paste0("sjoin_", ln)]]%||% 0L)
    }
  }

  # per-layer outputs & gap stats
  all_tile_files <- list.files(tile_root, pattern="\\.tif$", full.names=TRUE)
  outputs <- setNames(rep(NA_character_, nlyr), out_layername)
  merge_skipped <- setNames(rep(FALSE, nlyr), out_layername)
  gap_count_vec <- rep(NA_real_, nlyr)
  max_gap_vec   <- rep(NA_real_, nlyr)
  fs_used_vec   <- rep(NA_integer_, nlyr)
  gap_filled_vec<- rep(FALSE, nlyr)

  def_opts <- c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER","NUM_THREADS=ALL_CPUS","BLOCKXSIZE=256","BLOCKYSIZE=256")
  gdal_final <- unique(c(gdal_opts, def_opts))
  dt_final   <- safe_datatype(write_datatype)
  na_final   <- if (is.null(NAflag) || is.na(NAflag)) infer_default_naflag(dt_final) else NAflag
  if (startsWith(dt_final, "FLT") && !any(grepl("^PREDICTOR=", gdal_final))) gdal_final <- c(gdal_final, "PREDICTOR=2")

  for (i in seq_len(nlyr)) {
    lyr_name   <- out_layername[i]
    safe_lyr   <- sanitize_fs(lyr_name)
    final_path <- file.path(out_dir, out_filename[i])

    lyr_paths <- all_tile_files[stringr::str_detect(
      all_tile_files, paste0("_", stringr::fixed(safe_lyr), "\\.tif$")
    )]
    if (length(lyr_paths) == 0) {
      outputs[[lyr_name]] <- final_path
      next
    }

    if (isTRUE(skip_existing) && file.exists(final_path)) {
      merge_skipped[[lyr_name]] <- TRUE
      outputs[[lyr_name]] <- final_path
      next
    }

    # --- merge tiles to a VRT, tag CRS
    vrt <- terra::vrt(lyr_paths); names(vrt) <- lyr_name
    terra::crs(vrt) <- tmpl_wkt

    # --- write the raw mosaic to a tmp (same dtype/opts as final)
    ext  <- tools::file_ext(final_path); ext <- if (nzchar(ext)) ext else "tif"
    stem <- sub(sprintf("\\.%s$", ext), "", final_path)
    tmp_merge <- sprintf("%s._merge.%s", stem, ext)
    if (file.exists(tmp_merge)) try(unlink(tmp_merge), silent = TRUE)

    terra::writeRaster(vrt, filename = tmp_merge, overwrite = TRUE,
                       gdal = gdal_final, datatype = dt_final, NAflag = na_final)

    # --- ALIGN to the *full* template grid (handles buffer around presence)
    r_mos <- terra::rast(tmp_merge)
    r_aln <- if (terra::same.crs(r_mos, tmpl)) {
      # same CRS: resample onto template extent/origin (fills outer area with NA)
      terra::resample(
        r_mos, tmpl, method = "near",
        filename = if (isTRUE(terra_todisk)) tempfile(fileext = ".tif") else "",
        overwrite = TRUE
      )
    } else {
      # CRS differs (shouldn't normally happen here, but guard anyway)
      terra::project(
        r_mos, tmpl, method = "near",
        filename = if (isTRUE(terra_todisk)) tempfile(fileext = ".tif") else "",
        overwrite = TRUE
      )
    }
    terra::crs(r_aln) <- tmpl_wkt  # ensure exact WKT equality

    # --- GAP ANALYSIS (extents match template)
    gaps <- is.na(r_aln) & !is.na(tmpl)
    gap_count <- if (isTRUE(report_gaps) || isTRUE(report_gap_size) ||
                     isTRUE(plot_gaps)   || (isTRUE(fill_gaps) && identical(filter_size_cells, "auto"))) {
      terra::global(gaps, fun = "sum", na.rm = TRUE)[[1]]
    } else NA_real_
    max_gap <- NA_real_

    if (!is.na(gap_count) && gap_count > 0 &&
        (isTRUE(report_gap_size) || (isTRUE(fill_gaps) && identical(filter_size_cells, "auto")))) {
      data_mask <- terra::ifel(gaps, NA, 1)
      dist_r <- terra::distance(data_mask)
      max_gap <- terra::global(dist_r, fun = "max", na.rm = TRUE)[[1]]
      rm(dist_r); .maybe_gc()
    }

    if (isTRUE(report_gaps))     say(sprintf("[%s] Gap cells (inside template): %s", lyr_name, format(gap_count, big.mark=",")))
    if (isTRUE(report_gap_size) && !is.na(max_gap)) say(sprintf("[%s] Maximum gap width: %.3f (template units)", lyr_name, max_gap))

    # --- Optional Whitebox IDW fill (on the aligned, template-sized grid)
    fs_used <- NA_integer_
    if (isTRUE(fill_gaps) && !is.na(gap_count) && gap_count > 0) {
      if (!requireNamespace("whitebox", quietly = TRUE)) {
        warning(sprintf("[%s] Package 'whitebox' not installed; skipping gap filling.", lyr_name))
      } else {
        if (identical(filter_size_cells, "auto")) {
          pix <- mean(terra::res(tmpl))
          fs  <- ceiling((max_gap / pix) * 2)
          if (fs %% 2 == 0) fs <- fs + 1
          fs_used <- max(3L, as.integer(fs))
        } else {
          fs_used <- max(3L, as.integer(filter_size_cells))
          if (fs_used %% 2 == 0) fs_used <- fs_used + 1
        }

        tmp_in  <- tempfile(fileext = ".tif")
        tmp_fill<- tempfile(fileext = ".tif")
        terra::writeRaster(r_aln, filename = tmp_in, overwrite = TRUE,
                           gdal = c("COMPRESS=NONE","TILED=YES","BLOCKXSIZE=256","BLOCKYSIZE=256"))
        ok <- TRUE
        say(sprintf("[%s] Whitebox fill: filter=%d, weight=%.2f", lyr_name, fs_used, idw_weight))
        tryCatch({
          whitebox::wbt_fill_missing_data(
            input    = tmp_in,
            output   = tmp_fill,
            filter   = as.integer(fs_used),
            weight   = idw_weight,
            no_edges = FALSE,
            verbose_mode = FALSE
          )
        }, error = function(e){
          warning(sprintf("[%s] Whitebox fill failed: %s (keeping unfilled)", lyr_name, conditionMessage(e)))
          ok <<- FALSE
        })
        if (ok && file.exists(tmp_fill)) {
          r_aln <- terra::rast(tmp_fill)
          terra::crs(r_aln) <- tmpl_wkt
          gap_filled_vec[i] <- TRUE
        }
        unlink(c(tmp_in, tmp_fill), force = TRUE)
      }
    }

    # --- plots (use aligned raster) ---
    if ((isTRUE(plot_result) || isTRUE(plot_gaps)) && interactive()) {
      oldpar <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(oldpar), add = TRUE)

      two_panels <- isTRUE(plot_result) && isTRUE(plot_gaps)
      if (two_panels) graphics::par(mfrow = c(1, 2))

      if (isTRUE(plot_result)) {
        terra::plot(r_aln, main = paste0("landscape_function: ", lyr_name))
      }

      if (isTRUE(plot_gaps)) {
        if (!is.na(gap_count) && gap_count > 0) {
          # real gap map (1 = gap inside template)
          gap_r <- terra::ifel(is.na(r_aln) & !is.na(tmpl), 1, NA)
          terra::plot(gap_r, main = "Gaps (1 = NA inside template)", col = c(NA, "red"), legend = FALSE)
        } else {
          # show a clear "no gaps" panel
          base_no_gaps <- terra::ifel(!is.na(tmpl), 0, NA)
          terra::plot(base_no_gaps, main = "Gaps: none", col = "white", legend = FALSE)
          box()
        }
      }
    }

    # --- final write (atomic)
    tmp_final <- sprintf("%s._tmp.%s", stem, ext)
    if (file.exists(tmp_final)) try(unlink(tmp_final), silent = TRUE)
    names(r_aln) <- lyr_name
    terra::writeRaster(r_aln, filename = tmp_final, overwrite = TRUE,
                       gdal = gdal_final, datatype = dt_final, NAflag = na_final)
    if (file.exists(final_path)) try(unlink(final_path), silent = TRUE)
    file.rename(tmp_final, final_path)
    outputs[[lyr_name]] <- final_path

    # record stats
    gap_count_vec[i] <- gap_count
    max_gap_vec[i]   <- max_gap
    fs_used_vec[i]   <- fs_used

    # cleanup the merged mosaic file (keep only the aligned final)
    if (file.exists(tmp_merge)) try(unlink(tmp_merge), silent = TRUE)
  }


  elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  df <- data.frame(
    layer_name               = out_layername,
    output_path              = unname(unlist(outputs)),
    tiles_written            = layer_counts$tiles_written,
    tiles_skipped_existing   = layer_counts$tiles_skipped_existing,
    tiles_skipped_empty_crop = layer_counts$tiles_skipped_empty_crop,
    tiles_skipped_empty_join = layer_counts$tiles_skipped_empty_join,
    merge_skipped            = as.logical(unname(unlist(merge_skipped))),
    gap_count                = gap_count_vec,
    max_gap_distance         = max_gap_vec,
    filter_size_cells_used   = fs_used_vec,
    gap_filled               = gap_filled_vec,
    n_tiles                  = length(tiles_sf),
    n_zones                  = nrow(z_in),
    n_layers                 = nlyr,
    elapsed_sec              = elapsed_sec,
    stringsAsFactors = FALSE
  )

  attr(df, "tile_dir") <- if (keep_tiles) tile_root else NULL
  attr(df, "run_params") <- list(
    what = what, lm_args = lm_args, rasterize_engine = rasterize_engine,
    n_workers = n_workers, os_type = os_type, buffer_m = buffer_m,
    gdal_opts = gdal_final, write_datatype = safe_datatype(write_datatype),
    NAflag = if (is.null(NAflag) || is.na(NAflag)) infer_default_naflag(safe_datatype(write_datatype)) else NAflag,
    skip_existing = skip_existing,
    report_gaps = report_gaps, report_gap_size = report_gap_size,
    fill_gaps = fill_gaps, idw_weight = idw_weight, filter_size_cells = filter_size_cells,
    terra_memfrac = terra_memfrac, terra_tempdir = terra_tempdir, terra_todisk = terra_todisk
  )

  df
}
