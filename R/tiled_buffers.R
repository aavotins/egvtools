#' Create buffered tiles from point layers
#'
#' @description
#' Buffers point parquet files and writes **per-tile** outputs named:
#'   - Fixed radius: `<base>_r<radius>_<tileid>.parquet`
#'   - Per-feature radii (`radius_field`): `<base>_varradius_<tileid>.parquet`
#'
#' Where:
#' - **base** = prefix of input filename before the first underscore (e.g., "pts100" from `pts100_sauzeme.parquet`)
#' - **tileid** = unique value of `split_field` (default "tks50km") inside the input file
#'
#' The `buffer_mode` determines how radii are assigned:
#' - **"dense"**: Buffers the best-matching `pts100*.parquet` (prefers `pts100_sauzeme.parquet`)
#'   for each tile by `radii_dense` (default: 500, 1250, 3000, 10000 m).
#' - **"sparse"**: Uses a file to radius mapping. Default mapping:
#'     * `pts100_sauzeme.parquet` to c(500, 1250)
#'     * `pts300_sauzeme.parquet` to 3000
#'     * `pts1000_sauzeme.parquet` to 10000
#'   You can override this via `mapping_sparse` (named list or data.frame with columns `file`, `radius`).
#' - **"specified"**: You provide `points_path` and either
#'   `buffer_radius` (uniform, one or more fixed radii) **or** `radius_field`
#'   (numeric column with per-feature radii in meters). Outputs are still split by `split_field`.
#'
#' @param in_dir Character. Directory containing input points parquet files
#'   (default "./Templates/TemplateGridPoints"). Used in "dense" and "sparse" modes.
#' @param out_dir Character. Output directory for buffered tiles
#'   (default "./Templates/TemplateGridPoints/lapas").
#' @param buffer_mode Character. One of "dense", "sparse", "specified". Default "dense".
#' @param radii_dense Numeric vector of radii (m) used when `buffer_mode = "dense"`.
#'   Default c(500, 1250, 3000, 10000).
#' @param mapping_sparse Named list **or** data.frame describing file to radii for
#'   `buffer_mode = "sparse"`. Default:
#'   `list("pts100_sauzeme.parquet" = c(500, 1250),
#'         "pts300_sauzeme.parquet" = 3000,
#'         "pts1000_sauzeme.parquet" = 10000)`.
#'   If a data.frame is supplied, it must have columns `file` and `radius`.
#' @param points_path Character. Path to a single parquet file for
#'   `buffer_mode = "specified"`.
#' @param buffer_radius Numeric vector. Used in "specified" when you want fixed radii.
#'   Ignored if `radius_field` is provided.
#' @param radius_field Character or NULL. Column name in `points_path` that provides
#'   per-feature radii (meters) for "specified". If given, `buffer_radius` is ignored.
#' @param split_field Character. Field in the point data that defines tiles. Default "tks50km".
#' @param n_workers Integer. Parallel workers. Default `max(1L, parallel::detectCores())`.
#' @param os_type Optional character to force backend plan:
#'   "windows", "mac", "darwin", "linux", "slurm". Default NULL becomes auto-detect.
#' @param future_max_mem_gb Numeric. Max size of exported globals per worker (GiB).
#'   Sets `options(future.globals.maxSize = future_max_mem_gb * 1024^3)`. Default 4.
#' @param overwrite Logical. Overwrite existing outputs? Default FALSE.
#' @param quiet Logical. Suppress messages? Default FALSE.
#'
#' @return Invisibly returns a data.frame with columns:
#'   `input`, `tileid`, `mode`, `radius_m` (NA if `radius_field` is used),
#'   `radius_field`, `out_file`, `wrote`.
#'
#' @details
#' - Uses \pkg{sfarrow} for parquet I/O and \pkg{sf} for buffering.
#' - Jobs are created **per tile** (unique `split_field` value) and (when applicable) per radius.
#' - Workers open data from file paths (keeps RAM low).
#' - Files are written atomically (temp file then move).
#' - **Safety**: jobs are de-duplicated by content **and** by predicted output path, so concurrent
#'   workers never write to the same file.
#'
#' @seealso [tile_vector_grid()]
#' @source Zenodo grids/points example: https://zenodo.org/records/14277114
#'
#' @examples
#' \dontrun{
#' tiled_buffers(buffer_mode = "sparse", split_field = "tks50km")
#' }
#'
#' @importFrom fs dir_exists dir_create dir_ls file_exists path_file is_file file_delete file_move
#' @importFrom sfarrow st_read_parquet st_write_parquet
#' @importFrom sf st_buffer st_is_empty st_make_valid st_geometry st_geometry<- sf_use_s2
#' @importFrom furrr future_map
#' @importFrom future plan multisession multicore cluster sequential
#' @export
tiled_buffers <- function(
    in_dir         = "./Templates/TemplateGridPoints",
    out_dir        = "./Templates/TemplateGridPoints/lapas",
    buffer_mode    = "dense",
    radii_dense    = c(500, 1250, 3000, 10000),
    mapping_sparse = list(
      "pts100_sauzeme.parquet"  = c(500, 1250),
      "pts300_sauzeme.parquet"  = 3000,
      "pts1000_sauzeme.parquet" = 10000
    ),
    points_path    = NULL,
    buffer_radius  = NULL,
    radius_field   = NULL,
    split_field    = "tks50km",
    n_workers      = max(1L, parallel::detectCores()),
    os_type        = NULL,
    future_max_mem_gb = 4,
    overwrite      = FALSE,
    quiet          = FALSE
) {

  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("sf",      "reading points & buffering")
  .need_pkg("sfarrow", "GeoParquet IO")
  .need_pkg("fs",      "filesystem ops")
  # Parallel (only if used):
  if (n_workers > 1) { .need_pkg("future","parallel"); .need_pkg("furrr","parallel mapping") }

  if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)

  # ---- sink safety: snapshot & restore on exit (protects against stuck sinks) ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)

  say <- function(...) if (!quiet) cat(..., "\n")



  # ---- preserve caller's plan and options; restore on exit ----
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  old_opts <- options(future.globals.maxSize = as.numeric(future_max_mem_gb) * (1024^3))
  on.exit(options(old_opts), add = TRUE)

  # ---- choose parallel backend ----
  plan_auto <- function(os_detected) {
    if (!is.null(os_detected) && os_detected == "slurm") {
      cl <- tryCatch(parallel::makeCluster(n_workers), error = function(e) NULL)
      if (!is.null(cl)) {
        future::plan(future::cluster, workers = cl)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE, after = TRUE)
        if (!quiet) say("Using SLURM cluster backend (", n_workers, " workers).")
        return(invisible())
      }
    }
    sys <- tolower(Sys.info()[["sysname"]])
    if (!is.null(os_detected)) sys <- os_detected
    if (sys %in% c("windows","mac","darwin")) {
      future::plan(future::multisession, workers = n_workers)
      if (!quiet) say("Using multisession (", n_workers, " workers).")
    } else if (sys == "linux") {
      future::plan(future::multicore, workers = n_workers)
      if (!quiet) say("Using multicore on Linux (", n_workers, " workers).")
    } else {
      future::plan(future::sequential)
      if (!quiet) say("Falling back to sequential.")
    }
  }
  detect_os <- function() {
    if (!is.null(os_type)) return(tolower(os_type))
    if (nzchar(Sys.getenv("SLURM_JOB_ID"))) return("slurm")
    tolower(Sys.info()[["sysname"]])
  }
  plan_auto(detect_os())

  # ---- helpers (exported to workers via globals=TRUE) ----
  get_base <- function(fname) {
    parts <- strsplit(fname, "_", fixed = TRUE)[[1]]
    parts[1]
  }

  # Read only the split_field using sfarrow (columns=); fallback to full read
  list_tiles <- function(pq_file, split_field) {
    x <- tryCatch(
      sfarrow::st_read_parquet(pq_file, columns = split_field),
      error = function(e) {
        y <- sfarrow::st_read_parquet(pq_file)
        if (!split_field %in% names(y)) stop("split_field '", split_field, "' not found in: ", pq_file)
        y[split_field]
      }
    )
    if (!nrow(x)) return(character(0))
    unique(x[[split_field]])
  }

  read_points <- function(p) {
    sf::sf_use_s2(FALSE)
    x <- sfarrow::st_read_parquet(p)
    if (nrow(x)) {
      if (!all(sf::st_is_empty(x) == FALSE)) x <- x[!sf::st_is_empty(x), , drop = FALSE]
      x <- sf::st_make_valid(x)
    }
    x
  }

  write_parquet_atomic <- function(obj, out_file) {
    tmp <- paste0(out_file, ".tmp")
    sfarrow::st_write_parquet(obj, tmp)
    if (fs::file_exists(out_file)) fs::file_delete(out_file)
    fs::file_move(tmp, out_file)
    out_file
  }

  predict_out_r <- function(input, tileid, r, out_dir) {
    base <- get_base(fs::path_file(input))
    file.path(out_dir, sprintf("%s_r%s_%s.parquet", base, as.integer(r), as.character(tileid)))
  }
  predict_out_field <- function(input, tileid, out_dir) {
    base <- get_base(fs::path_file(input))
    file.path(out_dir, sprintf("%s_varradius_%s.parquet", base, as.character(tileid)))
  }

  # Buffer a single (file, tileid) with a constant radius
  do_buffer_const <- function(in_file, tileid, r, split_field, out_dir, overwrite) {
    pts <- read_points(in_file)
    if (!nrow(pts)) return(NULL)
    if (!split_field %in% names(pts)) stop("split_field '", split_field, "' not found in: ", in_file)
    pts_t <- pts[pts[[split_field]] == tileid, , drop = FALSE]
    if (!nrow(pts_t)) return(NULL)

    out_file <- predict_out_r(in_file, tileid, r, out_dir)
    if (fs::file_exists(out_file) && !overwrite) {
      return(data.frame(input = in_file, tileid = tileid, mode = "const",
                        radius_m = r, radius_field = NA_character_,
                        out_file = out_file, wrote = FALSE))
    }

    buf <- sf::st_buffer(pts_t, dist = r)

    write_parquet_atomic(buf, out_file)
    data.frame(input = in_file, tileid = tileid, mode = "const",
               radius_m = r, radius_field = NA_character_,
               out_file = out_file, wrote = TRUE)
  }

  # Buffer a single (file, tileid) using per-feature radius_field
  do_buffer_field <- function(in_file, tileid, field, split_field, out_dir, overwrite) {
    pts <- read_points(in_file)
    if (!nrow(pts)) return(NULL)
    if (!split_field %in% names(pts)) stop("split_field '", split_field, "' not found in: ", in_file)
    if (!field %in% names(pts)) stop("radius_field '", field, "' not found in: ", in_file)

    pts_t <- pts[pts[[split_field]] == tileid, , drop = FALSE]
    if (!nrow(pts_t)) return(NULL)

    dists <- pts_t[[field]]
    if (!is.numeric(dists)) stop("radius_field must be numeric (meters).")

    out_file <- predict_out_field(in_file, tileid, out_dir)
    if (fs::file_exists(out_file) && !overwrite) {
      return(data.frame(input = in_file, tileid = tileid, mode = "field",
                        radius_m = NA_real_, radius_field = field,
                        out_file = out_file, wrote = FALSE))
    }

    buf <- sf::st_buffer(pts_t, dist = dists)

    write_parquet_atomic(buf, out_file)
    data.frame(input = in_file, tileid = tileid, mode = "field",
               radius_m = NA_real_, radius_field = field,
               out_file = out_file, wrote = TRUE)
  }

  # ---- build jobs per mode (file by tile by radii/field) ----
  if (buffer_mode == "dense") {
    if (!fs::dir_exists(in_dir)) stop("in_dir not found: ", in_dir)
    all_files <- fs::dir_ls(in_dir, recurse = FALSE, type = "file")
    cand <- all_files[grepl("(?i)^pts100.*\\.parquet$", fs::path_file(all_files))]
    if (!length(cand)) cand <- all_files[grepl("(?i)^pts.*\\.parquet$", fs::path_file(all_files))]
    if (!length(cand)) stop("No pts*.parquet found in: ", in_dir)
    target <- cand[grepl("(?i)^pts100_sauzeme\\.parquet$", fs::path_file(cand))]
    if (!length(target)) target <- cand[1]

    if (!length(radii_dense)) stop("radii_dense must have at least one element.")
    tiles <- list_tiles(target, split_field)
    if (!length(tiles)) stop("No tiles found in split_field '", split_field, "' of: ", target)

    jobs <- expand.grid(input = target, tileid = tiles, r = radii_dense, stringsAsFactors = FALSE)
    mode_tag <- "dense"

  } else if (buffer_mode == "sparse") {
    if (!fs::dir_exists(in_dir)) stop("in_dir not found: ", in_dir)
    all_files <- fs::dir_ls(in_dir, recurse = FALSE, type = "file")
    files <- all_files[grepl("(?i)^pts.*\\.parquet$", fs::path_file(all_files))]
    if (!length(files)) stop("No pts*.parquet found in: ", in_dir)

    # normalize mapping to data.frame(file, radius)
    if (is.list(mapping_sparse) && !is.data.frame(mapping_sparse)) {
      map_df <- do.call(rbind, lapply(names(mapping_sparse), function(nm) {
        rad <- mapping_sparse[[nm]]
        if (length(rad) == 0L) return(NULL)
        data.frame(file = nm, radius = as.numeric(rad), stringsAsFactors = FALSE)
      }))
    } else if (is.data.frame(mapping_sparse)) {
      stopifnot(all(c("file","radius") %in% names(mapping_sparse)))
      map_df <- data.frame(file = as.character(mapping_sparse$file),
                           radius = as.numeric(mapping_sparse$radius),
                           stringsAsFactors = FALSE)
    } else {
      stop("mapping_sparse must be a named list or a data.frame with columns 'file' and 'radius'.")
    }
    if (!nrow(map_df)) stop("mapping_sparse resolved to zero rows.")

    # resolve mapping file names to actual paths in in_dir (case-insensitive)
    resolve_file <- function(pat) {
      hits <- files[tolower(fs::path_file(files)) == tolower(pat)]
      if (length(hits)) return(hits[1])
      NA_character_
    }
    map_df$path <- vapply(map_df$file, resolve_file, character(1))
    map_df <- map_df[!is.na(map_df$path), , drop = FALSE]
    if (!nrow(map_df)) stop("No mapping entries matched files in: ", in_dir)

    # enumerate tiles per file, then radii
    job_list <- list()
    for (k in seq_len(nrow(map_df))) {
      f   <- map_df$path[k]
      rad <- map_df$radius[k]
      tiles <- list_tiles(f, split_field)
      if (!length(tiles)) next
      job_list[[length(job_list) + 1L]] <-
        expand.grid(input = f, tileid = tiles, r = rad, stringsAsFactors = FALSE)
    }
    if (!length(job_list)) stop("No tiles found for mapping across files.")
    jobs <- do.call(rbind, job_list)
    mode_tag <- "sparse"

  } else if (buffer_mode == "specified") {
    if (is.null(points_path) || !fs::is_file(points_path)) {
      stop("When buffer_mode='specified', provide a valid points_path.")
    }
    tiles <- list_tiles(points_path, split_field)
    if (!length(tiles)) stop("No tiles found in split_field '", split_field, "' of: ", points_path)

    if (!is.null(radius_field)) {
      jobs <- data.frame(input = points_path, tileid = tiles, field = radius_field, stringsAsFactors = FALSE)
      mode_tag <- "specified_field"
    } else {
      if (is.null(buffer_radius) || !length(buffer_radius)) {
        stop("For 'specified' without radius_field, provide buffer_radius (numeric).")
      }
      jobs <- expand.grid(input = points_path, tileid = tiles, r = as.numeric(buffer_radius), stringsAsFactors = FALSE)
      mode_tag <- "specified_const"
    }

  } else {
    stop("buffer_mode must be one of: 'dense','sparse','specified'.")
  }

  # ---- SAFETY: dedup jobs (by content and by predicted out path) ----
  jobs <- unique(jobs)
  jobs$.out_file <- if ("r" %in% names(jobs)) {
    mapply(predict_out_r, jobs$input, jobs$tileid, jobs$r, MoreArgs = list(out_dir = out_dir), USE.NAMES = FALSE)
  } else {
    mapply(predict_out_field, jobs$input, jobs$tileid, MoreArgs = list(out_dir = out_dir), USE.NAMES = FALSE)
  }
  jobs <- jobs[!duplicated(jobs$.out_file), , drop = FALSE]

  # ---- run jobs (export helpers via globals=TRUE so workers can see them) ----
  opts <- furrr::furrr_options(globals = TRUE)
  if (mode_tag %in% c("dense","sparse","specified_const")) {
    out <- furrr::future_map(
      seq_len(nrow(jobs)),
      function(i) do_buffer_const(jobs$input[i], jobs$tileid[i], jobs$r[i],
                                  split_field = split_field, out_dir = out_dir, overwrite = overwrite),
      .options = opts
    )
  } else { # specified_field
    out <- furrr::future_map(
      seq_len(nrow(jobs)),
      function(i) do_buffer_field(jobs$input[i], jobs$tileid[i], jobs$field[i],
                                  split_field = split_field, out_dir = out_dir, overwrite = overwrite),
      .options = opts
    )
  }
  res <- do.call(rbind, out)

  if (!quiet) {
    wrote <- sum(res$wrote, na.rm = TRUE)
    say("tiled_buffers complete. Wrote ", wrote, " / ", nrow(res), " files at ", out_dir)
  }
  invisible(res)
}
