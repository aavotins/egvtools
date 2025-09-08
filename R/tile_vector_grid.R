#' Tile a vector grid into parquet tiles
#'
#' @description
#' Reads a grid parquet (e.g., `tikls100_sauzeme.parquet`) and writes per-tile parquet
#' files into `out_dir`. Output filenames are formed as:
#' **`[base]_[tilename].parquet`**, where:
#' - **base** = the substring of the input filename before the first underscore
#'   (e.g., `"tikls100"` from `tikls100_sauzeme.parquet`);
#' - **tilename** = the value from `tile_field` (e.g., `2434`), or an auto-generated
#'   chunk id (e.g., `tile_00001`) if `tile_field` is `NULL` or not present.
#'
#' Tiling can be controlled either by an existing `tile_field` (preferred) or by
#' automatic chunking using `chunk_size`.
#'
#' @param grid_path Character. Path to the source grid parquet
#'   (default `"./Templates/TemplateGrids/tikls100_sauzeme.parquet"`).
#' @param out_dir Character. Where to write tiles (default
#'   `"./Templates/TemplateGrids/lapas"`).
#' @param tile_field Character or `NULL`. If provided and exists in data (e.g., `"lapa"`),
#'   tiles are split by unique values of this field. If `NULL` or not found, uses `chunk_size`.
#'   Default `"lapa"`.
#' @param chunk_size Integer. Number of rows per tile when `tile_field` is `NULL` or missing.
#'   Default `50000L`.
#' @param overwrite Logical. Overwrite existing tile files? Default `FALSE`.
#' @param quiet Logical. Suppress messages? Default `FALSE`.
#'
#' @return Invisibly returns a `data.frame` with columns: `tile_id`, `n_rows`, `path`, `wrote`.
#'
#' @details
#' Uses {sfarrow} for Parquet I/O. CRS and geometry are preserved. Designed for stable,
#' idempotent outputs, deterministic filenames, and clean skipping of existing tiles.
#'
#' @seealso [download_vector_templates()], [tiled_buffers()]
#' @source Zenodo (grid example): https://zenodo.org/records/14277114
#'
#' @examples
#' \dontrun{
#' # Split by an existing field (e.g., "lapa")
#' tile_vector_grid(
#'   grid_path = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
#'   out_dir   = "./Templates/TemplateGrids/lapas",
#'   tile_field = "lapa"
#' )
#'
#' # Fallback to chunking if tile_field is NULL or missing
#' tile_vector_grid(
#'   grid_path  = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
#'   out_dir    = "./Templates/TemplateGrids/lapas",
#'   tile_field = NULL,
#'   chunk_size = 25000
#' )
#' }
#'
#' @importFrom sfarrow st_read_parquet st_write_parquet
#' @importFrom sf st_crs st_geometry_type
#' @importFrom fs dir_exists dir_create file_exists file_delete file_move path_file
#' @export
tile_vector_grid <- function(
    grid_path  = "./Templates/TemplateGrids/tikls100_sauzeme.parquet",
    out_dir    = "./Templates/TemplateGrids/lapas",
    tile_field = "lapa",
    chunk_size = 50000L,
    overwrite  = FALSE,
    quiet      = FALSE
) {
  
  
  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("sf",      "grid read/split")
  .need_pkg("sfarrow", "GeoParquet IO")
  .need_pkg("fs",      "filesystem ops")
  
  
  if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)
  if (!fs::file_exists(grid_path)) stop("grid_path not found: ", grid_path)
  
  # ---- sink safety: snapshot & restore on exit (protects against stuck sinks) ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)
  
  say <- function(...) if (!quiet) cat(..., "\n")
  
  
  # derive "base" = part of input filename before first "_"
  base <- {
    nm <- fs::path_file(grid_path)
    parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
    parts[1]
  }
  
  sf::sf_use_s2(FALSE)
  g <- sfarrow::st_read_parquet(grid_path)
  if (nrow(g) == 0L) stop("Grid is empty: ", grid_path)
  
  # Determine tiling strategy
  by_field <- !is.null(tile_field) && tile_field %in% names(g)
  
  if (by_field) {
    # Split by field (preferred)
    split_idx <- split(seq_len(nrow(g)), g[[tile_field]])
    tiles <- lapply(names(split_idx), function(k) {
      idx <- split_idx[[k]]
      list(tile_id = as.character(k), rows = idx)
    })
    if (!quiet) say("Tiling by '", tile_field, "' with ", length(tiles), " tiles.")
  } else {
    # Fallback to chunking by row blocks
    n <- nrow(g)
    ids <- ceiling(seq_len(n) / as.integer(chunk_size))
    split_idx <- split(seq_len(n), ids)
    tiles <- lapply(seq_along(split_idx), function(i) {
      idx <- split_idx[[i]]
      list(tile_id = sprintf("tile_%05d", i), rows = idx)
    })
    if (!quiet) say("Tiling by chunks of ~", chunk_size, " rows (", length(tiles), " tiles).")
  }
  
  res <- vector("list", length(tiles))
  for (i in seq_along(tiles)) {
    tile_id <- tiles[[i]]$tile_id
    idx     <- tiles[[i]]$rows
    if (!length(idx)) {
      res[[i]] <- data.frame(tile_id = tile_id, n_rows = 0L, path = NA_character_, wrote = FALSE)
      next
    }
    
    tile_sf <- g[idx, , drop = FALSE]
    out_file <- file.path(out_dir, paste0(base, "_", tile_id, ".parquet"))
    
    if (fs::file_exists(out_file) && !overwrite) {
      if (!quiet) say("Skip existing: ", fs::path_file(out_file))
      res[[i]] <- data.frame(tile_id = tile_id, n_rows = nrow(tile_sf), path = out_file, wrote = FALSE)
      next
    }
    
    tmp <- paste0(out_file, ".tmp")
    sfarrow::st_write_parquet(tile_sf, tmp)
    if (fs::file_exists(out_file)) fs::file_delete(out_file)
    fs::file_move(tmp, out_file)
    
    if (!quiet) say("Wrote: ", fs::path_file(out_file), " (", nrow(tile_sf), " rows)")
    res[[i]] <- data.frame(tile_id = tile_id, n_rows = nrow(tile_sf), path = out_file, wrote = TRUE)
  }
  
  out <- do.call(rbind, res)
  invisible(out)
}
