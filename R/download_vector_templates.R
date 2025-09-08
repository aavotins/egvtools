#' Download and unpack vector templates (points/grids/gpkg)
#'
#' @description
#' Downloads vector templates archive (by default from Zenodo) and places files into:
#' - Parquet grids into `grid_dir`
#' - Parquet points into `points_dir`
#' - GPKG into `gpkg_dir`
#' The function auto-classifies files by filename patterns: "tikls*.parquet" (grids),
#' "pts*.parquet" (points), and "vector_grids.gpkg" (GPKG).
#'
#' @param url Character. Default:
#'   "https://zenodo.org/api/records/14277114/files-archive"
#' @param grid_dir Character. Default "./Templates/TemplateGrids"
#' @param points_dir Character. Default "./Templates/TemplateGridPoints"
#' @param gpkg_dir Character. Default "./Templates"
#' @param overwrite Logical. Overwrite existing files? Default FALSE.
#' @param quiet Logical. Suppress progress messages? Default FALSE.
#'
#' @return Invisibly returns a list of the three dirs.
#'
#' @seealso [download_raster_templates()], [tile_vector_grid()]
#' @source Zenodo: https://zenodo.org/records/14277114
#'
#' @examples
#' \dontrun{
#' download_vector_templates()
#' }
#'
#' @importFrom curl curl_download
#' @importFrom utils unzip
#' @importFrom fs dir_exists dir_create file_exists file_delete file_copy path_file is_dir dir_ls path_rel
#' @export
download_vector_templates <- function(
    url        = "https://zenodo.org/api/records/14277114/files-archive",
    grid_dir   = "./Templates/TemplateGrids",
    points_dir = "./Templates/TemplateGridPoints",
    gpkg_dir   = "./Templates",
    overwrite  = FALSE,
    quiet      = FALSE
) {
  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("fs",    "dir creation & file moves")
  .need_pkg("utils", "download.file & unzip")
  
  # ---- sink safety: snapshot & restore on exit (protects against stuck sinks) ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)
  
  say <- function(...) if (!quiet) cat(..., "\n")
  
  for (d in c(grid_dir, points_dir, gpkg_dir)) if (!fs::dir_exists(d)) fs::dir_create(d, recurse = TRUE)
  
  tmp_zip <- tempfile(fileext = ".zip")
  on.exit(if (fs::file_exists(tmp_zip)) fs::file_delete(tmp_zip), add = TRUE, after = TRUE)
  
  ok <- FALSE
  for (i in 1:3) {
    try({
      curl::curl_download(url, destfile = tmp_zip, mode = "wb", quiet = quiet)
      ok <- TRUE
      break
    }, silent = TRUE)
    if (!quiet) say("Download attempt ", i, " failed, retrying...")
  }
  if (!ok) stop("Failed to download archive from: ", url)
  
  tmp_unzip <- tempfile("unz_")
  fs::dir_create(tmp_unzip)
  utils::unzip(tmp_zip, exdir = tmp_unzip)
  
  # Classify & copy
  files <- list.files(tmp_unzip, full.names = TRUE, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  for (f in files) {
    if (fs::is_dir(f)) next
    nm <- fs::path_file(f)
    if (grepl("^tikls.*\\.parquet$", nm, ignore.case = TRUE) || grepl("^tks93_50km.*\\.parquet$", nm, ignore.case = TRUE)) {
      fs::file_copy(f, file.path(grid_dir, nm), overwrite = overwrite)
    } else if (grepl("^pts.*\\.parquet$", nm, ignore.case = TRUE)) {
      fs::file_copy(f, file.path(points_dir, nm), overwrite = overwrite)
    } else if (tolower(nm) == "vector_grids.gpkg") {
      fs::file_copy(f, file.path(gpkg_dir, nm), overwrite = overwrite)
    }
  }
  
  if (!quiet) say("Vector templates ready:\n  grid_dir: ", grid_dir, "\n  points_dir: ", points_dir, "\n  gpkg_dir: ", gpkg_dir)
  invisible(list(grid_dir = grid_dir, points_dir = points_dir, gpkg_dir = gpkg_dir))
}
