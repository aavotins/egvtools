#' Download and unpack raster templates
#'
#' @description
#' Downloads a ZIP from Zenodo (or a user URL), unpacks it into a destination
#' directory, and removes the archive after a successful extraction. Writes to a
#' tempfile first and then moves to final destination for atomicity.
#'
#' @param url Character. Archive URL to download.
#'   Default: Zenodo (EGV raster templates) v2.0.0:
#'   "https://zenodo.org/api/records/14497070/files-archive".
#' @param out_dir Character. Destination directory for the unzipped rasters.
#'   Default: "./Templates/TemplateRasters".
#' @param overwrite Logical. Re-download and overwrite existing files/dirs?
#'   Default: FALSE.
#' @param quiet Logical. Suppress progress messages. Default: FALSE.
#'
#' @return Invisibly returns a named list with `out_dir`, `files_written`.
#'
#' @seealso [download_vector_templates()], [tile_vector_grid()]
#' @source Zenodo: https://doi.org/10.5281/zenodo.14497070
#'
#' @examples
#' \dontrun{
#' download_raster_templates()
#' }
#'
#' @importFrom curl curl_download
#' @importFrom utils unzip
#' @importFrom fs dir_exists dir_create file_exists file_delete
#' @export
download_raster_templates <- function(
    url     = "https://zenodo.org/api/records/14497070/files-archive",
    out_dir = "./Templates/TemplateRasters",
    overwrite = FALSE,
    quiet = FALSE
) {
  # deps
  .need_pkg <- function(p, why) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for %s. Please install it.", p, why), call. = FALSE)
    }
  }
  .need_pkg("fs",    "dir creation & file moves")
  .need_pkg("utils", "download.file & unzip")  # utils is base-recommended; this keeps the intent clear
  
  
  
  if (!fs::dir_exists(out_dir)) fs::dir_create(out_dir, recurse = TRUE)
  
  # ---- sink safety: snapshot & restore on exit (protects against stuck sinks) ----
  orig_out <- sink.number()
  orig_msg <- sink.number(type = "message")
  on.exit({
    while (sink.number(type = "message") > orig_msg) sink(type = "message")
    while (sink.number() > orig_out) sink()
  }, add = TRUE)
  
  say <- function(...) if (!quiet) cat(..., "\n")
  
  
  
  # Simple idempotency check: if dir non-empty and !overwrite, skip
  if (!overwrite) {
    existing <- list.files(out_dir, all.files = FALSE, recursive = TRUE, no.. = TRUE)
    if (length(existing) > 0L) {
      if (!quiet) say("Raster templates already present in '", out_dir, "'. Set overwrite=TRUE to re-download.")
      return(invisible(list(out_dir = out_dir, files_written = character(0))))
    }
  }
  
  tmp_zip <- tempfile(fileext = ".zip")
  on.exit(if (fs::file_exists(tmp_zip)) fs::file_delete(tmp_zip), add = TRUE, after = TRUE)
  
  # Robust download with up to 3 tries
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
  
  # Unzip to temp dir, then move
  tmp_unzip <- tempfile("unz_")
  fs::dir_create(tmp_unzip)
  utils::unzip(tmp_zip, exdir = tmp_unzip)
  
  # Move content into out_dir
  unz_files <- list.files(tmp_unzip, full.names = TRUE, all.files = TRUE, recursive = FALSE, no.. = TRUE)
  for (p in unz_files) {
    # If the archive has a top-level folder, move its children
    if (fs::is_dir(p)) {
      inner <- fs::dir_ls(p, recurse = TRUE, all = TRUE)
      fs::dir_create(out_dir, recurse = TRUE)
      fs::file_copy(inner, file.path(out_dir, fs::path_rel(inner, start = p)), overwrite = overwrite)
    } else {
      fs::file_copy(p, file.path(out_dir, fs::path_file(p)), overwrite = overwrite)
    }
  }
  
  written <- list.files(out_dir, full.names = TRUE, recursive = TRUE)
  if (!quiet) say("Downloaded and unpacked ", length(written), " files to: ", out_dir)
  invisible(list(out_dir = out_dir, files_written = written))
}
