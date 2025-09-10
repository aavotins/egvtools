.onLoad <- function(libname, pkgname) {
  op <- options(
    egvtools.future_plan = 'sequential',  # default parallel plan
    egvtools.progress    = TRUE           # whether to show progress bars
  )
  toset <- !(names(op) %in% names(options()))
  if (any(toset)) options(op[toset])
  invisible()
}

.onAttach <- function(libname, pkgname) {
  if (!interactive() || isFALSE(getOption("egvtools.whitebox_msg", TRUE))) return()

  nudge <- function(txt) packageStartupMessage(paste0("egvtools: ", txt))

  if (!requireNamespace("whitebox", quietly = TRUE)) {
    nudge("package 'whitebox' is not installed. Install with install.packages('whitebox').")
    return()
  }

  ok <- isTRUE(tryCatch(whitebox::wbt_init(), error = function(e) FALSE))
  if (!ok) {
    nudge("WhiteboxTools binary not found. Run whitebox::install_whitebox() or whitebox::wbt_install() to install or point whitebox::wbt_init() to existing installation. ")
    nudge("To silence this message: options(egvtools.whitebox_msg = FALSE)")
  }
}
