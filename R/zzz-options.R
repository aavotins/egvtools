.onLoad <- function(libname, pkgname) {
  op <- options(
    hiqbiodiv.future_plan = 'sequential',  # default parallel plan
    hiqbiodiv.progress    = TRUE           # whether to show progress bars
  )
  toset <- !(names(op) %in% names(options()))
  if (any(toset)) options(op[toset])
  invisible()
}
