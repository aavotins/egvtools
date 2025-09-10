# Package-private cache env (not exported)
.egvtools_cache <- new.env(parent = emptyenv())

# Internal helpers ----
.egv_cache_get <- function(name) {
  if (exists(name, envir = .egvtools_cache, inherits = FALSE)) {
    get(name, envir = .egvtools_cache, inherits = FALSE)
  } else {
    NULL
  }
}

.egv_cache_set <- function(name, value) {
  assign(name, value, envir = .egvtools_cache)
  invisible(value)
}

# Optional utilities (call manually if you want)
.egv_cache_drop  <- function(name) if (exists(name, envir = .egvtools_cache, inherits = FALSE))
  rm(list = name, envir = .egvtools_cache)
.egv_cache_clear <- function() rm(list = ls(envir = .egvtools_cache), envir = .egvtools_cache)
