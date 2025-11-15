#' Null-coalescing operator
#'
#' Returns the left-hand side if it is not `NULL`, otherwise returns the
#' right-hand side. Used internally for configuration defaults.
#'
#' @param x Primary value.
#' @param y Fallback value if `x` is `NULL`.
#'
#' @return `x` if non-`NULL`, otherwise `y`.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_seed_for_iter <- function(base_seed, iter_id) {
  
  return(as.integer(abs(Reduce(`+`, utf8ToInt(iter_id))) + base_seed))
}

get_arg_or_stop <- function(i, name) {
  if (length(args) >= i) return(args[[i]])
  stop(sprintf("[ERROR] A %s was not provided.", name))
}

check_path_exists <- function(path, label) {
  if (!file.exists(path)) {
    rel_path <- sub(".*(/?experiments/.*)", "\\1", path)
    stop(sprintf("[ERROR] %s does not exist at /%s", label, rel_path))
  }
}
