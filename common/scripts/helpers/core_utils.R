# common/scripts/helpers/core_utils.R

#' Get maximum and user-requested number of cores
#'
#' This detects the max available cores (from SLURM or local machine),
#' and allows the user to specify how many cores they want to use.
#'
#' @param requested_cores Optional user-specified number of cores to use. If NULL, use all available.
#'
#' @return A list with `max_cores` and `num_workers`
#' @export
get_core_config <- function(requested_cores = NULL) {
  
  mc_cores_env <- Sys.getenv("MC_CORES", unset = NA)
  
  max_cores <- if (!is.na(mc_cores_env)) {
    as.integer(mc_cores_env)
  } else {
    parallel::detectCores()
  }
  
  if (!is.null(requested_cores)) {
    if (is.na(requested_cores) || requested_cores <= 0) {
      stop("Invalid requested_cores value. Must be a positive integer.")
    }
    if (requested_cores > max_cores) {
      warning("Requested more cores than available; using max instead")
    }
  }
  
  num_workers <- min(requested_cores %||% max_cores, max_cores)
  
  list(
    max_cores = max_cores,
    num_workers = num_workers
  )
}
