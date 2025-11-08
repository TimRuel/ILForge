# applications/negative_binomial/process_rates_weighted_sum/naive_rates/scripts/helpers/config_utils.R

# -------------------------------------------------------------------------
# Configuration Utilities
# -------------------------------------------------------------------------

#' Recursively coerce numeric-looking strings in a list to numeric
#'
#' This utility scans through a configuration list (e.g., one loaded from YAML)
#' and converts any scalar character elements that look like numbers
#' (e.g., `"1e-12"`, `"1000"`, `"0.05"`) into true numeric values using
#' \code{as.numeric()}. Nested lists are handled recursively.
#'
#' It is particularly useful when YAML or JSON configuration files
#' contain numeric values wrapped in quotes, which are parsed as strings
#' by \pkg{yaml::read_yaml()} or similar functions.
#'
#' @param cfg A (potentially nested) list, typically representing a configuration
#'   object read from a YAML or JSON file.
#'
#' @return The same list structure, but with numeric-like strings coerced
#'   to numeric values where possible.
#'
#' @examples
#' cfg <- list(
#'   step_size = "0.25",
#'   fine_window = "1.0",
#'   nested = list(p_cutoff = "1e-12", num_workers = "4")
#' )
#' sanitize_config(cfg)
#'
#' @export
sanitize_config <- function(cfg) {
  for (nm in names(cfg)) {
    val <- cfg[[nm]]
    # Detect scalar strings that look numeric (including scientific notation)
    if (is.character(val) && length(val) == 1L && grepl("^[0-9eE.+-]+$", val)) {
      suppressWarnings(num_val <- as.numeric(val))
      if (!is.na(num_val)) cfg[[nm]] <- num_val
    } else if (is.list(val)) {
      cfg[[nm]] <- sanitize_config(val)
    }
  }
  cfg
}
