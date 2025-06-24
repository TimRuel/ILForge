# common/scripts/helpers/prompt_utils.R

should_proceed_or_abort <- function(prompt, force = FALSE, non_interactive = FALSE) {
  if (force) return(TRUE)
  if (non_interactive) return(FALSE)
  
  cat(paste0(prompt, " [y/N]: "))
  flush.console()  # ensure prompt shows immediately
  
  response <- tolower(trimws(readLines(con = "stdin", n = 1)))
  response %in% c("y", "yes")
}

