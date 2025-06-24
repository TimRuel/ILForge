append_run_log <- function(metadata, log_path) {
  # Convert to tibble
  new_row <- as_tibble(metadata)
  
  # Load existing log if it exists
  if (file.exists(log_path)) {
    existing_log <- read_csv(log_path, show_col_types = FALSE)
    
    # Convert date/time columns in new_row to match existing_log
    datetime_cols <- intersect(c("timestamp_start", "timestamp_end"), colnames(existing_log))
    for (col in datetime_cols) {
      if (col %in% names(new_row)) {
        new_row[[col]] <- as.POSIXct(new_row[[col]], tz = "UTC")
      }
    }
    
    # Bind and save
    full_log <- bind_rows(existing_log, new_row)
  } else {
    # If no existing log, initialize all date/time fields as POSIXct
    for (col in c("timestamp_start", "timestamp_end")) {
      if (col %in% names(new_row)) {
        new_row[[col]] <- as.POSIXct(new_row[[col]], tz = "UTC")
      }
    }
    full_log <- new_row
  }
  
  # Save back to CSV
  write_csv(full_log, log_path)
}

save_run_metadata <- function(metadata, output_dir) {
  # Ensure metadata directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save as both CSV and RDS
  metadata_df <- as_tibble(lapply(metadata, \(x) if (is.null(x)) NA else x))
  write_csv(metadata_df, here(output_dir, "metadata.csv"))
  saveRDS(metadata, here(output_dir, "metadata.rds"))
}

