save_iter_metadata <- function(metadata, output_dir) {
  # Ensure metadata directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save as both CSV and RDS
  metadata_df <- as_tibble(lapply(metadata, \(x) if (is.null(x)) NA else x))
  write_csv(metadata_df, here(output_dir, "metadata.csv"))
  saveRDS(metadata, here(output_dir, "metadata.rds"))
}

