#!/usr/bin/env Rscript

# Suppress messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

message("✅ Serial R script started")

# Minimal computation
x <- 1:10
y <- x^2 + rnorm(10)
result <- data.frame(x = x, y = y)

message("Result:")
print(result)

# Write output file
out_file <- here("test_serial_output.rds")
saveRDS(result, out_file)
message("✅ Output saved to: ", out_file)

message("Serial R script completed successfully")
