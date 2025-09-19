#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

# -------------------------------
# Minimal multicore test
# -------------------------------
message("Starting multicore test...")

plan(multicore, workers = 4)   # adjust to match your SLURM cpus
res <- future_lapply(1:10, function(x) {
  Sys.getpid()   # return PID of worker process
})

plan(sequential) # clean up plan

message("Results:")
print(res)

message("Unique worker PIDs used: ", length(unique(unlist(res))))
