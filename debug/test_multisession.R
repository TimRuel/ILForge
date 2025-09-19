#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(doFuture)
  library(foreach)
  library(future)
})

registerDoFuture()

args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args) > 0) as.integer(args[1]) else 1
message("Running replication ", run_id)
set.seed(1000 + run_id)

run_branch_side <- function(i) {
  Sys.sleep(runif(1, 0, 0.1))
  data.frame(i = i, result = rnorm(1))
}

compute_IL_branch <- function(j) {
  Sys.sleep(runif(1, 0, 0.05))
  rnorm(1)
}

# Outer parallel loop
N_outer <- 10
plan(multisession, workers = 4)
res_outer <- foreach(i = 1:N_outer, .combine = rbind) %dopar% {
  run_branch_side(i)
}
plan(sequential)

message("Outer loop results:")
print(res_outer)

# Inner parallel loop
N_inner <- 20
plan(multisession, workers = 4)
res_inner <- foreach(j = 1:N_inner, .combine = c) %dopar% {
  compute_IL_branch(j)
}
plan(sequential)

message("Inner loop results:")
print(res_inner)
