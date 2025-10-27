#!/usr/bin/env Rscript
# 01-generate-population.R
# This script is run ONCE to create the shared population data
# inside the specific results folder.

print("Loading libraries")
library(dplyr)

# --- Read command line arguments ---
args <- commandArgs(trailingOnly = TRUE)

csv_path         <- if (length(args) >= 1) args[1] else "2022ACS"
pop_params_path  <- if (length(args) >= 2) args[2] else "params-pop-1"
miss_params_path <- if (length(args) >= 3) args[3] else "missingness_params-1"
n_samp           <- if (length(args) >= 4) as.integer(args[4]) else 2
L                <- if (length(args) >= 5) as.integer(args[5]) else 10
n_chains         <- if (length(args) >= 6) as.integer(args[6]) else 2
n_burnin         <- if (length(args) >= 7) as.integer(args[7]) else 500
n_iter           <- if (length(args) >= 8) as.integer(args[8]) else 2000

print("--- Population Generation Script ---")
print(paste("Using parameters to build folder path..."))

# --- Create the *exact* same folder name ---
folder_name <- sprintf("pop-%s-miss-%s-n_samp-%d-L-%d-data-%s-chains-%d-burnin-%d-iter-%d", 
                       pop_params_path, miss_params_path, n_samp, L, csv_path, n_chains, n_burnin, n_iter)

print(paste("Aggregating results from:", folder_name))

# Initialize empty lists
SAMPLES <- vector("list", n_samp)
MICE_RESULTS <- vector("list", n_samp)
JAGS_RESULTS <- vector("list", n_samp)

for (i in 1:n_samp) {
  print(paste("Loading sample", i))
  
  # Define file paths
  sample_file <- file.path(folder_name, sprintf("sample_result_%d.rds", i))
  mice_file   <- file.path(folder_name, sprintf("mice_result_%d.rds", i))
  jags_file   <- file.path(folder_name, sprintf("jags_result_%d.rds", i))
  
  # Check if files exist (in case a job failed)
  if (file.exists(sample_file) && file.exists(mice_file) && file.exists(jags_file)) {
    SAMPLES[[i]]      <- readRDS(sample_file)
    MICE_RESULTS[[i]] <- readRDS(mice_file)
    JAGS_RESULTS[[i]] <- readRDS(jags_file)
  } else {
    warning(paste("Missing .rds files for sample", i, "- check logs."))
  }
}

# Save the final aggregated lists
print("Saving aggregated SAMPLES list...")
samples_save_path <- file.path(folder_name, "samples_list.rds")
saveRDS(SAMPLES, file = samples_save_path)

print("Saving aggregated MICE_RESULTS list...")
mice_results_save_path <- file.path(folder_name, "mice_results_list.rds")
saveRDS(MICE_RESULTS, file = mice_results_save_path)

print("Saving aggregated JAGS_RESULTS list...")
jags_results_save_path <- file.path(folder_name, "jags_results_list.rds")
saveRDS(JAGS_RESULTS, file = jags_results_save_path)

print("All aggregation complete!")