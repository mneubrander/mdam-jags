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

print(paste("Creating results folder:", folder_name))
dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)

# --- Data Generation ---
print("Loading source data...")
set.seed(456) # Use the fixed seed

df <- read.csv(paste0("data/", csv_path, ".csv"))
N_pop <- nrow(df)
param_list_pop <- readRDS(paste0("params/", pop_params_path, ".rds"))

print("Generating population data...")
pop_data <- data.frame(id = 1:N_pop)
pop_data$W <- 10 * df$PWGTP
pop_data$sampling_prob <- 1 / pop_data$W

with(as.list(param_list_pop), {
  pop_data$x1 <<- rbinom(N_pop, 1, plogis(a1 + b1w * pop_data$W))
  pop_data$x2 <<- rbinom(N_pop, 1, plogis(a2 + b21 * pop_data$x1))
  pop_data$x3 <<- rbinom(N_pop, 1, plogis(a3 + b31 * pop_data$x1 + b32 * pop_data$x2))
  pop_data$x4 <<- rbinom(N_pop, 1, plogis(a4 + b41 * pop_data$x1 + b42 * pop_data$x2 + b43 * pop_data$x3))
  pop_data$x5 <<- rnorm(N_pop, a5 + b51 * pop_data$x1 + b52 * pop_data$x2 + b53 * pop_data$x3 + b54 * pop_data$x4, s5)
  pop_data$x6 <<- rnorm(N_pop, a6 + b61 * pop_data$x1 + b62 * pop_data$x2 + b63 * pop_data$x3 +  b64 * pop_data$x4 + b65 * pop_data$x5, s6)
})

# --- Save to the specific results folder ---
pop_save_path <- file.path(folder_name, "population_data.rds")

print(paste("Saving population data to:", pop_save_path))
saveRDS(pop_data, file = pop_save_path)

print("Population data generated successfully.")