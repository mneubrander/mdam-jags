#!/usr/bin/env Rscript

########################################################################
##############################  SETUP ##################################

#### Read params from sbatch command line
# Command line arguments order:
# 1: csv_path (e.g., "2022ACS")
# 2: pop_params_path (e.g., "params-pop-1") 
# 3: miss_params_path (e.g., "missingness_params-1")
# 4: n_samp (e.g., 2)
# 5: L (e.g., 10)
# 6: n_chains (e.g., 2)
# 7: n_burnin (e.g., 500)
# 8: n_iter (e.g., 2000)

# Load required libraries
print("Loading libraries")
library(reshape2) # For melt
library(ggplot2)
library(tidyr) # For pivot_longer
library(dplyr) # For group_by and summarize
library(scales) # For percent_format
library(mice)
library(tidyverse)
library(parallel)
library(jagsUI)
print("Libraries loaded")

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables with defaults for interactive testing
csv_path         <- if (length(args) >= 1) args[1] else "2022ACS"
pop_params_path  <- if (length(args) >= 2) args[2] else "params-pop-1"
miss_params_path <- if (length(args) >= 3) args[3] else "missingness_params-1"
n_samp           <- if (length(args) >= 4) as.integer(args[4]) else 2
L                <- if (length(args) >= 5) as.integer(args[5]) else 10
n_chains         <- if (length(args) >= 6) as.integer(args[6]) else 2
n_burnin         <- if (length(args) >= 7) as.integer(args[7]) else 500
n_iter           <- if (length(args) >= 8) as.integer(args[8]) else 2000

# Print inputs for logging
print(paste("n_samp:", n_samp))
print(paste("L:", L))
print(paste("n_chains:", n_chains))
print(paste("n_burnin:", n_burnin))
print(paste("n_iter:", n_iter))
print(paste("CSV path:", csv_path))
print(paste("Pop params path:", pop_params_path))
print(paste("Missingness params path:", miss_params_path))

# get seed from command line when batch the code
set.seed(456)

df <- read.csv(paste0("data/",csv_path, ".csv"))
N_pop <- nrow(df)
param_list_pop <- readRDS(paste0("params/",pop_params_path,".rds"))
params_miss <- readRDS(paste0("params/", miss_params_path,".rds"))

# create folder for results
folder_name <- sprintf("pop-%s-miss-%s-n_samp-%d-L-%d-data-%s-chains-%d-burnin-%d-iter-%d", 
                       pop_params_path, miss_params_path, n_samp, L, csv_path, n_chains, n_burnin, n_iter)

print(paste("Creating results folder:", folder_name))
dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)

# Generate population data sequentially
pop_data <- data.frame(id = 1:N_pop)
pop_data$W <- 10*df$PWGTP
pop_data$sampling_prob <- 1/pop_data$W

with(as.list(param_list_pop), {
  pop_data$x1 <<- rbinom(N_pop, 1, plogis(a1 + b1w * pop_data$W))
  pop_data$x2 <<- rbinom(N_pop, 1, plogis(a2 + b21 * pop_data$x1))
  pop_data$x3 <<- rbinom(N_pop, 1, plogis(a3 + b31 * pop_data$x1 + b32 * pop_data$x2))
  pop_data$x4 <<- rbinom(N_pop, 1, plogis(a4 + b41 * pop_data$x1 + b42 * pop_data$x2 + b43 * pop_data$x3))
  pop_data$x5 <<- rnorm(N_pop, a5 + b51 * pop_data$x1 + b52 * pop_data$x2 + b53 * pop_data$x3 + b54 * pop_data$x4, s5)
  pop_data$x6 <<- rnorm(N_pop, a6 + b61 * pop_data$x1 + b62 * pop_data$x2 + b63 * pop_data$x3 +  b64 * pop_data$x4 + b65 * pop_data$x5, s6)
})

pop_save_path <- file.path(folder_name, "population_data.rds")
print(paste("Saving population data to:", pop_save_path))
saveRDS(pop_data, file = pop_save_path)

T1_known = sum(pop_data$x1)
T2_known = sum(pop_data$x2)
T3_known = sum(pop_data$x3)
T4_known = sum(pop_data$x4)
T5_known = sum(pop_data$x5)
T6_known = sum(pop_data$x6)

# Initialize function for JAGS will be defined inside process_sample

# Function to process a single sample
process_sample <- function(i) {
  
  print(paste("Processing sample", i, "on process", Sys.getpid()))
  
  ##############################################################
  ######################### GET SAMPLE #########################
  ##############################################################
  
  set.seed(i)
  sample_indices <- which(runif(N_pop) < pop_data$sampling_prob)
  sample_data <- pop_data[sample_indices, ]

  T1_samp = sum(sample_data$x1 * 1/sample_data$sampling_prob)
  T2_samp = sum(sample_data$x2 * 1/sample_data$sampling_prob)
  T3_samp = sum(sample_data$x3 * 1/sample_data$sampling_prob)
  T4_samp = sum(sample_data$x4 * 1/sample_data$sampling_prob)
  T5_samp = sum(sample_data$x5 * 1/sample_data$sampling_prob)
  T6_samp = sum(sample_data$x6 * 1/sample_data$sampling_prob)

  # Add missingness
  # Generate Item Nonresponse with Vectorized, Full Models 
  # Create a matrix of the predictor variables for efficient calculations

  X_matrix <- as.matrix(sample_data[, c("x1", "x2", "x3", "x4", "x5", "x6")])
  
  # Generate missingness indicators
  logit_p_miss_x1 <- params_miss$nr_int_1 + X_matrix %*% params_miss$nr_beta_1
  p_miss_x1 <- plogis(logit_p_miss_x1)
  is_missing_x1 <- rbinom(nrow(sample_data), 1, p_miss_x1)
  sample_data$x1[is_missing_x1 == 1] <- NA
  
  logit_p_miss_x2 <- params_miss$nr_int_2 + X_matrix %*% params_miss$nr_beta_2
  p_miss_x2 <- plogis(logit_p_miss_x2)
  is_missing_x2 <- rbinom(nrow(sample_data), 1, p_miss_x2)
  sample_data$x2[is_missing_x2 == 1] <- NA
  
  # MAR: x3, x4, x5, x6
  logit_p_miss_x3 <- params_miss$nr_int_3 + X_matrix[, -3] %*% params_miss$nr_beta_3
  p_miss_x3 <- plogis(logit_p_miss_x3)
  is_missing_x3 <- rbinom(nrow(sample_data), 1, p_miss_x3)
  sample_data$x3[is_missing_x3 == 1] <- NA
  
  logit_p_miss_x4 <- params_miss$nr_int_4 + X_matrix[, -4] %*% params_miss$nr_beta_4
  p_miss_x4 <- plogis(logit_p_miss_x4)
  is_missing_x4 <- rbinom(nrow(sample_data), 1, p_miss_x4)
  sample_data$x4[is_missing_x4 == 1] <- NA
  
  logit_p_miss_x5 <- params_miss$nr_int_5 + X_matrix[, -5] %*% params_miss$nr_beta_5
  p_miss_x5 <- plogis(logit_p_miss_x5)
  is_missing_x5 <- rbinom(nrow(sample_data), 1, p_miss_x5)
  sample_data$x5[is_missing_x5 == 1] <- NA
  
  logit_p_miss_x6 <- params_miss$nr_int_6 + X_matrix[, -6] %*% params_miss$nr_beta_6
  p_miss_x6 <- plogis(logit_p_miss_x6)
  is_missing_x6 <- rbinom(nrow(sample_data), 1, p_miss_x6)
  sample_data$x6[is_missing_x6 == 1] <- NA
  
  sample_result <- sample_data |> select(W, x1, x2, x3, x4, x5, x6) 
  
  ##############################################################
  ######################### MICE IMPS! #########################
  ##############################################################
  
  sample_data_for_mice <- sample_data |> select(W, x1, x2, x3, x4, x5, x6) 

  sample_data_for_mice$r1 = is_missing_x1
  sample_data_for_mice$r2 = is_missing_x2
  sample_data_for_mice$r3 = is_missing_x3
  sample_data_for_mice$r4 = is_missing_x4
  sample_data_for_mice$r5 = is_missing_x5
  sample_data_for_mice$r6 = is_missing_x6

  total_imputations_mice <- L
  
  sample_imp = mice(sample_data_for_mice, m = total_imputations_mice)
  mice_ALL_complete = complete(sample_imp, "long")
  
  start_index <- (total_imputations_mice - L) + 1
  end_index <- total_imputations_mice
  imputations_to_save <- start_index:end_index
  
  imputations_list <- lapply(imputations_to_save, function(j) complete(sample_imp, j))
  
  mice_result <- imputations_list

  ##############################################################
  ######################### JAGS IMPS! #########################
  ##############################################################
  
  # estimate variance
  tmp <- complete(mice(sample_data |> select(W, x1, x2, x3, x4, x5, x6), m = 1, seed = 1))
  V_x1_pop_HT <- sum((tmp$x1 / (1/tmp$W))^2 * (1 - (1/tmp$W)))
  V_x2_pop_HT <- sum((tmp$x2 / (1/tmp$W))^2 * (1 - (1/tmp$W)))
  
  X <- as.matrix(sample_data[, c("x1", "x2", "x3", "x4", "x5", "x6")])

  # prepare jags data
  
  jags_data <- list()
  jags_data$X <- X 
  jags_data$N <- nrow(jags_data$X)
  jags_data$W <- sample_data$W
  jags_data$R <- 1*!is.na(X)
  
  jags_data$T1_known <- sum(pop_data$x1)
  jags_data$T2_known <- sum(pop_data$x2)
  jags_data$V1_known <- V_x1_pop_HT
  jags_data$V2_known <- V_x2_pop_HT
  
  # Initialize function for JAGS
  inits_function <- function() {
    mice_imp <- mice(jags_data$X, m = 1, maxit = 5, printFlag = FALSE)
    X_mice_complete <- complete(mice_imp)
    X_mice_numeric <- do.call(cbind, lapply(X_mice_complete, as.numeric))
    
    X_init <- matrix(NA, nrow = nrow(jags_data$X), ncol = ncol(jags_data$X))
    
    is_missing <- is.na(jags_data$X)
    
    X_init[is_missing] <- X_mice_numeric[is_missing]
    
    return(list(X = X_init))
  }
  
  # jags model setup
  params_to_save <- c(
  paste0("alpha", 1:6), "beta1w",
  paste0("beta2", 1),
  paste0("beta3", 1:2),
  paste0("beta4", 1:3),
  paste0("beta5", 1:4),
  paste0("beta6", 1:5),
  paste0("sigma", 5:6),
  paste0("tau",   5:6),
  
  # Non-response model priors
  paste0("nr_int_", 1:6),
  paste0("nr_beta_1_", 1:6),
  paste0("nr_beta_2_", 1:6),
  paste0("nr_beta_3_", 1:5),
  paste0("nr_beta_4_", 1:5),
  paste0("nr_beta_5_", 1:5),
  paste0("nr_beta_6_", 1:5),
  
  # totals
  "T1_imputed_total", "T2_imputed_total", "T3_imputed_total", 
  "T4_imputed_total", "T5_imputed_total", "T6_imputed_total")

  #### run sampler - fit params

  # The jags() call with command line arguments
  jags_fit_item <- jagsUI::jags(
    data = jags_data,
    inits = inits_function,
    parameters.to.save = params_to_save,
    model.file = "item_nonresponse_model.bug",
    n.chains = n_chains,
    n.adapt = 1000,
    n.iter = n_iter,
    n.burnin = n_burnin,
    parallel = TRUE
  )
  
  summary_matrix <- jags_fit_item$summary
  summary_df <- as.data.frame(summary_matrix)
  summary_df$parameter <- rownames(summary_matrix)
  summary_df <- summary_df[, c(ncol(summary_df), 1:(ncol(summary_df)-1))]
  
  jags_filename <- sprintf("jags_summary_samp_%d.csv", i)
  jags_save_path <- file.path(folder_name, jags_filename)
  
  print(paste("Saving JAGS summary for sample", i, "to:", jags_save_path))
  write.csv(summary_df, file = jags_save_path, row.names = FALSE)
  
  # --- Generate L Imputations ---
  # Now that the model is converged, run L more iterations 
  # from where it left off and ONLY save the "X" matrix.
  jags_imputations <- update(
    jags_fit_item, 
    parameters.to.save = "X",
    n.iter = L,
    n.thin = 1
  )
  
  imputed_X_array <- jags_imputations$sims.list$X
  
  # Total imputations generated is (L * n_chains)
  n_total_generated <- dim(imputed_X_array)[1] 
  
  W_col <- sample_data$W
  imputations_list <- list()
  
  # Get L evenly spaced indices from the total pool
  # This "thins" the samples and gives you exactly L
  indices_to_save <- round(seq(1, n_total_generated, length.out = L))

  for (j in 1:L) { # Loop exactly L times
    
    # Get the j-th index from our spaced-out list
    imp_index <- indices_to_save[j] 
    
    imp_matrix <- imputed_X_array[imp_index, , ] 
    imp_df <- as.data.frame(imp_matrix)
    colnames(imp_df) <- c("x1", "x2", "x3", "x4", "x5", "x6")
    imp_df$W <- W_col
    imp_df <- imp_df[, c("W", "x1","x2", "x3", "x4", "x5", "x6")]
    imputations_list[[j]] <- imp_df
  }
  
  jags_result <- imputations_list
  print(paste("Stored", L, "JAGS imputations for sample", i)) # Prints L
  
  #############################################################################
  ##################### PLOT GENERATION #######################################
  #############################################################################
  
  print(paste("Generating plots for sample", i))

  pdf_save_path <- file.path(folder_name, sprintf("imputation_plots_samp_%d.pdf", i))

  pdf(pdf_save_path, width = 12, height = 16)
  par(mfrow = c(6, 2), mar = c(4, 4, 3, 1)) 
  
  current_mice_imputations <- mice_result
  
  sum_stats_list <- lapply(current_mice_imputations, function(imp_df) {
    data.frame(
      sx1 = sum(imp_df$W * imp_df$x1, na.rm = TRUE),
      sx2 = sum(imp_df$W * imp_df$x2, na.rm = TRUE),
      sx3 = sum(imp_df$W * imp_df$x3, na.rm = TRUE),
      sx4 = sum(imp_df$W * imp_df$x4, na.rm = TRUE),
      sx5 = sum(imp_df$W * imp_df$x5, na.rm = TRUE),
      sx6 = sum(imp_df$W * imp_df$x6, na.rm = TRUE)
    )
  })
  
  sum_stats <- dplyr::bind_rows(sum_stats_list)
  
  data1_mice <- sum_stats$sx1
  data2_mice <- sum_stats$sx2
  data3_mice <- sum_stats$sx3
  data4_mice <- sum_stats$sx4
  data5_mice <- sum_stats$sx5
  data6_mice <- sum_stats$sx6
  
  data1_jags <- jags_fit_item$sims.list$T1_imputed_total
  data2_jags <- jags_fit_item$sims.list$T2_imputed_total
  data3_jags <- jags_fit_item$sims.list$T3_imputed_total
  data4_jags <- jags_fit_item$sims.list$T4_imputed_total
  data5_jags <- jags_fit_item$sims.list$T5_imputed_total
  data6_jags <- jags_fit_item$sims.list$T6_imputed_total
  
  # --- 5. Generate Plots in Pairs ---
  
  # --- Pair 1: x1 ---
  plot_xlim1 <- range(c(data1_mice, data1_jags, T1_known, T1_samp), na.rm = TRUE)
  hist(data1_mice, xlim = plot_xlim1, xlab = "Sum of x1", main = "MICE: sx1")
  abline(v = T1_known, col = "red", lwd = 2)
  abline(v = T1_samp, col = "blue", lwd = 2)
  
  hist(data1_jags, xlim = plot_xlim1, col = "lightblue", xlab = "Imputed Total T1", main = "JAGS: T1")
  abline(v = T1_known, col = "red", lwd = 2)
  abline(v = T1_samp, col = "blue", lwd = 2)
  
  # --- Pair 2: x2 ---
  plot_xlim2 <- range(c(data2_mice, data2_jags, T2_known, T2_samp), na.rm = TRUE)
  hist(data2_mice, xlim = plot_xlim2, xlab = "Sum of x2", main = "MICE: sx2")
  abline(v = T2_known, col = "red", lwd = 2)
  abline(v = T2_samp, col = "blue", lwd = 2)
  
  hist(data2_jags, xlim = plot_xlim2, col = "lightblue", xlab = "Imputed Total T2", main = "JAGS: T2")
  abline(v = T2_known, col = "red", lwd = 2)
  abline(v = T2_samp, col = "blue", lwd = 2)
  
  # --- Pair 3: x3 ---
  plot_xlim3 <- range(c(data3_mice, data3_jags, T3_known, T3_samp), na.rm = TRUE)
  hist(data3_mice, xlim = plot_xlim3, xlab = "Sum of x3", main = "MICE: sx3")
  abline(v = T3_known, col = "red", lwd = 2)
  abline(v = T3_samp, col = "blue", lwd = 2)
  
  hist(data3_jags, xlim = plot_xlim3, col = "lightblue", xlab = "Imputed Total T3", main = "JAGS: T3")
  abline(v = T3_known, col = "red", lwd = 2)
  abline(v = T3_samp, col = "blue", lwd = 2)
  
  # --- Pair 4: x4 ---
  plot_xlim4 <- range(c(data4_mice, data4_jags, T4_known, T4_samp), na.rm = TRUE)
  hist(data4_mice, xlim = plot_xlim4, xlab = "Sum of x4", main = "MICE: sx4")
  abline(v = T4_known, col = "red", lwd = 2)
  abline(v = T4_samp, col = "blue", lwd = 2)
  
  hist(data4_jags, xlim = plot_xlim4, col = "lightblue", xlab = "Imputed Total T4", main = "JAGS: T4")
  abline(v = T4_known, col = "red", lwd = 2)
  abline(v = T4_samp, col = "blue", lwd = 2)
  
  # --- Pair 5: x5 ---
  plot_xlim5 <- range(c(data5_mice, data5_jags, T5_known, T5_samp), na.rm = TRUE)
  hist(data5_mice, xlim = plot_xlim5, xlab = "Sum of x5", main = "MICE: sx5")
  abline(v = T5_known, col = "red", lwd = 2)
  abline(v = T5_samp, col = "blue", lwd = 2)
  
  hist(data5_jags, xlim = plot_xlim5, col = "lightblue", xlab = "Imputed Total T5", main = "JAGS: T5")
  abline(v = T5_known, col = "red", lwd = 2)
  abline(v = T5_samp, col = "blue", lwd = 2)
  
  # --- Pair 6: x6 ---
  plot_xlim6 <- range(c(data6_mice, data6_jags, T6_known, T6_samp), na.rm = TRUE)
  hist(data6_mice, xlim = plot_xlim6, xlab = "Sum of x6", main = "MICE: sx6")
  abline(v = T6_known, col = "red", lwd = 2)
  abline(v = T6_samp, col = "blue", lwd = 2)
  
  hist(data6_jags, xlim = plot_xlim6, col = "lightblue", xlab = "Imputed Total T6", main = "JAGS: T6")
  abline(v = T6_known, col = "red", lwd = 2)
  abline(v = T6_samp, col = "blue", lwd = 2)
  
  # --- 6. Close PDF Device ---
  dev.off()
  
  print(paste("Plots saved to:", pdf_save_path))
  
  # Return results for this sample
  return(list(
    sample = sample_result,
    mice_result = mice_result,
    jags_result = jags_result
  ))
}

# Determine number of cores to use
# Be conservative: if JAGS uses n_chains cores internally, 
# use fewer cores for the outer loop to avoid oversubscription
max_cores <- detectCores()
if (n_chains > 1) {
  # Use fewer cores for outer loop when JAGS is using multiple chains
  outer_cores <- min(max(1, max_cores %/% n_chains), n_samp)
} else {
  # Use more cores when JAGS is single-threaded
  outer_cores <- min(max_cores - 1, n_samp)
}

print(paste("Using", outer_cores, "cores for parallel processing"))
print(paste("JAGS will use", n_chains, "chains internally"))

# Run samples sequentially for debugging
results_list <- lapply(1:n_samp, process_sample)

# Extract results from parallel processing
SAMPLES <- list()
MICE_RESULTS <- list()
JAGS_RESULTS <- list()

for (i in 1:n_samp) {
  SAMPLES[[i]] <- results_list[[i]]$sample
  MICE_RESULTS[[i]] <- results_list[[i]]$mice_result
  JAGS_RESULTS[[i]] <- results_list[[i]]$jags_result
}

# Save the results lists to the folder
print("Saving SAMPLES list...")
samples_save_path <- file.path(folder_name, "samples_list.rds")
saveRDS(SAMPLES, file = samples_save_path)
print(paste("SAMPLES saved to:", samples_save_path))

print("Saving MICE_RESULTS list...")
mice_results_save_path <- file.path(folder_name, "mice_results_list.rds")
saveRDS(MICE_RESULTS, file = mice_results_save_path)
print(paste("MICE_RESULTS saved to:", mice_results_save_path))

print("Saving JAGS_RESULTS list...")
jags_results_save_path <- file.path(folder_name, "jags_results_list.rds")
saveRDS(JAGS_RESULTS, file = jags_results_save_path)
print(paste("JAGS_RESULTS saved to:", jags_results_save_path))

print("All results saved successfully!")
print(paste("Script completed at:", Sys.time()))
