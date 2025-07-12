# 01.calculate_power.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs:
# - 03.ICC_posteriors/dat_filter1.rds
#
# Outputs:
# - 04.retrospective_power/outputs/retro_power.rds
#
# Purpose:
# Calculates power on historical studies. Uses two methods, one based on z-test and another using DRpower Bayesian approach. Saves calculated power to file for plotting later.
#
# ------------------------------------------------------------------
library(tidyverse)
#devtools::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)

# read in pre-filtered data from previous script output
dat <- readRDS("03.ICC_posteriors/outputs/dat_filter1.rds")

# split into distinct analysis groups
dat_list <- dat |>
  group_split(study_ID, COUNTRY_NAME, ADMIN1_NAME, YEAR_START)
n_groups <- length(dat_list)

# calculate power in 3 ways
pow_z_Dplugin <- rep(NA, n_groups)
pow_z_Dest <- rep(NA, n_groups)
pow_Bayes <- rep(NA, n_groups)

for (i in 1:n_groups) {
  message(sprintf("%s of %s", i, n_groups))
  
  p <- 0.08
  p_thresh <- 0.05
  N <- dat_list[[i]]$HRP2_TESTED
  N_mean <- sum(N^2) / sum(N)
  Deff <- 1 + (N_mean - 1)*0.05
  
  pow_z_Dplugin[i] <- 1e2 * pnorm(abs(p - p_thresh) / sqrt(1.5 * p*(1 - p) / sum(N)) - qnorm(1 - 0.05/2))
  pow_z_Dest[i] <- 1e2 * pnorm(abs(p - p_thresh) / sqrt(Deff * p*(1 - p) / sum(N)) - qnorm(1 - 0.05/2))
  
  pow_Bayes[i] <- get_power_threshold(N = N, reps = 1e3, silent = TRUE)$power
}

# save results to file
ret <- data.frame(pow_z_Dplugin, pow_z_Dest, pow_Bayes)
saveRDS(ret, file = "04.retrospective_power/outputs/retro_power.rds")


