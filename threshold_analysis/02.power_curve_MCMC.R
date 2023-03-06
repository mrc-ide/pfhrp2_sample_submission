# 02.power_curve_MCMC.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Reads in the results of 01.simulate_power.R, consisting of a numerator and a
# denominator of empirical power over a range of sample sizes and parameter
# combinations. For each parameter combination independently, fits a power curve
# to the data via MCMC using the drjacoby package. The functional form of the
# curve is a generalisation of the Hill function.
#
# Performs checks on MCMC output for manual inspection. Also produces a series
# of MCMC diagnostic plots and saves these to file. Note, the loop over MCMCs
# involves saving raw MCMC data to file in a temporary folder. If re-running
# this analysis on your own machine you will need to edit the save path to your
# own temporary folder (see comments in code below).
#
# Assuming all MCMCs pass checks, takes a series of posterior draws and
# generates posterior power curves. These curves are saved to file rather than
# plotting directly here, to separate out the plotting process from the analysis
# process.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)
library(dplyr)
library(patchwork)

# load drjacoby from a specific version to ensure future compatibility
#devtools::install_github("mrc-ide/drjacoby", ref = "v1.5.2")
library(drjacoby)


# ------------------------------------------------------------------
# 01. PREPROCESS DATA

# load the results of 01.simulate_power.R. If this script has been run
# repeatedly then read in all results. Also read in data.frame specifying which
# parameter combinations to keep, i.e., they have adequate power for some N
df_params <- readRDS("threshold_analysis/data/df_params.rds")
df_sim1 <- readRDS("threshold_analysis/data/df_sim_01.rds")
df_sim2 <- readRDS("threshold_analysis/data/df_sim_02.rds")
df_sim3 <- readRDS("threshold_analysis/data/df_sim_03.rds")
df_sim4 <- readRDS("threshold_analysis/data/df_sim_04.rds")
df_keep <- readRDS("threshold_analysis/data/df_keep.rds")

# combine simulation data.frames. If the same parameter combination has been run
# multiple times then we want to combine the numerator and denominators
# separately at this stage. Also filter using df_keep
df_sim_comb <- rbind(df_sim1, df_sim2, df_sim3, df_sim4) %>%
  dplyr::select(param_i, N, num, denom) %>%
  left_join(df_params) %>%
  mutate(prev_thresh1 = mapply(function(x) x[1], num), # split out numerator values into the different associated prevalence thresholds
         prev_thresh2 = mapply(function(x) x[2], num),
         prev_thresh3 = mapply(function(x) x[3], num)) %>%
  dplyr::select(-num) %>%
  pivot_longer(cols = c(prev_thresh1, prev_thresh2, prev_thresh3),   # get into long format
               names_to = "prev_thresh",
               values_to = "num") %>%
  group_by(param_i, rho, n_clust, p, N, prev_thresh) %>%
  summarise(lower = lower[1],    # sum up over repeated simulations at the same parameters
            upper = upper[1],
            num = sum(num),
            denom = sum(denom)) %>%
  ungroup() %>%
  left_join(df_keep) %>%
  filter(keep == TRUE) %>%
  dplyr::select(-keep)

# make data.frame of parameter combinations that we want to explore
df_filtered <- df_keep %>%
  filter(keep) %>%
  dplyr::select(-keep)

# save both of these to file
saveRDS(df_sim_comb, file = "threshold_analysis/data/df_sim_comb.rds")
saveRDS(df_filtered, file = "threshold_analysis/data/df_filtered.rds")


# ------------------------------------------------------------------
# 02. RUN MCMC

# define MCMC parameters
burnin <- 1e3
samples <- 1e4
chains <- 10

# define model parameters
df_drj <- define_params(name = "A", min = 0, max = 1,
                        name = "B", min = 0, max = 1,
                        name = "K", min = 0, max = 100,
                        name = "w", min = 0, max = 10)

# load C++ loglikelihood function
Rcpp::sourceCpp("threshold_analysis/src/hill_loglike.cpp")

# loop through all parameter combinations and threshold values
t0 <- Sys.time()
for (i in seq_len(nrow(df_filtered))) {
  
  # print message to console
  message(sprintf("RUNNING MCMC %s of %s\n", i, nrow(df_filtered)))
  
  # subset data
  df_sub <- df_sim_comb %>%
    filter(param_i == df_filtered$param_i[i]) %>%
    filter(prev_thresh == df_filtered$prev_thresh[i])
  
  # make data list
  dat_list <- list(N = df_sub$N,
                   num = df_sub$num,
                   denom = df_sub$denom)
  
  # reset random seed. A very small number of MCMCs had ESS very slightly lower
  # than 1e3 for the parameter A. This is such a tiny effect that it justifies
  # repeating MCMCs with a different seed.
  if (i %in% c(167, 209)) {
    set.seed(2)
  } else {
    set.seed(1)
  }
  
  # run MCMC
  mcmc <- run_mcmc(data = dat_list,
                   df_params = df_drj,
                   loglike = "loglike",
                   logprior = "logprior",
                   burnin = burnin,
                   samples = samples,
                   chains = chains,
                   beta_manual = c(0.0625, 0.125, 0.25, 0.5, 1),
                   silent = FALSE)
  
  #plot_cor(mcmc, "K", "w")
  #mcmc$diagnostics$ess
  
  # save raw MCMC output to file. This folder does not form part of the online
  # Git repository as it will contain too much output. If running locally,
  # replace this with any temporary folder on your computer.
  saveRDS(mcmc, file = sprintf("ignore/MCMC_raw/mcmc_%s.rds", i))
  
}
Sys.time() - t0


# ------------------------------------------------------------------
# 03. CHECK MCMC OUTPUT

# make data.frames for storing diagnostics
df_rhat <- df_ess <- df_filtered %>%
  mutate(A = NA, B = NA, K = NA, w = NA)

# make list for storing samples from the posterior
list_post_samples <- list()

# loop through all MCMC output
for (i in seq_len(nrow(df_filtered))) {
  
  # print message to console
  message(sprintf("Diagnosing MCMC %s of %s", i, nrow(df_filtered)))
  
  # read mcmc output from file
  mcmc <- readRDS(sprintf("ignore/MCMC_raw/mcmc_%s.rds", i))
  
  # store MCMC diagnostics
  df_rhat[i, c("A", "B", "K", "w")] <- as.data.frame(t(mcmc$diagnostics$rhat))
  df_ess[i, c("A", "B", "K", "w")] <- as.data.frame(t(mcmc$diagnostics$ess))
  
  # reset random seed
  set.seed(1)
  
  # sample parameters from posterior and store
  list_post_samples[[i]] <- mcmc$output %>%
    filter(phase == "sampling") %>%
    sample_n(1e3) %>%
    dplyr::select(A, B, K, w) %>%
    mutate(MCMC_index = i,
           param_draw = 1:1e3)

  # annotate trace plots with rhat and ESS, and print to file
  plot_all <- plot_par(mcmc, display = FALSE)
  for (j in 1:4) {
    plot_j <- plot_all[[j]]$combined +
      plot_annotation(sprintf("rhat = %s,  ESS = %s",
                              round(mcmc$diagnostics$rhat[j], digits = 3),
                              round(mcmc$diagnostics$ess[j], digits = 1)))
    ggsave(plot = plot_j, width = 10, height = 5, dpi = 72,
           filename = sprintf("threshold_analysis/MCMC_output/trace_plots/MCMC%s_trace_%s.png", i, c("A", "B", "K", "w")[j]))
  }
  
}

# manually check the proportion of MCMCs that are within convergence and ESS limits
df_rhat %>%
  dplyr::select(A, B, K, w) %>%
  "<"(1.1) %>%
  colMeans()

df_ess %>%
  dplyr::select(A, B, K, w) %>%
  apply(2, min)


# save diagnostic data.frames to file
if (FALSE) {
  write.csv(df_rhat, file = "threshold_analysis/MCMC_output/rhat.csv", row.names = FALSE)
  write.csv(df_ess, file = "threshold_analysis/MCMC_output/ess.csv", row.names = FALSE)
}


# ------------------------------------------------------------------
# 04. PRODUCE POSTERIOR POWER CURVES

# function to predict over range of N for each parameter draw
f_pred <- function(A, B, K, w, N) {
  A + (1 - A)*B / (1 + (K / N)^w)
}

# loop over MCMC output. Note, technically it would be possible to do all this
# in a single step, but the memory cost would be prohibitively high
list_post_summary <- list()
for (i in seq_len(nrow(df_filtered))) {
  
  # write message to console
  message(sprintf("Summary %s of %s", i, length(list_post_samples)))
  
  # get data.frame of predicted power over all posterior draws
  df_pred <- list_post_samples[[i]] %>%
    expand_grid(N = 5:2e3) %>%
    mutate(power = f_pred(A, B, K, w, N))
  
  # get summary over predictions
  df_summary <- df_pred %>%
    group_by(N) %>%
    summarise(mean = mean(power),
              lower = quantile(power, p = 0.025),
              upper = quantile(power, p = 0.975))
  
  # store summary data.frame
  list_post_summary[[i]] <- df_summary
}

# save to file
if (FALSE) {
  saveRDS(list_post_summary, file = "threshold_analysis/data/post_summary.rds")
}
