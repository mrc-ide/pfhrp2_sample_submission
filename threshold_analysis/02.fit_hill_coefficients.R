# 02.fit_hill_coefficients.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Inputs:
# df_sim.rds, produced by 01.simulate_power.R
#
# Outputs:
# df_hill_fits.rds: fitted coefficients for each parameter combination
#
# Purpose:
# Fits a hill function to the curve of power vs. N for each parameter
# combination. The hill function is a 4 parameter model and is fitted via
# maximum likelihood.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# load DRpower package from local version
# TODO - replace with tagged online version once ready
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# hill function
get_model_power <- function(A, B, K, w, N) {
  #A + (B - A) / (1 + K*N)^(w/K)
  #A + (B - A) / (1 + K*N)^w
  #A + (B - A) / (1 + (K / N))^w
  A + (B - A)*(1 - exp(-K*N))^w
}

# run optim a single time to get maximum likelihood
get_param_optim0 <- function(N, n_pos, reps, theta0 = c(0, 0, 1, 0)) {
  
  ret <- optim(par = theta0,
        fn = function(theta) {
          A <- 1 / (1 + exp(-theta[1]))  # some values are transformed to allow fitting to occur over the real line
          B <- 1 / (1 + exp(-theta[2]))
          K <- exp(theta[3])
          w <- 3 / (1 + exp(-theta[4]))
          
          power_mod <- get_model_power(A, B, K, w, N)
          ll <- sum(dbinom(n_pos, size = reps, prob = power_mod, log = TRUE))
          
          -ll
        })
  
  return(ret)
}

# run optimisation multiple times from random initial conditions
get_param_optim <- function(N, n_pos, reps, random_restart = 10) {
  
  optim_list <- list()
  for (i in 1:random_restart) {
    optim_list[[i]] <- get_param_optim0(N, n_pos, reps, theta0 = rnorm(4, sd = 1))
  }
  w <- which.min(mapply(function(x) x$value, optim_list))
  
  return(optim_list[[w]])
}

# ------------------------------------------------------------------

# read in simulation results
df_sim <- readRDS("threshold_analysis/outputs/df_sim.rds") %>%
  filter(prior_ICC_shape2 == 1) %>%
  filter(n_clust > 1)

df_sim <- df_sim %>%
  filter(ICC == 0.05) %>%
  filter(prev_thresh == 0.05)

# get unique parameter groups and add placeholder for hill function parameters
df_hill_fits <- df_sim %>%
  group_by(n_clust, prevalence, ICC, alpha, prior_prev_shape1, prior_prev_shape2,
           prior_ICC_shape1, prior_ICC_shape2, CrI_type, n_intervals, round_digits,
           rejection_threshold, seed, prev_thresh) %>%
  summarise(group = cur_group_id()) %>%
  ungroup() %>%
  mutate(A = NA,
         B = NA,
         K = NA,
         w = NA)

# split data into list of data.frames for each group
tmp <- df_sim %>%
  left_join(df_hill_fits)
df_list <- split(x = tmp, f = tmp$group)  

# manually set i to look at specific plots
if (FALSE) {
  i <- which(df_hill_fits$n_clust == 11 &
               abs(df_hill_fits$prevalence - 0.06) < 1e-6 &
               df_hill_fits$ICC == 0.0 &
               df_hill_fits$prev_thresh == 0.05)
}

# fit hill function parameters for each group by maximum likelihood
for (i in seq_along(df_list)) {
  message(sprintf("%s of %s", i, length(df_list)))
  
  N <- df_list[[i]]$N
  reps <- df_list[[i]]$reps
  n_pos <- df_list[[i]]$power_mean * reps
  
  # find maximum likelihood vector of parameters
  set.seed(1)
  par_optim <- get_param_optim(N, n_pos, reps, random_restart = 100)
  
  # extract parameter values on final scale
  A <- 1 / (1 + exp(-par_optim$par[1]))
  B <- 1 / (1 + exp(-par_optim$par[2]))
  K <- exp(par_optim$par[3])
  w <- 3 / (1 + exp(-par_optim$par[4]))
  
  # optional plot
  if (F) {
    df_list[[i]] %>%
      ggplot() + theme_bw() +
      geom_errorbar(aes(x = N, y = power_mean, ymin = power_lower, ymax = power_upper)) +
      #geom_point(aes(x = N, y = n_pos / reps)) +
      geom_line(aes(x = N, y = power), col = "red",
                data = data.frame(N = 5:2e3,
                                  power = get_model_power(A, B, K, w, 5:2e3))) +
      ylim(c(0, 1)) + xlim(c(0, 2000)) +
      ylab("Power")
  }
  
  # store values
  df_hill_fits$A[i] <- A
  df_hill_fits$B[i] <- B
  df_hill_fits$K[i] <- K
  df_hill_fits$w[i] <- w
}

# drop groups column
df_hill_fits <- df_hill_fits %>%
  select(-group)

# save fitted hill function parameters
if (FALSE) {
  saveRDS(df_hill_fits, file = "threshold_analysis/outputs/df_hill_fits.rds")
}
