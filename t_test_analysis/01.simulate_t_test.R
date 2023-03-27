# 01.simulate_t_test.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Compares power for the t-test vs. the DRpower method. For the t-test, power is
# calculated in two ways: 1) by simulation (empirical power), 2) using the
# analytical formula. This is because the analytical formula cannot fully
# propagate uncertainty (i.e. error in the estimate of prevalence is present in
# the estimate of the variance, etc.), meaning the empirical power only
# converges on this result for large cluster numbers. The power by the DRpower
# method is estimated in the usual way, by simulation.
#
# Results are saved to file to be visualised elsewhere.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# load DRpower package from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# define parameters
N <- 50
p_thresh <- 0.05
rho <- 0.05
reps <- 1e4
t_conf <- 0.25

# create data.frame with parameter combinations to explore, and with empty slots
# for all future outputs
df_pow <- expand_grid(n_clust = c(2, 10),
                      p = seq(0.05, 0.2, 0.005)) %>%
  mutate(t_power = NA,
         t_num = NA,
         t_denom = reps,
         t_analytical = NA,
         DR_power = NA,
         DR_lower = NA,
         DR_upper = NA)

# get empirical power of t-test for all combinations
for (i in 1:nrow(df_pow)) {
  message(sprintf("%s of %s", i, nrow(df_pow)))
  
  # simulate from beta-binomial
  sim_mat <- rbbinom_reparam(n_clust = df_pow$n_clust[i] * reps,
                             N = N, p = df_pow$p[i], rho = rho) %>%
    "/"(N) %>%
    matrix(nrow = reps)
  
  # identify simulations with zero variance, in which case we cannot perform
  # t-test
  zero_var <- (apply(sim_mat, 1, var) == 0)
  
  # otherwise perform one-sample t-test
  p_vec <- sim_mat[!zero_var,] %>%
    apply(1, function(x) {
      t.test(x, mu = p_thresh, alternative = "greater")$p.value
    })
  
  # calculate empirical power
  df_pow$t_num[i] <- sum(p_vec < t_conf)
  df_pow$t_power[i] <- mean(p_vec < t_conf)
}

# fill in column for analytical solution to t-test power
df_pow <- df_pow %>%
  mutate(t_analytical = power.t.test(n = n_clust, delta = p - p_thresh,
                                     sd = sqrt(p*(1 - p) * (1 / N + rho)),
                                     type = "one.sample", sig.level = 0.05,
                                     alternative = "one.sided")$power)

# get empirical power of DRpower method for all combinations
t0 <- Sys.time()
for (i in 1:nrow(df_pow)) {
  message(sprintf("%s of %s", i, nrow(df_pow)))
  
  DR_est <- get_power(rep(N, df_pow$n_clust[i]), prevalence = df_pow$p[i],
                      ICC = rho, prior_ICC_shape2 = 19, reps = reps)
  df_pow$DR_power[i] <- DR_est$power
  df_pow$DR_lower[i] <- DR_est$lower
  df_pow$DR_upper[i] <- DR_est$upper
}
Sys.time() - t0

# save results to file
if (FALSE) {
  saveRDS(df_pow, file = "t_test_analysis/data/df_pow.rds")
}

