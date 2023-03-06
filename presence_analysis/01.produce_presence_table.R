# 01.produce_presence_table.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Produces sample size table for the "presence" analysis, in which our null
# hypothesis is that there are zero deletions in our region.
#
# This is a relatively simple analysis that requires no simulation. Power can be
# calculated exactly, and the optimal sample size can be determined by brute
# force.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# make data.frame of all parameter combinations that we want to explore
df_params <- expand_grid(rho = c(0.05, 0.2),
                         p = seq(0.01, 0.05, 0.01),
                         n_clust = 2:10,
                         N = NA)

# defined fixed terms
target_power <- 0.8
N_vec <- 5:2e3

# loop through all combinations
for (i in seq_len(nrow(df_params))) {
  
  # calculate power over range of N
  alpha <- df_params$p[i] * (1 / df_params$rho[i] - 1)
  beta <- (1 - df_params$p[i]) * (1 / df_params$rho[i] - 1)
  log_prob_zero <- df_params$n_clust[i] * (lgamma(alpha + beta) - lgamma(beta) + lgamma(N_vec + beta) - lgamma(N_vec + alpha + beta))
  power <- 1 - exp(log_prob_zero)
  
  # find first N that exceeds target power
  if (any(power > target_power)) {
    df_params$N[i] <- N_vec[which(power > target_power)[1]]
  }
}

# get into wide format for each value of rho
df_rho_0.05 <- df_params %>%
  filter(rho == 0.05) %>%
  select(p, n_clust, N) %>%
  pivot_wider(names_from = p, values_from = N)

df_rho_0.20 <- df_params %>%
  filter(rho == 0.20) %>%
  select(p, n_clust, N) %>%
  pivot_wider(names_from = p, values_from = N)

df_rho_0.05
df_rho_0.20

# save to file
write.csv(df_rho_0.05, file = "presence_analysis/tables/sample_size_ICC_05.csv", row.names = FALSE)
write.csv(df_rho_0.20, file = "presence_analysis/tables/sample_size_ICC_20.csv", row.names = FALSE)
