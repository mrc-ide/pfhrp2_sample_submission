# 04.sample_size.R
#
# Author: Bob Verity
# Date: 2023-06-02
#
# Inputs:
# df_hill_fits.rds, produced by 02.fit_hill_coefficients.R
#
# Outputs:
# df_power.rds: power curves for all parameter combinations
# 
# Purpose:
# Produces power curves from fitted hill coefficients and saves results and
# plots to file.
#
# ------------------------------------------------------------------

library(tidyverse)

# hill function
get_model_power <- function(A, B, K, w, N) {
  #A + (B - A) / (1 + K*N)^(w/K)
  A + (B - A)*(1 - exp(-K*N))^w
}

# load parameter fits from file
df_hill_fits <- readRDS("threshold_analysis/outputs/df_hill_fits.rds")

# make power curves for all parameter combinations
df_power <- df_hill_fits %>%
  expand_grid(N = 5:2e3) %>%
  mutate(power = get_model_power(A, B, K, w, N))

# define target power
target_power <- 0.8

# find minimum sample size to achieve target power
df_ss <- df_power %>%
  group_by(n_clust, prevalence, ICC, alpha, prior_prev_shape1, prior_prev_shape2,
           prior_ICC_shape1, prior_ICC_shape2, CrI_type, n_intervals, round_digits,
           rejection_threshold, seed, prev_thresh) %>%
  summarise(N = ifelse(any(power >= target_power), N[which(power >= target_power)[1]], NA)) %>%
  ungroup() %>%
  mutate(N = ifelse(N < 2e3, N, NA),
         N = ifelse(N > 5, N, 5))

# checks on sample sizes:
# ideally should be decreasing with n_clust
df_ss %>%
  group_by_at(vars(-c(n_clust, N_opt))) %>%
  summarise(pass = all(diff(N_opt) <= 0)) %>%
  filter(!pass) %>%
  nrow()
  #left_join(df_ss) %>%
  #pivot_wider(names_from = n_clust, values_from = N) %>%
  #View()

# ideally should be decreasing with prevalence (above threshold)
df_ss %>%
  mutate(N = ifelse(is.na(N), 2e3, N)) %>%
  group_by_at(vars(-c(prevalence, N))) %>%
  summarise(pass = all(diff(N) <= 0)) %>%
  filter(!pass) %>%
  nrow()
  #left_join(df_ss) %>%
  #pivot_wider(names_from = prevalence, values_from = N) %>%
  #View()

# ideally should be increasing with prev_thresh
df_ss %>%
  mutate(N = ifelse(is.na(N), 2e3, N)) %>%
  group_by_at(vars(-c(prev_thresh, N))) %>%
  summarise(pass = all(diff(N) >= 0)) %>%
  filter(!pass) %>%
  nrow()
  #left_join(df_ss) %>%
  #pivot_wider(names_from = ICC, values_from = N) %>%
  #View()

# write to file
if (FALSE) {
  saveRDS(df_ss, file = "threshold_analysis/outputs/df_ss.rds")
}

df_ss %>%
  filter(ICC == 0.05) %>%
  filter(prior_ICC_shape2 == 1) %>%
  filter(prev_thresh == 0.05) %>%
  select(n_clust, prevalence, N) %>%
  pivot_wider(names_from = prevalence, values_from = N) %>%
  View()

# filter and pivot wide
tab1 <- df_ss %>%
  filter(ICC == 0.05) %>%
  filter(prior_ICC_shape2 == 1) %>%
  filter(prev_thresh == 0.05) %>%
  select(n_clust, prevalence, N) %>%
  pivot_wider(names_from = prevalence, values_from = N)

View(tab1)

tab2 <- df_ss %>%
  filter(ICC == 0.04) %>%
  filter(prior_ICC_shape2 == 1) %>%
  filter(prev_thresh == 0.1) %>%
  select(n_clust, prevalence, N) %>%
  pivot_wider(names_from = prevalence, values_from = N)

if (FALSE) {
  write.csv(rbind(tab1, tab2), file = "threshold_analysis/outputs/df_tables.csv",
            quote = FALSE, row.names = FALSE)
}
