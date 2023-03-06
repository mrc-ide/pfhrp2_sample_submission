# 03.plot_power_curve.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Reads in the results of 02.power_curve_MCMC and produces plots of fitted powe
# curves. Also determines the appropriate sample size to hit a target power and
# writes table to file.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)
library(dplyr)
library(cowplot)

# ------------------------------------------------------------------
# 01. GET FITTED SAMPLE SIZE

# fixed parameters
target_power <- 0.8

# fixed plotting parameters
ribbon_col <- "royalblue3"
ribbon_alpha <- 0.5
error_width <- 2
point_size <- 0.5

# load fitted curves and associated data
list_post_summary <- readRDS("threshold_analysis/data/post_summary.rds")
df_sim_comb <- readRDS("threshold_analysis/data/df_sim_comb.rds")
df_keep <- readRDS("threshold_analysis/data/df_keep.rds")
df_filtered <- readRDS("threshold_analysis/data/df_filtered.rds")

# append df_filtered with final fitted sample size
df_filtered <- df_filtered %>%
  mutate(N_fit = NA,
         N_max = NA)
for (i in seq_along(list_post_summary)) {
  post_mean <- list_post_summary[[i]]$mean
  N_fit <- list_post_summary[[i]]$N[which(post_mean > target_power)[1]]
  df_filtered$N_fit[i] <- N_fit
  df_filtered$N_max[i] <- max(min(2*N_fit, 2e3), 100)
}

# get into wide format and save each table separately
table_final <- list()
for (j in 1:3) {
  table_final[[j]] <- df_filtered %>%
    arrange(p) %>%
    filter(prev_thresh == sprintf("prev_thresh%s", j)) %>%
    dplyr::select(n_clust, p, N_fit) %>%
    pivot_wider(names_from = p, values_from = N_fit) %>%
    arrange(n_clust)
}

# save tables to file
for (j in 1:3) {
  write.csv(table_final[[j]], file = sprintf("threshold_analysis/power_curve_output/sample_size_thresh%s.csv", j),
            row.names = FALSE)
}


# ------------------------------------------------------------------
# 02. PLOT POWER CURVES

# get Clopper-Pearson 95% CIs on empirical power numerator and denominator
# counts
df_CI <- df_sim_comb %>%
  mutate(CI_lower = qbeta(p = 0.025, shape1 = num, shape2 = denom - num + 1),
         CI_upper = qbeta(p = 0.975, shape1 = num + 1, shape2 = denom - num))

# make single long data.frame of all power curves
for (i in seq_along(list_post_summary)) {
  list_post_summary[[i]] <- df_filtered[i,] %>%
    bind_cols(list_post_summary[[i]])
}
df_curves <- list_post_summary %>%
  bind_rows() %>%
  filter(N <= N_max)

# produce separate plot for each prevalence threshold
for (j in 1:3) {
  
  # subset data to this prevalence threshold and subset N to within plotting range
  df_sub <- df_CI %>%
    filter(prev_thresh == sprintf("prev_thresh%s", j)) %>%
    left_join(df_filtered) %>%
    filter(N <= N_max)
  
  # produce grid plot of all power curves
  plot1 <- df_keep %>%
    left_join(df_curves) %>%
    filter(prev_thresh == sprintf("prev_thresh%s", j)) %>%
    ggplot() + theme_bw() +
    geom_ribbon(aes(x = N, ymin = lower, ymax = upper), fill = ribbon_col, alpha = ribbon_alpha) +
    geom_line(aes(x = N, y = mean), color = ribbon_col) +
    geom_errorbar(aes(x = N, ymin = CI_lower, ymax = CI_upper), width = error_width, data = df_sub) +
    geom_point(aes(x = N, y = num / denom), size = point_size, data = df_sub) +
    geom_hline(yintercept = target_power, linetype = "dotted") +
    geom_vline(aes(xintercept = N_fit), linetype = "dashed", data = df_sub) +
    ylim(c(0, 1)) +
    facet_wrap(vars(n_clust, p), nrow = 11, scales = "free_x") +
    theme(strip.background = element_rect(colour = "white", fill = "black"),
          strip.text = element_text(face = "bold", colour = "white"))
  
  # save plot to file
  ggsave(plot = plot1, filename = sprintf("threshold_analysis/power_curve_output/curve_thresh%s.pdf", j),
         width = 35, height = 35)
}

