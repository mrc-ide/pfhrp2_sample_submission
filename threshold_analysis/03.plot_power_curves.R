# 03.plot_power_curves.R
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

# load packages
library(tidyverse)

# load DRpower package from local version
# TODO - replace with tagged online version once ready
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# hill function
get_model_power <- function(A, B, K, w, N) {
  #A + (B - A) / (1 + K*N)^(w/K)
  A + (B - A)*(1 - exp(-K*N))^w
}

# ------------------------------------------------------------------
# make curves and save data to file

# load objects from file
df_sim <- readRDS("threshold_analysis/outputs/df_sim.rds") %>%
  #filter(ICC == 0.05) %>%
  filter(prior_ICC_shape2 == 1)# %>%
  #filter(prev_thresh == 0.05)

df_hill_fits <- readRDS("threshold_analysis/outputs/df_hill_fits.rds")


for (i in 1:nrow(df_hill_fits)) {
  message(sprintf("%s of %s", i, nrow(df_hill_fits)))
  
  # subset data
  df_sim_sub <- df_sim %>%
    filter(ICC == df_hill_fits$ICC[i]) %>%
    filter(prevalence == df_hill_fits$prevalence[i]) %>%
    filter(prev_thresh == df_hill_fits$prev_thresh[i]) %>%
    filter(n_clust == df_hill_fits$n_clust[i])
  
  df_power_sub <- df_hill_fits %>%
    filter(ICC == df_hill_fits$ICC[i]) %>%
    filter(prevalence == df_hill_fits$prevalence[i]) %>%
    filter(prev_thresh == df_hill_fits$prev_thresh[i]) %>%
    filter(n_clust == df_hill_fits$n_clust[i]) %>%
    expand_grid(N = 5:2e3) %>%
    mutate(power = get_model_power(A, B, K, w, N))
  
  # make plot
  df_power_sub %>%
    ggplot() + theme_bw() +
    geom_hline(yintercept = 0.8, col = "blue", alpha = 0.5) +
    geom_errorbar(aes(x = N, ymin = power_lower, ymax = power_upper), data = df_sim_sub) +
    geom_line(aes(x = N, y = power), col = "red") +
    ylim(c(0, 1)) + ylab("Power") +
    ggtitle(sprintf("prior_ICC_shape2 = %s, ICC = %s\nprevalence = %s, prev_thresh = %s\nn_clust = %s",
                    df_hill_fits$prior_ICC_shape2[i], df_hill_fits$ICC[i],
                    df_hill_fits$prevalence[i], df_hill_fits$prev_thresh[i],
                    df_hill_fits$n_clust[i])) + xlim(c(0, 300))
  
  # save to file
  filename <- sprintf("ICC_%s_prev_thresh_%s_prevalence_%s_n_clust_%s",
                      format(df_hill_fits$ICC[i], nsmall = 2),
                      format(df_hill_fits$prev_thresh[i], nsmall = 2),
                      format(df_hill_fits$prevalence[i], nsmall = 2),
                      sprintf("%02d", df_hill_fits$n_clust[i]))
  ggsave(sprintf("threshold_analysis/outputs/power_plots/%s.pdf", filename), 
         width = 8, height = 4)
}
