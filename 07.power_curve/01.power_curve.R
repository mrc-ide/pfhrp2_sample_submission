# 01.power_curve.R
#
# Author: Bob Verity
# Date: 2025-07-14
#
# Inputs: none
#
# Outputs: none
#
# Purpose:
# Extracts and views power estimates from DRpower package. Produces an example power curve.
#
# This script simply points to the correct version of DRpower to get the data frame of results - all the simulation and calculation is done inside the DRpower package. Check out the R_ignore folder of that repos for fully reproducible workflow on how to produce this table.
#
# ------------------------------------------------------------------
library(tidyverse)
#devtools::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)
library(grid)  # for arrow()

# choose filtering parameters
filter_n_clust <- 8

# filter simulation results
df_plot <- df_sim |>
  filter(ICC == 0.05) |>
  filter(prev_thresh == 0.05) |>
  filter((prevalence > 0.099) & (prevalence < 0.1001)) |>
  filter(n_clust == filter_n_clust)

# get chosen sample size from other precomputed table (calculated using linear interpolation)
ss <- df_ss |>
  filter(ICC == 0.05) |>
  filter(prev_thresh == 0.05) |>
  filter((prevalence > 0.099) & (prevalence < 0.1001)) |>
  filter(n_clust == filter_n_clust) |>
  pull(N_opt)

# plot curve
df_plot |>
  ggplot() + theme_bw() +
  geom_pointrange(aes(x = N, y = power, ymin = lower, ymax = upper)) +
  geom_line(aes(x = N, y = power)) +
  geom_hline(yintercept = 80, linetype = "dashed") +
  geom_segment(aes(x = ss, xend = ss, y = 80, yend = 0),
               arrow = arrow(length = unit(0.2, "inches")),
               color = "red", linewidth = 0.5) +
  annotate("text",
           x = ss + 2, y = 10,  # adjust position as needed
           label = ss,
           color = "red", size = 4, hjust = 0) +
  xlim(c(0, 100)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ggtitle(sprintf("Power Curve for %s Clusters", filter_n_clust)) +
  xlab("Sample size (per cluster)") + ylab("Power")

# save to file
ggsave("07.power_curve/outputs/power_curve.png")
