# 01.make_posterior.R
#
# Author: Bob Verity
# Date: 2024-06-10
#
# Inputs: (none)
#
# Outputs: outputs/posterior_plot.png
#
# Purpose:
# Creates and example posterior distribution plot from DRpower.
#
# ------------------------------------------------------------------
library(tidyverse)

# install and load DRpower
#devtools::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)

# make plot
n_clust <- 4
n <- rep(3, n_clust)
N <- rep(100, n_clust)
plot_prevalence(n = n,
                N = N,
                prev_range = c(0, 0.25))

# save
ggsave("02.posterior_example/outputs/posterior_plot.png", scale = 0.7)
