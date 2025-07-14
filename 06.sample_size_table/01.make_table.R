# 01.make_table.R
#
# Author: Bob Verity
# Date: 2025-07-14
#
# Inputs: none
#
# Outputs: none
#
# Purpose:
# Extracts and views sample size table from DRpower package. This script simply points to the correct version of DRpower to get the table - all the table simulation and calculation is done inside the DRpower package. Check out the R_ignore folder of that repos for fully reproducible workflow on how to produce this table.
#
# ------------------------------------------------------------------
library(tidyverse)
#devtools::install_github("mrc-ide/DRpower@v1.0.3")
library(DRpower)

# filter and reformat precomputed sample size table from DRpower package
df_ss |>
  filter(ICC == 0.05) |>
  filter(prev_thresh == 0.05) |>
  filter((prevalence >= 0.08) & (prevalence <= 0.151)) |>
  filter(n_clust <= 15) |>
  select(prevalence, n_clust, N_opt) |>
  pivot_wider(names_from = prevalence, values_from = N_opt)


