# 01.simulate_power.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Inputs:
# None
#
# Outputs:
# df_sim.rds: estimates of empirical power and other summaries from repeated
# simulation over a wide range of parameter combinations
#
# Purpose:
# Produce binomial estimates of power (number of times correctly rejected null
# hypothesis vs. total number of trials) for a range of parameter combinations
#
# ------------------------------------------------------------------

# load packages
library(dplyr)

# load DRpower package from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# TODO - copy over code that was run on cluster

# save to file
if (FALSE) {
  saveRDS(df_sim, file = "threshold_analysis/outputs/df_sim.rds")
}
