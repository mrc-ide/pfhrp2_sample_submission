# sanity_check_power.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Performs a sanity check on the sample sizes calculated in
# produce_presence_table.R to ensure power does indeed hit the target. Does this
# via simulation.
#
# ------------------------------------------------------------------

# load packages
library(extraDistr)
library(dplyr)

# define parameters
N <- 8
n_clust <- 5
p <- 0.05
rho <- 0.05
reps <- 1e4

# simulate data from beta-binomial model
alpha <- p*(1 / rho - 1)
beta <- (1 - p)*(1 / rho - 1)
sim_post <- rbbinom(n_clust*reps, N, alpha, beta) %>%
  matrix(ncol  = n_clust)

# count proportion of times at least one positive result over all clusters
mean(rowSums(sim_post) > 0)
