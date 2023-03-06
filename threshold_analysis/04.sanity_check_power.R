# 04.sanity_check_power.R
#
# Author: Bob Verity
# Date: 2023-03-06
#
# Purpose:
# Performs a sanity check on the sample sizes calculated in 03.plot_power.R to
# ensure power does indeed hit the target.
#
# ------------------------------------------------------------------

# define parameters
n_clust <- 10
N <- 30
p <- 0.10
reps <- 1e2

# get empirical power
t0 <- Sys.time()
get_power(N = rep(N, n_clust),
          prevalence = p,
          ICC = 0.05,
          prior_ICC_shape2 = 19,
          reps = reps)
Sys.time() - t0
