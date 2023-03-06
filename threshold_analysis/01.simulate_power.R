# 01.simulate_power.R
#
# Author: Bob Verity
# Date: 2023-03-03
#
# Purpose:
# Estimates empirical power by simulating data from the generative model (i.e.,
# the beta-binomial) and performing the exact analysis on simulated data that
# would be carried out on real data. Analysis is carried out via the
# get_prevalence() function, with the output being the probability of being
# above the set prevalence threshold. If this probability is above the rejection
# threshold (fixed at 0.95) then reject the null hypothesis. Repeat this process
# many times over a range of sample sizes and parameter combinations.
#
# Note that the get_prevalence() function can analyse the same data at multiple
# prevalence thresholds, and this is much more efficient that doing it
# separately. For this reason, a given parameter combination is always analysed
# at all possible thresholds, and results are stored in a list column within the
# final data.frame.
#
# ------------------------------------------------------------------

# load packages
library(dplyr)

# load DRpower package from local version
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# define or load parameter combinations
if (FALSE) {
  df_params <- expand_grid(rho = 0.05,
                           n_clust = c(2:10, 15, 20),
                           p = 6:20 / 100,
                           lower = 5,
                           upper = 2e3) %>%
    mutate(param_i = seq_along(rho))
  saveRDS(df_params, file = "threshold_analysis/data/df_params.rds")
}
if (FALSE) {
  df_params <- readRDS("threshold_analysis/data/df_params.rds")
}

# define fixed parameters
n_seq <- 11
reps <- 1e3
p_thresh <- c(0.05, 0.08, 0.10)
prior_ICC_shape2 <- 19

# loop through parameter combinations
t0 <- Sys.time()
res_list <- list()
for (param_i in 1:nrow(df_params)) {
  
  # skip if not within max
  if (is.na(df_params$upper[param_i])) {
    next
  }
  
  # write progress to console
  message(sprintf("param = %s/%s, rho = %s, n_clust = %s, p = %s", param_i, nrow(df_params),
                  df_params$rho[param_i], df_params$n_clust[param_i], df_params$p[param_i]))
  
  # define sequence of N values to explore
  N_lower <- df_params$lower[param_i]
  N_upper <- df_params$upper[param_i]
  N_seq <- unique(round(seq(N_lower, N_upper, l = n_seq)))
  
  # initialise dataframe for storing simulation results of this param_i
  df_sim <- data.frame(param_i = param_i,
                       N = N_seq)
  
  # explore range of N
  power_num <- list()
  for (j in seq_along(N_seq)) {
    
    # simulate repeatedly from lower limit
    power_num[[j]]  <- rep(0, 3)
    for (i in 1:reps) {
      n <- rbbinom_reparam(n = df_params$n_clust[param_i],
                           N = N_seq[j],
                           p = df_params$p[param_i],
                           rho = df_params$rho[param_i])
      
      for (int_i in seq(10, 100, 10)) {
        p_est <- try(get_prevalence(n = n, N = N_seq[j],
                                    prev_thresh = p_thresh,
                                    prior_prev_shape1 = 1, prior_prev_shape2 = 1,
                                    prior_ICC_shape1 = 1, prior_ICC_shape2 = prior_ICC_shape2,
                                    return_type = list(mean_on = FALSE,
                                                       median_on = FALSE,
                                                       CrI_on = FALSE,
                                                       thresh_on = TRUE,
                                                       full_on = FALSE),
                                    n_intervals = int_i, debug_on = FALSE,
                                    use_cpp = TRUE), silent = T)
        if (class(p_est) != "try-error") {
          break()
        }
      }
      
      # add to running numerator
      power_num[[j]] <- power_num[[j]] + (p_est$prob_above_threshold[[1]] > 0.95)
    }
  }
  
  # store results
  df_sim$num <- power_num
  df_sim$denom <- reps
  
  res_list[[param_i]] <- df_sim
}
print(Sys.time() - t0)

# convert to data.frame
res_df <- res_list %>%
  bind_rows()

# merge back with parameters and drop NAs
df_loess <- df_params %>%
  left_join(res_df) %>%
  filter(!is.na(N))

# save to file
if (FALSE) {
  saveRDS(df_loess, file = "threshold_analysis/data/df_sim_04.rds")
}
