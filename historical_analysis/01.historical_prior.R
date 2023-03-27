# 01.historical_prior.R
#
# Author: Bob Verity
# Date: 2023-02-28
#
# Purpose:
# Reads in historical data on pfhrp2/3 deletions and analyses data to arrive at
# a suitable prior on the ICC for future studies. Writes credible intervals on
# the ICC to file.
#
# ------------------------------------------------------------------

# takes a discritised version of the posterior and computes the high density
# interval by brute force. Only works on unimodal distributions
get_HDI_manual <- function(x, y, alpha = 0.05) {
  
  # get rough area of each interval via trapezoidal rule and normalise area under
  # curve
  post_area <- (y[-1] + y[-length(y)]) / 2
  post_area <- post_area / sum(post_area)
  
  # find the y-threshold corresponding to the HDI at threshold alpha
  area_sort <- sort(post_area, decreasing = TRUE)
  area_cs <- cumsum(area_sort)
  w <- which(area_cs > (1 - alpha))[1]
  yt <- area_sort[w]
  
  # find the range of x-values implied by this threshold
  w2 <- which(post_area > yt)
  ret <- c(lower = min(x[w2]),
           upper = max(x[-1][w2]))
  return(ret)
}

# other packages
library(dplyr)

# load working version of DRpower
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Sample design/DRpower")

# ------------------------------------------------------------------

# read in data
dat1 <- data.table::fread("historical_analysis/data/included_data.csv")

# subset to ADMIN1 regions that have data for 3 or more lower geographical units
# in the same year
multi_site <- dat1 %>%
  group_by(Paper_name, WHO_country_name, WHO_ADMIN1, Year) %>%
  summarise(unique_count = n()) %>%
  ungroup() %>%
  filter(unique_count >= 3) %>%
  mutate(study_ID = seq_along(unique_count))

# merge dataset and finalise
dat_final <- multi_site %>%
  left_join(dat1)

# get final IDs
study_IDs <- unique(dat_final$study_ID)

# loop through studies and get full posterior in each case
x <- seq(0, 1, l = 1001)
post_mat <- matrix(NA, length(study_IDs), length(x))
for (i in seq_along(study_IDs)) {
  
  # subset to this study and get raw posterior ICC
  dat_sub <- dat_final %>%
    filter(study_ID == study_IDs[i])
  post_df <- get_ICC(n = dat_sub$Count_deletions, N = dat_sub$Tested,
                     return_type = list(mean_on = FALSE, median_on = FALSE, CrI_on = FALSE, full_on = TRUE),
                     n_intervals = 100)
  post_mat[i,] <- post_df$full[[1]]
  
}

# fix numerical errors
post_mat[post_mat < 0] <- 0

# get the posterior median and CrI for all studies
df_CrI <- multi_site %>%
  mutate(plot_name = paste(Paper_name, WHO_country_name, WHO_ADMIN1, sep = ": "),
         median = NA, lower = NA, upper = NA) %>%
  dplyr::select(plot_name, median, lower, upper)
for (i in seq_along(study_IDs)) {
  y <- post_mat[i,] / sum(post_mat[i,])
  df_CrI$median[i] <- x[which(cumsum(y) > 0.5)[1]]
  HDI <- get_HDI_manual(x, post_mat[i,])
  df_CrI$lower[i] <- HDI["lower"]
  df_CrI$upper[i] <- HDI["upper"]
}

# get the product of the posterior over all studies
post_combined <- apply(post_mat, 2, prod)
post_median <- x[which(cumsum(post_combined / sum(post_combined)) > 0.5)[1]]
HDI_combined <- get_HDI_manual(x, post_combined)
df_CrI_2 <- rbind(df_CrI, data.frame(plot_name = "Combined",
                                     median = post_median,
                                     lower = HDI_combined["lower"],
                                     upper = HDI_combined["upper"]))

# quick plot of posteriors
plot(x, x, type = 'n', ylim = c(0, 0.025), xlab = "ICC", ylab = "prob")
for (i in seq_len(nrow(post_mat))) {
  lines(x, post_mat[i,] / sum(post_mat[i,]), col = "#00000050")
}
lines(x, post_combined / sum(post_combined), col = 2)

# save to file
saveRDS(df_CrI_2, file = "historical_analysis/data/df_CrI.rds")
