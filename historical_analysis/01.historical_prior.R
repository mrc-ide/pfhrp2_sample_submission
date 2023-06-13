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

# experiment with dropping Peru
#multi_site <- multi_site %>%
#  filter(WHO_country_name != "Peru")

# merge dataset and finalise
dat_final <- multi_site %>%
  left_join(dat1)

# loop through studies and get full posterior in each case
x <- seq(0, 1, l = 10001)
post_mat <- matrix(NA, nrow(multi_site), length(x))
l <- list()
for (i in seq_along(multi_site$study_ID)) {
  
  # subset to this study and get raw posterior ICC
  dat_sub <- dat_final %>%
    filter(study_ID == multi_site$study_ID[i])
  post_df <- get_ICC(n = dat_sub$Count_deletions, N = dat_sub$Tested,
                     prior_ICC_shape2 = 1, post_median_on = TRUE, post_full_on = TRUE,
                     n_intervals = 100, post_full_breaks = x)
  
  # save results
  l[[i]] <- post_df
  post_mat[i,] <- post_df$post_full[[1]]
}

# merge results with data
df_CrI <- multi_site %>%
  bind_cols(bind_rows(l)) %>%
  mutate(plot_name = paste(Paper_name, WHO_country_name, WHO_ADMIN1, sep = ": ")) %>%
  dplyr::select(plot_name, MAP, post_median, CrI_lower, CrI_upper)

# add row for combined results
post_combined <- apply(post_mat, 2, prod)
post_combined <- post_combined / sum(post_combined)
MAP_combined <- x[which.max(post_combined)]
post_median_combined <- x[which(cumsum(post_combined) > 0.5)[1]]
HDI_combined <- get_HDI_manual(x, post_combined)
df_CrI_combined <- rbind(df_CrI, data.frame(plot_name = "Combined",
                                            MAP = MAP_combined,
                                            post_median = post_median_combined,
                                            CrI_lower = HDI_combined["lower"],
                                            CrI_upper = HDI_combined["upper"]))

df_CrI_combined

plot(x, post_combined, type = 'l', xlim = c(0, 0.2))
abline(v = 0.04, lty = 2)

# save to file
saveRDS(df_CrI_combined, file = "historical_analysis/data/df_CrI_combined.rds")
