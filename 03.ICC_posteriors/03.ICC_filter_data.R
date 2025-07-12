# 03.ICC_filter_data.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs:
# - 03.ICC_posteriors/outputs/dat_filter1.rds
#
# Outputs:
# - 03.ICC_posteriors/outputs/dat_ICCfilter.rds
#
# Purpose:
# Reads in pre-filtered data. Applies additional filters needed to ensure robust ICC estimation. Saves filtered data to file for plotting.
#
# ------------------------------------------------------------------
library(tidyverse)

# read in data produced by previous script
dat <- readRDS("03.ICC_posteriors/outputs/dat_filter1.rds")

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# must have at least 10 samples per site
dat <- dat |>
  filter(HRP2_TESTED >= 10)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# must have at least one pfhrp2 deletion per site
dat <- dat |>
  filter(HRP2_NUM_DELETION > 0)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# count sites per ADMIN1, run in the same year
dat_nsites <- dat |>
  group_by(COUNTRY_NAME, ADMIN1_NAME, YEAR_START) |>
  summarise(NSITES_COUNTRY = n()) |>
  ungroup()

# must have at least 3 sites per ADMIN1
dat <- dat |>
  left_join(dat_nsites) |>
  filter(NSITES_COUNTRY >= 3) |>
  select(-NSITES_COUNTRY)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# write to file
saveRDS(dat, file = "03.ICC_posteriors/outputs/dat_ICCfilter.rds")

