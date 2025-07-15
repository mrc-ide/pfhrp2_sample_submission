# 02.prefilter_data.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs:
# - 03.ICC_posteriors/data/MTM_PFHRP23_GENE_DELETIONS_20231127_edited.xlsx, data from WHO Malaria Threats Map
# - 03.ICC_posteriors/data/continent_key.csv, a key mapping countries to continents
# - 03.ICC_posteriors/data/shp_combined.rds, shapefile including all countries in the data
#
# Outputs:
# - 03.ICC_posteriors/outputs/dat_filter1.rds
#
# Purpose:
# Creates and example posterior distribution plot from DRpower.
#
# ------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(sf)

# ------------------------------------------------------------------
# BASIC FILTERS

# read in raw data
#dat <- openxlsx::read.xlsx("03.ICC_posteriors/data/MTM_PFHRP23_GENE_DELETIONS_20231127_edited.xlsx", sheet = 2)
dat <- openxlsx::read.xlsx("03.ICC_posteriors/data/MTM_PFHRP23_GENE_DELETIONS_20250714_edited.xlsx", sheet = 2)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# filter by continent
continent_key <- read.csv("03.ICC_posteriors/data/continent_key.csv")
dat <- dat |>
  left_join(continent_key) |>
  filter(CONTINENT_NAME %in% c("Africa", "Asia", "South America"))

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# filter by symptomatic
dat <- dat |>
  filter(PATIENT_TYPE == "Symptomatic")

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# filter by survey type
dat <- dat |>
  filter(SURVEY_TYPE %in% c("Convenience survey", "Cross-sectional prospective survey"))

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# ------------------------------------------------

# count how many studies and sites need to be discarded
dat_discard <- filter(dat, discard)
#View(dat_discard) # look for reasons for discard

length(unique(dat_discard$study_ID)) # count unique studies
nrow(distinct(dat_discard, study_ID, SITE_NAME)) # count unique sites

# filter out discarded studies
dat <- filter(dat, !discard)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# clean up columns
dat <- dat |>
  mutate(LATITUDE = as.numeric(LATITUDE),
         LONGITUDE = as.numeric(LONGITUDE),
         HRP2_TESTED = as.numeric(HRP2_TESTED),
         HRP2_PROPORTION_DELETION = as.numeric(HRP2_PROPORTION_DELETION),
         HRP2_NUM_DELETION = round(HRP2_TESTED * HRP2_PROPORTION_DELETION)) |>
  select(CONTINENT_NAME, COUNTRY_NAME, SITE_NAME, LONGITUDE, LATITUDE, YEAR_START,
         YEAR_END, HRP2_TESTED, HRP2_NUM_DELETION, CITATION_URL, study_ID)

nrow(dat)

# combine data collected in exact same location (lat/lon) in same year
dat <- dat |>
  group_by(CONTINENT_NAME, COUNTRY_NAME, SITE_NAME, LONGITUDE, LATITUDE, YEAR_START,
           YEAR_END, CITATION_URL, study_ID) |>
  summarise(HRP2_TESTED = sum(HRP2_TESTED),
            HRP2_NUM_DELETION = sum(HRP2_NUM_DELETION)) |>
  ungroup()

# the number of rows in the data decreases, but the number of unique studies and
# sites remains the same (by definition, because site names are identical)
nrow(dat)
length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# read in sf object containing all the relevant countries
shp_combined <- readRDS("03.ICC_posteriors/data/shp_combined.rds")

# make an sf version of the site coordinates
dat_sf <- data.frame(
  x = dat$LONGITUDE,
  y = dat$LATITUDE) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 4326)

# find which polygon each point intersects
intersect_vec <- mapply(function(x) x[1], st_intersects(dat_sf, shp_combined))

# there should be no NA values
any(is.na(intersect_vec))

# merge the ADMIN1 unit back with the original data
dat <- dat |>
  mutate(COUNTRY = shp_combined$COUNTRY[intersect_vec],
         ADMIN1_NAME = shp_combined$NAME_1[intersect_vec]) |>
  select(CONTINENT_NAME, COUNTRY_NAME, ADMIN1_NAME,SITE_NAME, LONGITUDE, LATITUDE, 
         YEAR_START, YEAR_END, HRP2_TESTED, HRP2_NUM_DELETION, CITATION_URL, study_ID)

# merge with manually added entries
dat_add <- read.csv("03.ICC_posteriors/data/additional_data.csv")

dat <- bind_rows(dat, dat_add)

length(unique(dat$study_ID)) # count unique studies
nrow(distinct(dat, study_ID, SITE_NAME)) # count unique sites

# save data to file
saveRDS(dat, file = "03.ICC_posteriors/outputs/dat_filter1.rds")
