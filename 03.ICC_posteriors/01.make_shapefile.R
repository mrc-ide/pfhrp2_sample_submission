# 01.make_shapefile.R
#
# Author: Bob Verity
# Date: 2025-07-11
#
# Inputs: Multiple shapefiles from GADM in the GADM_shapefiles directory
#
# Outputs: 03.ICC_posteriors/data/shp_combined.rds
#
# Purpose:
# Reads in all shapefiles, combines and tidies up before saving as a single object
#
# ------------------------------------------------------------------
library(tidyverse)
library(sf)

# Define the parent directory
parent_dir <- "03.ICC_posteriors/data/GADM_shapefiles"

# Get all subdirectories within the parent directory
subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)

# Read all shapefiles into a named list
shapefiles <- lapply(subdirs, function(dir) {
  st_read(dir, quiet = TRUE)
}) |>
  bind_rows()

# Fix invalid geometries
shapefiles <- st_make_valid(shapefiles)

# save shapefile
saveRDS(shapefiles, file = "03.ICC_posteriors/data/shp_combined.rds")

