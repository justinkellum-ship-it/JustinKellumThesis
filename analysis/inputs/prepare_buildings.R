#!/usr/bin/env Rscript
# prepare_buildings.R ------------------------------------------------------------
# Makes data/external/buildings_oahu.gpkg from the FEMA / ORNL USA Structures inventory
# for Hawaiʻi (HI_Structures.gdb, https://gis-fema.hub.arcgis.com/pages/usa-structures):
#   * the footprints of Honolulu County (FIPS 15003: Oʻahu), 207,031 polygons;
#   * height_m = the inventory's HEIGHT attribute (metres; lidar-derived where present),
#     with HEIGHT_FILL_M (5 m) for the 46,926 footprints that have no recorded height;
#   * height_recorded = whether the inventory gives a height (FALSE for those 46,926;
#     1,263 other footprints have a recorded height of exactly 5 m, which the flag keeps apart);
#   * area_m2 = the polygon area in WGS 84 / UTM zone 4N (EPSG:32604), the CRS of the file.
# The analysis scripts clip heights to HEIGHT_RANGE_M when they read the file; the flag is
# used by check B of 07_sensitivity_checks.R.
#
# Usage: Rscript analysis/inputs/prepare_buildings.R <path to HI_Structures.gdb>
# The output is committed; run this only to rebuild it from the inventory.
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) file.path(dirname(.f[1]), "..") else ".." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE)

gdb <- commandArgs(trailingOnly = TRUE)[1]
stopifnot("give the path to HI_Structures.gdb" = !is.na(gdb) && dir.exists(gdb))
out <- file.path(EXT, BUILDINGS_FILE)

message("reading the inventory ...")
x <- st_read(gdb, query = "SELECT BUILD_ID, HEIGHT, FIPS FROM HI_Structures WHERE FIPS = '15003'", quiet = TRUE)
message(nrow(x), " footprints in Honolulu County; ", sum(is.na(x$HEIGHT)), " without a recorded height")
x <- st_transform(x, 32604)
b <- x %>% transmute(height_m = coalesce(HEIGHT, HEIGHT_FILL_M), height_recorded = !is.na(HEIGHT),
                     area_m2 = as.numeric(st_area(st_geometry(x))))
st_write(b, out, layer = "buildings_oahu", delete_dsn = file.exists(out), quiet = TRUE)
message("wrote ", out, ": ", nrow(b), " footprints")
