#!/usr/bin/env Rscript
# extract_island_polygon.R -------------------------------------------------------
# Makes data/external/oahu_gshhs_f.geojson: the Oʻahu land polygon from the GSHHG
# full-resolution shoreline, level 1 (land), Wessel & Smith (1996), downloaded as the
# shapefile distribution from https://www.soest.hawaii.edu/pwessel/gshhg/ and unpacked
# into data/external/gshhg/ (GSHHS_shp/f/GSHHS_f_L1.shp).  The polygon is the one that
# contains the centre of Honolulu; it is kept in WGS 84 with the shapefile's attributes.
# Distance to the coast in 01_site_predictors.R is measured to the boundary of this
# polygon, so Pearl Harbor counts as coast.
#
# Usage: Rscript analysis/inputs/extract_island_polygon.R
# The output is committed; run this only to rebuild it from the shoreline data.
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) file.path(dirname(.f[1]), "..") else ".." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")
suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE)

shp <- file.path(EXT, "gshhg", "GSHHS_shp", "f", "GSHHS_f_L1.shp")
stopifnot("GSHHS_f_L1.shp not found under data/external/gshhg" = file.exists(shp))
honolulu <- st_sfc(st_point(c(-157.858, 21.307)), crs = 4326)        # a point in the Honolulu core
# read only the polygons of the Oʻahu window, then keep the one that contains the point
land <- st_read(shp, wkt_filter = st_as_text(st_as_sfc(st_bbox(c(xmin = -158.4, ymin = 21.2, xmax = -157.6, ymax = 21.8), crs = st_crs(4326)))), quiet = TRUE)
oahu <- land[lengths(st_intersects(land, honolulu)) > 0, ]
stopifnot(nrow(oahu) == 1)
out <- file.path(EXT, COAST_FILE)
if (file.exists(out)) invisible(file.remove(out))
st_write(oahu, out, driver = "GeoJSON", quiet = TRUE)
message("wrote ", out, ": one polygon, area ", round(as.numeric(st_area(st_transform(oahu, CRS_M))) / 1e6), " km²")
