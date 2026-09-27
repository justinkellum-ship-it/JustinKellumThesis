# helpers.R ------------------------------------------------------------------
# Shared paths, constants and small functions for the thesis analysis in R.
#
# Every script in this folder starts with `source("helpers.R")` (relative to
# the analysis/ folder).  The project root is the folder that contains analysis/,
# data/ and results/; it is found from the location of this file, or can be forced
# with the environment variable THESIS_ROOT (used by run_all.R and when the scripts
# are run from the R Markdown walkthrough).
#
# Packages used here:
#   dplyr, tidyr, tibble, readr, purrr, stringr  (tidyverse data handling)
#   lubridate                                    (rounding of time stamps)
#   jsonlite                                     (results file)
# ---------------------------------------------------------------------------

# the scripts contain UTF-8 text (the Hawaiian okina, degree signs, Greek letters); make
# sure the session reads and draws them as such even under a plain "C" locale.  This file
# itself is plain ASCII so that it can be sourced before the locale is set.
if (!isTRUE(l10n_info()[["UTF-8"]])) suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
options(encoding = "UTF-8")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(purrr)
  library(stringr)
  library(lubridate)
  library(jsonlite)
})

# ---- project root -----------------------------------------------------------
.find_root <- function() {
  env <- Sys.getenv("THESIS_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env))
  is_root <- function(d) dir.exists(file.path(d, "analysis")) && dir.exists(file.path(d, "data"))
  # script run with Rscript: --file=path/to/analysis/xx.R or analysis/inputs/xx.R; walk up
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f) == 1) {
    d <- normalizePath(dirname(f))
    for (i in 1:3) { if (is_root(d)) return(d); d <- dirname(d) }
  }
  # sourced interactively from the analysis/ folder or from the project root
  if (file.exists("helpers.R") && is_root("..")) return(normalizePath(".."))
  if (is_root(".")) return(normalizePath("."))
  stop("Cannot find the project root: set THESIS_ROOT or run from analysis/ or the project folder")
}
ROOT <- .find_root()
DATA <- file.path(ROOT, "data")
RAW  <- file.path(DATA, "raw")          # the study's own measurements (logger readings, sites)
EXT  <- file.path(DATA, "external")     # data from other providers (weather series, shoreline, DEM, C-CAP, FEMA ...)
PROC <- file.path(DATA, "processed")    # built by 00_prepare_inputs.R
DER  <- file.path(DATA, "descriptors")  # built by 01_site_predictors.R and the checks (site_elevation.csv: inputs/fetch_site_elevation.py)
OUT  <- file.path(ROOT, "results")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
PY   <- file.path(ROOT, "replication", "results")   # the independent Python implementation writes here
for (d in c(PROC, DER, TAB, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- study configuration ----------------------------------------------------
# Everything that ties the pipeline to this particular network is set here (see
# AGENTS.md for how to set it for another one).  The scripts read these values;
# none of them repeats a number that lives here.  The alternative values tried by
# the sensitivity steps (other thresholds, the 3 m building height, the calm-core
# nights) are the checks' own parameters and stay in those scripts.
BIN_MINUTES <- 30        # logging interval: every reading is snapped to this grid (15, 30, 60 ...)
NIGHT_START <- 18        # local hour at which a night begins (HST; local standard time, no DST)
NIGHT_END   <- 6         # local hour at which it ends; a diurnal date runs 06:00-06:00
NIGHT_HOURS <- (NIGHT_END - NIGHT_START) %% 24              # 12
NIGHT_BINS  <- as.integer(NIGHT_HOURS * 60 / BIN_MINUTES)    # 24 bins at 30 min
WIND_MAX_KMH  <- 10      # 10 m wind speed threshold (weather series) for a "calm" bin
CLOUD_MAX_PCT <- 25      # total cloud cover threshold (weather series) for a "clear" bin
MIN_FRAC_STEPS <- 0.75   # share of the night bins that must be calm AND clear
MIN_SENSORS_PER_NIGHT <- 20      # network completeness for a usable reference
MIN_FRAC_SENSOR_NIGHT <- 0.75    # sensor-night completeness: share of the night bins present ...
MIN_STEPS_PER_SENSOR_NIGHT <- as.integer(round(MIN_FRAC_SENSOR_NIGHT * NIGHT_BINS))   # ... 18 of 24 bins

# the districts, in the order they are reported: read from the site table when it exists
REGIONS <- if (file.exists(file.path(RAW, "sites.csv"))) {
  unique(readr::read_csv(file.path(RAW, "sites.csv"), show_col_types = FALSE, col_types = readr::cols(.default = "c"))$region)
} else c("Honolulu", "Ewa")
REGION_TERM <- paste0("region", REGIONS[1])   # the district coefficient: the first district against the second (the reference level)
RADII   <- c(50, 100, 200)      # source-area radii (m); 02_models.R adopts one of them by AIC
PARSIMONIOUS <- c("imperv", "height")   # the reduced pooled model (district + these descriptors) that the sensitivity steps refit

# geospatial inputs (data/external): a projected CRS in metres, the land-cover class
# rasters (one set of tiles per class, matched by a pattern in the file name, with the
# pixel value that marks the class), the building footprints (a layer with a height_m
# column; buildings without a height get HEIGHT_FILL_M) and the coastline polygon
CRS_M <- 6634                    # NAD83(PA11) / UTM zone 4N: the projection of the C-CAP tiles
LANDCOVER_DIR <- "ccap"
LANDCOVER <- list(impervious = list(pattern = "impervious", value = 1),
                  tree       = list(pattern = "canopy",     value = 1),   # C-CAP canopy mask: 1 = tree, 2 = shrub
                  water      = list(pattern = "water",      value = 1))
BUILDINGS_FILE <- "buildings_oahu.gpkg"
HEIGHT_FILL_M  <- 5              # height given to footprints without a recorded height
HEIGHT_RANGE_M <- c(0.5, 60)     # heights are clipped to this range
COAST_FILE     <- "oahu_gshhs_f.geojson"
PRED    <- c("imperv", "tree", "bldg", "water", "height", "svf_point", "aspect", "coast_km")
LABEL   <- c(imperv = "Impervious surface fraction (%)", tree = "Tree canopy fraction (%)",
             bldg = "Building footprint fraction (%)", water = "Water surface fraction (%)",
             height = "Mean building height (m)", svf_point = "Sky view factor (building, at sensor)",
             aspect = "Canyon aspect ratio (H/W)", coast_km = "Distance to coast (km)")
SHORT   <- c(imperv = "impervious", tree = "tree canopy", bldg = "footprint", water = "water",
             height = "height", svf_point = "sky view", aspect = "aspect ratio", coast_km = "coast")
# model formulas (right-hand sides, without the random intercept that fit_lmm adds)
rhs_full        <- function(extra = "region") paste(c(extra, paste0(PRED, "_z")), collapse = " + ")
rhs_parsimonious <- function(extra = "region") paste(c(extra, paste0(PARSIMONIOUS, "_z")), collapse = " + ")
# natural-unit conversion of z-scored slopes: (multiplier, label)
UNIT <- list(imperv = list(10, "per +10 pp"), tree = list(10, "per +10 pp"), bldg = list(10, "per +10 pp"),
             water = list(10, "per +10 pp"), height = list(1, "per +1 m"), svf_point = list(-0.1, "per -0.1"),
             aspect = list(0.1, "per +0.1"), coast_km = list(1, "per +1 km"))

# ---- small utilities --------------------------------------------------------
# time stamps in the record are local standard time (HST, no daylight saving);
# they are handled as UTC-labelled clock times so that no offset is ever applied
read_stamps <- function(x) {
  # readr may already have parsed the column; a midnight stamp printed through
  # as.character() would lose its time part, so POSIXct input is only re-labelled
  if (inherits(x, "POSIXct")) return(force_tz(with_tz(x, "UTC"), "UTC"))
  ymd_hms(x, tz = "UTC", quiet = TRUE)
}
fmt_stamp   <- function(x) format(x, "%Y-%m-%d %H:%M:%S")

write_table <- function(df, name, dir = TAB) {
  readr::write_csv(df, file.path(dir, paste0(name, ".csv")), na = "")
  invisible(df)
}

zscore_frame <- function(df, cols, ref = NULL) {
  # z-scores with the mean and SD of `ref` (default: df itself), so that
  # coefficients of subsets stay comparable with the pooled model
  if (is.null(ref)) ref <- df
  sds <- means <- setNames(numeric(length(cols)), cols)
  for (cc in cols) {
    mu <- mean(ref[[cc]]); s <- sd(ref[[cc]])
    df[[paste0(cc, "_z")]] <- (df[[cc]] - mu) / s
    sds[cc] <- s; means[cc] <- mu
  }
  list(df = df, sds = sds, means = means)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
