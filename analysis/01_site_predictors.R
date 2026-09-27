#!/usr/bin/env Rscript
# 01_site_predictors.R ---------------------------------------------------------
# Two- and three-dimensional urban-form descriptors for every logger site at
# three source-area radii (50, 100, 200 m), from public inputs only:
#
#   NOAA C-CAP 2021 1 m land-cover masks (Hawaiʻi tiles, EPSG:6634)
#       -> impervious, tree-canopy and water fractions of the circle
#   FEMA / ORNL USA Structures footprints with heights (Oʻahu subset)
#       -> building footprint fraction, area-weighted mean building height,
#          canyon aspect ratio H/W from the mean building spacing
#   3-D ray casting on the same footprints
#       -> sky view factor at the sensor (1.5 m above ground), 36 azimuths,
#          200 m horizon search; SVF = mean cos^2(beta) (Dozier & Frew 1990)
#   GSHHG full-resolution shoreline (Wessel & Smith 1996)
#       -> distance from the sensor to the ocean
#
# Output: data/descriptors/site_predictors_multiscale.csv, one row per site x radius.
#
# Packages: sf (vector geometry: st_transform, st_buffer, st_intersects,
#           st_intersection, st_distance), terra (raster access: vrt, extract),
#           dplyr/purrr (bookkeeping).  All geometry is handled in NAD83(PA11) /
#           UTM zone 4N (EPSG:6634), the projection of the C-CAP tiles, so that
#           distances and areas are in metres.
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(sf); library(terra) })
Sys.setenv(PROJ_NETWORK = "OFF"); sf_proj_network(FALSE)   # no datum-grid downloads
sf_use_s2(FALSE)

CRS_M    <- 6634        # NAD83(PA11) / UTM zone 4N (metres) - the C-CAP projection
SENSOR_Z <- 1.5         # sensor height above ground (m)
N_AZ     <- 36          # azimuths for the sky view factor
HORIZON_R <- 200        # horizon search radius for the SVF (m)
QUADSEGS <- 16          # 64-vertex circles (same discretisation as the Python version)

# ---- inputs -------------------------------------------------------------------
ccap_mask <- function(kind) {
  # one virtual mosaic per C-CAP mask; the tiles do not overlap
  files <- list.files(file.path(EXT, "ccap"), pattern = paste0(kind, ".*\\.tif$"), full.names = TRUE)
  vrt(files, filename = file.path(tempdir(), paste0("ccap_", kind, ".vrt")), overwrite = TRUE)
}

raster_fraction <- function(r, geom, values) {
  # share (%) of the 1 m pixels whose centre lies inside `geom` that equal each of `values`
  # (pixels with the file's no-data value 0 are read as NA and count towards the total)
  v <- terra::extract(r, vect(geom), touches = FALSE)[[2]]
  n <- length(v)
  if (n == 0) return(list(frac = setNames(rep(NA_real_, length(values)), values), n = 0L))
  list(frac = setNames(vapply(values, function(k) 100 * sum(v == k, na.rm = TRUE) / n, numeric(1)), values), n = n)
}

building_metrics <- function(geom, bldg) {
  # descriptors of the footprints intersecting the circle `geom` (projected metres)
  area <- as.numeric(st_area(geom))
  idx <- st_intersects(geom, bldg)[[1]]
  empty <- list(bldg_frac = 0, zH = 0, zH_max = 0, n_bldg = 0L, hw_ratio = 0)
  if (length(idx) == 0) return(empty)
  sub <- bldg[idx, ]
  clipped <- suppressWarnings(st_intersection(st_geometry(sub), geom))
  a <- as.numeric(st_area(clipped))
  keep <- a > 0
  a <- a[keep]; h <- sub$height_m[keep]
  if (!length(a)) return(empty)
  bldg_area <- sum(a)
  zH <- sum(a * h) / bldg_area
  n <- length(a)
  # canyon aspect ratio H/W from the mean building spacing (Stewart & Oke 2012
  # style approximation)
  spacing <- sqrt(area / n) - sqrt(bldg_area / n)
  W <- max(6, spacing)
  list(bldg_frac = 100 * bldg_area / area, zH = zH, zH_max = max(h), n_bldg = n, hw_ratio = min(5, zH / W))
}

sky_view_factor <- function(pt, bldg, z0 = SENSOR_Z, R = HORIZON_R, n_az = N_AZ) {
  # SVF of a horizontal surface at height z0: mean over n_az azimuths of cos^2(beta),
  # beta = elevation angle of the highest building obstruction along the azimuth
  idx <- st_intersects(st_buffer(pt, R, nQuadSegs = QUADSEGS), bldg)[[1]]
  if (length(idx) == 0) return(1)
  cand <- bldg[idx, ]
  xy <- st_coordinates(pt)
  total <- 0
  for (k in seq_len(n_az) - 1) {
    th <- 2 * pi * k / n_az
    ray <- st_sfc(st_linestring(rbind(xy, xy + R * c(sin(th), cos(th)))), crs = st_crs(pt))
    hit <- cand[st_intersects(ray, cand)[[1]], ]
    beta <- 0
    if (nrow(hit)) {
      for (j in seq_len(nrow(hit))) {
        seg <- suppressWarnings(st_intersection(ray, st_geometry(hit)[j]))
        d <- max(1, as.numeric(st_distance(pt, seg)))
        beta <- max(beta, atan(max(hit$height_m[j] - z0, 0) / d))
      }
    }
    total <- total + cos(beta)^2
  }
  total / n_az
}

main <- function() {
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), show_col_types = FALSE)
  sites <- sn %>% distinct(sensor_id, .keep_all = TRUE) %>% select(sensor_id, region, latitude, longitude) %>%
    arrange(region, sensor_id)
  pts <- st_transform(st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)

  message("Loading building footprints ...")
  bldg <- st_transform(st_read(file.path(EXT, "buildings_oahu.gpkg"), quiet = TRUE), CRS_M) %>%
    mutate(height_m = pmin(pmax(coalesce(height_m, 3), 0.5), 60))
  # keep only the footprints that can matter (within 250 m of any site: the largest
  # radius plus a margin) so that the repeated spatial queries below stay fast
  near <- unique(unlist(st_intersects(st_buffer(st_geometry(pts), max(RADII, HORIZON_R) + 50), bldg)))
  bldg <- bldg[sort(near), ]
  message(nrow(bldg), " footprints within reach of the sites")
  coast <- st_transform(st_read(file.path(EXT, "oahu_gshhs_f.geojson"), quiet = TRUE), CRS_M)
  coast_line <- st_boundary(st_union(coast))
  masks <- list(impervious = ccap_mask("impervious"), canopy = ccap_mask("canopy"), water = ccap_mask("water"))

  rows <- list()
  for (i in seq_len(nrow(pts))) {
    p <- pts[i, ]
    svf <- sky_view_factor(st_geometry(p), bldg)
    coast_d <- as.numeric(st_distance(st_geometry(p), coast_line))
    for (R in RADII) {
      g <- st_buffer(st_geometry(p), R, nQuadSegs = QUADSEGS)
      imp <- raster_fraction(masks$impervious, g, 1)$frac[["1"]]
      can <- raster_fraction(masks$canopy, g, c(1, 2))$frac
      wat <- raster_fraction(masks$water, g, 1)$frac[["1"]]
      bm <- building_metrics(g, bldg)
      rows[[length(rows) + 1]] <- tibble(
        sensor_id = p$sensor_id, region = p$region, radius_m = R,
        imperv = imp, tree = can[["1"]], water = wat,
        bldg = bm$bldg_frac, height = bm$zH, height_max = bm$zH_max, n_bldg = bm$n_bldg,
        aspect = bm$hw_ratio, svf_point = svf, coast_km = coast_d / 1000)
    }
    message(sprintf("  %-8s svf=%.2f coast=%6.0f m", p$sensor_id, svf, coast_d))
  }
  out <- bind_rows(rows)
  write_csv(out, file.path(DER, "site_predictors_multiscale.csv"), na = "")
  message("wrote ", nrow(out), " rows")
  print(out %>% group_by(region, radius_m) %>%
          summarise(across(c(imperv, tree, bldg, height, svf_point, aspect), ~ round(mean(.x), 2)), .groups = "drop"))
  invisible(out)
}

if (sys.nframe() == 0) main()
