#!/usr/bin/env Rscript
# 07_sensitivity_checks.R -------------------------------------------------------
# Two site-level sensitivity checks on the models of 02_models.R, at the adopted
# 100 m radius, with the same z-scoring (pooled SDs of the 400 sensor-nights),
# the same random-intercept structure (1 | sensor_id), likelihood-ratio tests on
# ML fits and REML estimates with Satterthwaite p-values.
#
#   A  HNL08.  The warmest Honolulu site, variable from night to night and not
#      explained by its descriptors; a generator next to the station was likely the
#      source of its extra warmth.  The models are refitted without its seven
#      sensor-nights.  Its readings also entered the Honolulu network median; the
#      check reports how far the median moves when the logger is left out of it.
#
#   B  Buildings without a recorded height.  About a quarter of the Oʻahu buildings
#      have no HEIGHT in the FEMA/ORNL USA Structures inventory; in the footprint file
#      used here (data/external/buildings_oahu.gpkg) they were given 5 m during preprocessing.
#      For each site the share of building area within 100 m that carries the fill
#      value is computed, and the models are refitted (i) with 3 m instead of 5 m for
#      these buildings (height and canyon aspect ratio recomputed at 100 m) and
#      (ii) without the sites where such buildings make up more than half of the
#      building area.
#
#   C  Sensor positions inside building footprints.  Five recorded logger positions
#      fall inside a footprint polygon (a post under an eave or a covered walkway, or
#      a position error of a few metres), so the ray-casting treats the building as
#      an obstruction 1 m away in every direction and returns a sky view factor
#      near zero.  The sky view factor of these sites is recomputed with the
#      position moved to 1 m outside the nearest footprint edge, and the models
#      that contain the sky view factor are refitted.
#
# Outputs: results/sensitivity_checks.json,
#          results/tables/sensitivity_*.csv,
#          data/descriptors/site_default_height_share.csv, data/descriptors/site_svf_footprint_check.csv
#
# Packages: as 02_models.R (lme4, lmerTest, performance, car), plus sf for part B.
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "02_models.R"), encoding = "UTF-8")     # fit_lmm(), lrt(), coef_table(), r2_nakagawa(), ...

SC <- list()
DEFAULT_HEIGHT_M <- 5            # the repeated value in the footprint layer
DEFAULT_SHARE_MAX <- 0.5         # part B: sites above this share of default-height building area are left out
# CRS_M comes from helpers.R
QUADSEGS <- 16                   # 64-vertex circles, as in 01_site_predictors.R

load_analysis <- function() {
  d <- read_csv(file.path(TAB, "analysis_dataset.csv"), col_types = cols(night_date = col_character(), .default = col_guess())) %>%
    mutate(region = factor(region, levels = c("Ewa", "Honolulu")))
  stopifnot(nrow(d) == 400)
  res <- read_json(file.path(OUT, "thesis_results.json"))
  list(d = d, sds = unlist(res$z_sds), res = res)
}

nat_row <- function(ct, term) {
  r <- ct[ct$term == term, ]
  if (!nrow(r)) return(list(est = NA, se = NA, p = NA))
  list(est = r$estimate_natural, se = r$se_natural, p = r$p, est_z = r$estimate)
}

# the set of models whose results the thesis reports, fitted to a data set `d`
model_set <- function(d, sds, res) {
  out <- list(n_sites = n_distinct(d$sensor_id), n_obs = nrow(d))
  # pooled full model (district + eight descriptors)
  m <- fit_lmm(d, rhs_full("region"), reml = TRUE); m_ml <- fit_lmm(d, rhs_full("region"), reml = FALSE)
  ct <- coef_table(m, sds); r2 <- r2_nakagawa(m)
  out$pooled_full <- list(coefs = setNames(lapply(paste0(PRED, "_z"), function(t) nat_row(ct, t)), PRED),
                          region = nat_row(ct, "regionHonolulu"), r2m = r2$r2m, r2c = r2$r2c,
                          sd_site = sqrt(r2$var_site), sd_resid = sqrt(r2$var_resid))
  # joint district x descriptor interaction test
  f_all <- paste0(rhs_full("region"), " + ", paste0("region:", PRED, "_z", collapse = " + "))
  l <- lrt(fit_lmm(d, f_all, reml = FALSE), m_ml)
  out$joint_interaction <- list(chi2 = l$stat, df = l$df, p = l$p)
  # district models: null (ICC), full, and the best subset reported in the thesis
  for (reg in REGIONS) {
    g <- d %>% filter(region == reg)
    n0 <- r2_nakagawa(fit_lmm(g, "1", reml = TRUE))
    mf <- fit_lmm(g, rhs_full(NULL), reml = TRUE); ctf <- coef_table(mf, sds); r2f <- r2_nakagawa(mf)
    bs <- res$regional[[reg]]$best_subset$predictors
    fb <- paste(paste0(strsplit(bs, " \\+ ")[[1]], "_z"), collapse = " + ")
    mb <- fit_lmm(g, fb, reml = TRUE); ctb <- coef_table(mb, sds)
    # is the same subset still ranked first among all subsets of <= 3 descriptors?
    cands <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(PRED, k, simplify = FALSE)), recursive = FALSE))
    rr <- bind_rows(lapply(cands, function(cc) {
      rhs <- if (length(cc)) paste(paste0(cc, "_z"), collapse = " + ") else "1"
      tibble(predictors = if (length(cc)) paste(cc, collapse = " + ") else "(null)", aic = aic_ml(fit_lmm(g, rhs, reml = FALSE)))
    })) %>% arrange(aic) %>% mutate(delta = aic - min(aic))
    out[[reg]] <- list(
      n_sites = n_distinct(g$sensor_id), icc = n0$var_site / (n0$var_site + n0$var_resid),
      sd_site_null = sqrt(n0$var_site),
      full = list(r2m = r2f$r2m, sd_site = sqrt(r2f$var_site),
                  coefs = setNames(lapply(paste0(PRED, "_z"), function(t) nat_row(ctf, t)), PRED)),
      # marginal R2 of the best subset from its ML fit, as reported in 02_models.R (step 7)
      best_subset = list(predictors = bs, r2m = r2_nakagawa(fit_lmm(g, fb, reml = FALSE))$r2m,
                         coefs = setNames(lapply(strsplit(fb, " \\+ ")[[1]], function(t) nat_row(ctb, t)),
                                          strsplit(bs, " \\+ ")[[1]]),
                         rank_now = which(rr$predictors == bs), top_now = rr$predictors[1],
                         delta_aic_now = rr$delta[rr$predictors == bs]))
  }
  out
}

flat <- function(ms, label) {
  # one row per reported coefficient, for the CSV tables
  rows <- list()
  add <- function(model, term, x) rows[[length(rows) + 1]] <<- tibble(check = label, model = model, term = term,
                                                                      estimate = x$est, se = x$se, p = x$p)
  for (p in PRED) add("pooled full", p, ms$pooled_full$coefs[[p]])
  add("pooled full", "district (Honolulu)", ms$pooled_full$region)
  for (reg in REGIONS) {
    for (p in PRED) add(paste(reg, "full"), p, ms[[reg]]$full$coefs[[p]])
    for (p in names(ms[[reg]]$best_subset$coefs)) add(paste(reg, "best subset"), p, ms[[reg]]$best_subset$coefs[[p]])
  }
  bind_rows(rows)
}

# ---- A: HNL08 -------------------------------------------------------------------------------------------
check_hnl08 <- function(d, sds, res) {
  base <- model_set(d, sds, res)
  sub  <- model_set(d %>% filter(sensor_id != "HNL08"), sds, res)
  # the Honolulu network median with and without HNL08 on the analysed nights
  hh <- read_csv(file.path(PROC, "halfhourly_calm_clear.csv"), show_col_types = FALSE,
                 col_types = cols(night_date = col_character(), .default = col_guess())) %>%
    filter(region == "Honolulu", is_night)
  nights <- unique(d$night_date[d$region == "Honolulu"])
  hh <- hh %>% filter(night_date %in% nights)
  med <- hh %>% group_by(night_date, time_bin) %>%
    summarise(m_all = median(temp_c), m_wo = median(temp_c[sensor_id != "HNL08"]), .groups = "drop")
  shift <- med %>% mutate(diff = m_wo - m_all)
  SC$hnl08 <<- list(with = base, without = sub,
                    median_shift = list(max_abs_bin = max(abs(shift$diff)), mean_bin = mean(shift$diff),
                                        night_mean_range = range((shift %>% group_by(night_date) %>% summarise(x = mean(diff)))$x)),
                    site = d %>% filter(sensor_id == "HNL08") %>% summarise(dT_mean = mean(dT_night), dT_sd = sd(dT_night), n = n()) %>% as.list())
  write_table(bind_rows(flat(base, "all 74 sites"), flat(sub, "without HNL08")), "sensitivity_hnl08")
}

# ---- B: default building heights ------------------------------------------------------------------------
default_height_share <- function() {
  f_out <- file.path(DER, "site_default_height_share.csv")
  gpkg <- file.path(EXT, BUILDINGS_FILE)
  if (!file.exists(gpkg)) {
    stopifnot(file.exists(f_out))          # the footprint layer is not in the repository; use the stored shares
    return(read_csv(f_out, show_col_types = FALSE))
  }
  suppressPackageStartupMessages(library(sf))
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), show_col_types = FALSE)
  sites <- sn %>% distinct(sensor_id, .keep_all = TRUE) %>% select(sensor_id, region, latitude, longitude)
  pts <- st_transform(st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)
  bldg <- st_transform(st_read(gpkg, quiet = TRUE), CRS_M)
  all_default <- abs(bldg$height_m - DEFAULT_HEIGHT_M) < 1e-6
  isl <- tibble(n_footprints = nrow(bldg), share_footprints_default = mean(all_default),
                share_area_default = sum(bldg$area_m2[all_default]) / sum(bldg$area_m2))
  rows <- lapply(seq_len(nrow(pts)), function(i) {
    g <- st_buffer(st_geometry(pts[i, ]), 100, nQuadSegs = QUADSEGS)
    idx <- st_intersects(g, bldg)[[1]]
    if (!length(idx)) return(tibble(sensor_id = pts$sensor_id[i], region = pts$region[i], n_bldg = 0L,
                                    area_bldg_m2 = 0, share_default = NA_real_))
    sub <- bldg[idx, ]
    a <- as.numeric(st_area(suppressWarnings(st_intersection(st_geometry(sub), g))))
    dflt <- abs(sub$height_m - DEFAULT_HEIGHT_M) < 1e-6
    tibble(sensor_id = pts$sensor_id[i], region = pts$region[i], n_bldg = length(idx),
           area_bldg_m2 = sum(a), share_default = sum(a[dflt]) / sum(a))
  })
  out <- bind_rows(rows)
  write_csv(out, f_out, na = "")
  write_csv(isl, file.path(DER, "island_default_height_share.csv"))
  out
}

height_fill_variant <- function(fill_m = 3) {
  # height and canyon aspect ratio at 100 m with the buildings that have no recorded
  # height (5 m in the footprint file) set to `fill_m` instead
  f_out <- file.path(DER, sprintf("site_height_fill%dm_100m.csv", fill_m))
  gpkg <- file.path(EXT, BUILDINGS_FILE)
  if (!file.exists(gpkg)) { stopifnot(file.exists(f_out)); return(read_csv(f_out, show_col_types = FALSE)) }
  suppressPackageStartupMessages(library(sf))
  sp <- new.env(); source(file.path(.here, "01_site_predictors.R"), local = sp, encoding = "UTF-8")
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), show_col_types = FALSE)
  sites <- sn %>% distinct(sensor_id, .keep_all = TRUE) %>% select(sensor_id, latitude, longitude)
  pts <- st_transform(st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)
  bldg <- st_transform(st_read(gpkg, quiet = TRUE), CRS_M) %>%
    mutate(height_m = ifelse(abs(height_m - DEFAULT_HEIGHT_M) < 1e-6, fill_m, height_m),
           height_m = pmin(pmax(coalesce(height_m, 3), 0.5), 60))
  out <- bind_rows(lapply(seq_len(nrow(pts)), function(i) {
    bm <- sp$building_metrics(st_buffer(st_geometry(pts[i, ]), 100, nQuadSegs = QUADSEGS), bldg)
    tibble(sensor_id = pts$sensor_id[i], height = bm$zH, aspect = bm$hw_ratio)
  }))
  write_csv(out, f_out)
  out
}

check_default_heights <- function(d, sds, res) {
  sh <- default_height_share()
  isl_f <- file.path(DER, "island_default_height_share.csv")
  isl <- if (file.exists(isl_f)) as.list(read_csv(isl_f, show_col_types = FALSE)) else NULL
  heavy <- sh$sensor_id[!is.na(sh$share_default) & sh$share_default > DEFAULT_SHARE_MAX]
  # building-area-weighted share of default heights within 100 m, pooled over each district's sites
  by_region <- sh %>% group_by(region) %>%
    summarise(area_weighted_share = sum(area_bldg_m2 * coalesce(share_default, 0)) / sum(area_bldg_m2),
              mean_site_share = mean(share_default, na.rm = TRUE),
              n_sites_over_half = sum(share_default > DEFAULT_SHARE_MAX, na.rm = TRUE), .groups = "drop")
  base <- model_set(d, sds, res)
  sub <- model_set(d %>% filter(!sensor_id %in% heavy), sds, res)
  # (i) 3 m instead of 5 m: height and aspect ratio recomputed and re-standardized (pooled, 400 rows)
  v3 <- height_fill_variant(3)
  d3 <- d %>% select(-height, -aspect) %>% left_join(v3, by = "sensor_id")
  z <- zscore_frame(d3, c("height", "aspect")); d3 <- z$df
  sds3 <- sds; sds3["height"] <- z$sds["height"]; sds3["aspect"] <- z$sds["aspect"]
  fill3 <- model_set(d3, sds3, res)
  SC$default_heights <<- list(default_height_m = DEFAULT_HEIGHT_M, share_threshold = DEFAULT_SHARE_MAX,
                              island = isl, by_region = pred_rows(by_region), sites_left_out = heavy,
                              with = base, without = sub, fill_3m = fill3)
  write_table(bind_rows(flat(base, "all 74 sites"), flat(sub, "without default-height sites"),
                        flat(fill3, "3 m instead of 5 m for buildings without a recorded height")), "sensitivity_default_heights")
}


# ---- C: sensor positions inside building footprints -----------------------------------------------------
svf_footprint_table <- function() {
  f_out <- file.path(DER, "site_svf_footprint_check.csv")
  gpkg <- file.path(EXT, BUILDINGS_FILE)
  if (!file.exists(gpkg)) { stopifnot(file.exists(f_out)); return(read_csv(f_out, show_col_types = FALSE)) }
  suppressPackageStartupMessages(library(sf))
  sp <- new.env()                                                            # sky_view_factor() and its constants,
  source(file.path(.here, "01_site_predictors.R"), local = sp, encoding = "UTF-8")   # kept out of the global environment
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), show_col_types = FALSE)
  sites <- sn %>% distinct(sensor_id, .keep_all = TRUE) %>% select(sensor_id, region, latitude, longitude)
  pts <- st_transform(st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)
  bldg <- st_transform(st_read(gpkg, quiet = TRUE), CRS_M) %>% mutate(height_m = pmin(pmax(coalesce(height_m, 3), 0.5), 60))
  inside <- st_intersects(pts, bldg)
  rows <- lapply(which(lengths(inside) > 0), function(i) {
    p <- st_geometry(pts[i, ]); poly <- st_geometry(bldg[inside[[i]][1], ])
    edge <- st_cast(st_nearest_points(p, st_boundary(poly)), "POINT")[2]          # nearest point on the footprint edge
    xy0 <- st_coordinates(p); xy1 <- st_coordinates(edge); v <- (xy1 - xy0) / sqrt(sum((xy1 - xy0)^2))
    moved <- st_sfc(st_point(as.numeric(xy1 + v)), crs = CRS_M)                    # 1 m beyond the edge, away from the interior
    near <- bldg[st_intersects(st_buffer(moved, sp$HORIZON_R + 5), bldg)[[1]], ]
    tibble(sensor_id = pts$sensor_id[i], region = pts$region[i], footprint_height_m = bldg$height_m[inside[[i]][1]],
           distance_to_edge_m = sqrt(sum((xy1 - xy0)^2)), moved_inside = lengths(st_intersects(moved, bldg)) > 0,
           svf_moved = sp$sky_view_factor(moved, near))
  })
  out <- bind_rows(rows)
  write_csv(out, f_out, na = "")
  out
}

check_svf_footprints <- function(d, sds, res) {
  fp <- svf_footprint_table()
  orig <- d %>% distinct(sensor_id, svf_point)
  fp <- fp %>% left_join(orig, by = "sensor_id") %>% rename(svf_original = svf_point)
  d2 <- d %>% left_join(fp %>% select(sensor_id, svf_moved), by = "sensor_id") %>%
    mutate(svf_point = ifelse(is.na(svf_moved), svf_point, svf_moved)) %>% select(-svf_moved)
  # re-standardise the sky view factor on the corrected values (pooled over the 400 rows, as in 02_models.R)
  z <- zscore_frame(d2, "svf_point"); d2 <- z$df; sds2 <- sds; sds2["svf_point"] <- z$sds["svf_point"]
  single <- function(g, base, s_) {
    rhs <- if (base == "1") "svf_point_z" else paste(base, "svf_point_z", sep = " + ")
    m <- fit_lmm(g, rhs, reml = TRUE); l <- lrt(fit_lmm(g, rhs, reml = FALSE), fit_lmm(g, base, reml = FALSE))
    ct <- coef_table(m, s_); r <- ct[ct$term == "svf_point_z", ]
    list(est = r$estimate_natural, se = r$se_natural, p_lrt = l$p, r2m = r2_nakagawa(m)$r2m)
  }
  singles <- function(dd, s_) {
    list(pooled = single(dd, "region", s_), Honolulu = single(dd %>% filter(region == "Honolulu"), "1", s_),
         Ewa = single(dd %>% filter(region == "Ewa"), "1", s_))
  }
  SC$svf_footprints <<- list(sites = pred_rows(fp),
                             range_other_sites = range(orig$svf_point[!orig$sensor_id %in% fp$sensor_id]),
                             single_original = singles(d, sds), single_corrected = singles(d2, sds2),
                             with = model_set(d, sds, res), corrected = model_set(d2, sds2, res))
  write_table(bind_rows(flat(SC$svf_footprints$with, "original sky view factor"),
                        flat(SC$svf_footprints$corrected, "sky view factor recomputed for the five sites")), "sensitivity_svf_footprints")
}

main <- function() {
  a <- load_analysis()
  message("A: HNL08"); check_hnl08(a$d, a$sds, a$res)
  message("B: default building heights"); check_default_heights(a$d, a$sds, a$res)
  message("C: sensor positions inside footprints"); check_svf_footprints(a$d, a$sds, a$res)
  write_json(SC, file.path(OUT, "sensitivity_checks.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
  message("done")
}

if (sys.nframe() == 0) main()
