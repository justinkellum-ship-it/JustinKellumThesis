#!/usr/bin/env Rscript
# 06_elevation_check.R -----------------------------------------------------------
# Terrain sensitivity check: does site elevation account for the effects that the
# thesis attributes to urban form, and in particular for the distance-to-coast
# effect in Honolulu, where the shore-to-valley gradient runs up-slope?
#
# Elevation is the bare-earth elevation of each logger position from the USGS
# 3DEP 1 m lidar DEM, read through the USGS Elevation Point Query Service
# (data/descriptors/site_elevation.csv; the file records the raster id and the date
# of retrieval).  It is added to the models of 02_models.R at the adopted
# 100 m radius, with the same z-scoring, the same random-intercept structure
# (1 | sensor_id), the same likelihood-ratio tests on ML fits and the same
# REML estimates with Satterthwaite p-values.
#
# Steps
#   1  descriptives: elevation by district, correlation with the eight descriptors
#      and with site-mean dT, variance inflation with elevation added
#   2  elevation alone (pooled with district, Honolulu, ʻEwa)
#   3  coast vs elevation: coast alone, elevation alone, both together
#   4  the adopted models with elevation added (parsimonious pooled, full pooled,
#      district full and best-subset models): LRT, AIC, coefficient shifts
#   5  district x elevation interaction
#   6  best subsets re-run with elevation in the candidate pool (<= 3 terms)
#   7  leave-one-site-out cross-validation with and without elevation
#
# Outputs: results/elevation_check.json, results/tables/elevation_*.csv,
#          results/figures/T20_elevation.png|pdf
#
# Packages: as 02_models.R (lme4, lmerTest, performance, car) plus ggplot2 and
#           patchwork for the figure.
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "02_models.R"), encoding = "UTF-8")     # helpers: fit_lmm(), lrt(), coef_table(), r2_nakagawa(), ...
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

E <- list()
UNIT$elev_m <- list(10, "per +10 m")                           # natural-unit conversion for the new term
LABEL["elev_m"] <- "Elevation (m)"
PRED9 <- c(PRED, "elev_m")

# ---- data ---------------------------------------------------------------------------------
load_elev <- function() {
  # the analysis data set of 02_models.R (400 sensor-nights, 100 m descriptors already z-scored
  # on the pooled SDs) joined to the site elevations, which are z-scored the same way
  d <- read_csv(file.path(TAB, "analysis_dataset.csv"), col_types = cols(night_date = col_character(), .default = col_guess())) %>%
    mutate(region = factor(region, levels = c("Ewa", "Honolulu")))
  el <- read_csv(file.path(DER, "site_elevation.csv"), show_col_types = FALSE) %>% select(sensor_id, elev_m)
  stopifnot(all(unique(d$sensor_id) %in% el$sensor_id))
  d <- d %>% inner_join(el, by = "sensor_id")
  stopifnot(nrow(d) == 400)
  res <- read_json(file.path(OUT, "thesis_results.json"))
  sds <- unlist(res$z_sds); means <- unlist(res$z_means)
  z <- zscore_frame(d, "elev_m")                                # pooled mean and SD over the 400 rows, as for the others
  d <- z$df; sds["elev_m"] <- z$sds["elev_m"]; means["elev_m"] <- z$means["elev_m"]
  list(d = d, sds = sds, means = means, res = res)
}
scopes <- function(d) list(pooled = d, Honolulu = d %>% filter(region == "Honolulu"), Ewa = d %>% filter(region == "Ewa"))
base_rhs <- function(scope) if (scope == "pooled") "region" else "1"
with_term <- function(rhs, term) if (rhs %in% c("1", "")) term else paste(rhs, term, sep = " + ")
nat <- function(ct, term) { r <- ct[ct$term == term, ]; if (!nrow(r)) return(list(est = NA, se = NA, p = NA))
  list(est = r$estimate_natural, se = r$se_natural, p = r$p, est_z = r$estimate, se_z = r$se) }

# ---- step 1: descriptives ------------------------------------------------------------------
step1_descr <- function(d) {
  sites <- d %>% group_by(region, sensor_id) %>%
    summarise(dT_site = mean(dT_night), across(all_of(PRED9), first), .groups = "drop")
  desc <- sites %>% group_by(region) %>%
    summarise(n = n(), min = min(elev_m), q25 = quantile(elev_m, .25), median = median(elev_m), mean = mean(elev_m),
              q75 = quantile(elev_m, .75), max = max(elev_m), sd = sd(elev_m), .groups = "drop") %>% mutate(region = as.character(region))
  corr <- bind_rows(lapply(c("pooled", REGIONS), function(scope) {
    s <- if (scope == "pooled") sites else sites %>% filter(region == scope)
    tibble(scope = scope, variable = c(PRED, "dT_site"),
           r = sapply(c(PRED, "dT_site"), function(v) cor(s$elev_m, s[[v]])))
  }))
  vif <- bind_rows(lapply(c("pooled", REGIONS), function(scope) {
    s <- if (scope == "pooled") sites else sites %>% filter(region == scope)
    sz <- as.data.frame(scale(s[PRED9])); sz$y <- rnorm(nrow(sz))
    v <- car::vif(lm(y ~ ., data = sz))
    tibble(scope = scope, variable = PRED9, vif_with_elevation = as.numeric(v[PRED9]))
  }))
  write_table(desc, "elevation_descriptives"); write_table(corr, "elevation_correlations"); write_table(vif, "elevation_vif")
  E$descriptives <<- pred_rows(desc)
  E$correlations <<- pred_rows(corr)
  E$vif <<- pred_rows(vif)
  E$sites <<- pred_rows(sites %>% mutate(region = as.character(region)) %>% select(region, sensor_id, dT_site, elev_m, coast_km))
  sites
}

# ---- steps 2-4: models with elevation --------------------------------------------------------
fit_pair <- function(g, rhs, sds) list(ml = fit_lmm(g, rhs, reml = FALSE), reml = fit_lmm(g, rhs, reml = TRUE),
                                       ct = coef_table(fit_lmm(g, rhs, reml = TRUE), sds))
model_row <- function(scope, label, rhs, g, sds, null_ml = NULL) {
  f <- fit_pair(g, rhs, sds); r2 <- r2_nakagawa(f$reml)
  co <- nat(f$ct, "coast_km_z"); el <- nat(f$ct, "elev_m_z"); im <- nat(f$ct, "imperv_z"); he <- nat(f$ct, "height_z")
  tibble(scope = scope, model = label, rhs = rhs, k = n_fixed(f$ml), aic = aic_ml(f$ml), r2m = r2$r2m, r2c = r2$r2c,
         coast_per_km = co$est, coast_se = co$se, coast_p = co$p,
         elev_per_10m = el$est, elev_se = el$se, elev_p = el$p,
         imperv_per_10pp = im$est, imperv_p = im$p, height_per_m = he$est, height_p = he$p,
         lrt_chi2 = if (is.null(null_ml)) NA_real_ else lrt(f$ml, null_ml)$stat,
         lrt_df = if (is.null(null_ml)) NA_real_ else lrt(f$ml, null_ml)$df,
         lrt_p = if (is.null(null_ml)) NA_real_ else lrt(f$ml, null_ml)$p,
         converged = !length(f$reml@optinfo$conv$lme4$messages))
}

step2to4_models <- function(d, sds, res) {
  rows <- list(); coefs <- list()
  for (scope in c("pooled", REGIONS)) {
    g <- scopes(d)[[scope]]; b <- base_rhs(scope)
    null_ml <- fit_lmm(g, b, reml = FALSE)
    # step 2: elevation alone; step 3: coast alone, and the two together
    m_el <- with_term(b, "elev_m_z"); m_co <- with_term(b, "coast_km_z"); m_both <- with_term(b, "coast_km_z + elev_m_z")
    rows[[length(rows) + 1]] <- model_row(scope, "null", b, g, sds)
    rows[[length(rows) + 1]] <- model_row(scope, "elevation only", m_el, g, sds, null_ml)
    rows[[length(rows) + 1]] <- model_row(scope, "coast only", m_co, g, sds, null_ml)
    f_co <- fit_lmm(g, m_co, reml = FALSE); f_el <- fit_lmm(g, m_el, reml = FALSE); f_both <- fit_lmm(g, m_both, reml = FALSE)
    r <- model_row(scope, "coast + elevation", m_both, g, sds, null_ml)
    r$lrt_elev_given_coast_chi2 <- lrt(f_both, f_co)$stat; r$lrt_elev_given_coast_p <- lrt(f_both, f_co)$p
    r$lrt_coast_given_elev_chi2 <- lrt(f_both, f_el)$stat; r$lrt_coast_given_elev_p <- lrt(f_both, f_el)$p
    rows[[length(rows) + 1]] <- r
    # step 4: the adopted models, without and with elevation
    adopted <- if (scope == "pooled") {
      list(c("parsimonious (district + impervious + height)", "region + imperv_z + height_z"),
           c("full (district + eight descriptors)", rhs_full("region")))
    } else {
      bs <- res$regional[[scope]]$best_subset$predictors
      list(c(sprintf("best subset (%s)", gsub("coast_km", "coast", bs)), paste(paste0(strsplit(bs, " \\+ ")[[1]], "_z"), collapse = " + ")),
           c("full (eight descriptors)", rhs_full(NULL)))
    }
    for (a in adopted) {
      m0 <- fit_lmm(g, a[2], reml = FALSE)
      rows[[length(rows) + 1]] <- model_row(scope, a[1], a[2], g, sds, null_ml)
      r <- model_row(scope, paste(a[1], "+ elevation"), paste(a[2], "+ elev_m_z"), g, sds, null_ml)
      l <- lrt(fit_lmm(g, paste(a[2], "+ elev_m_z"), reml = FALSE), m0)
      r$lrt_elev_given_model_chi2 <- l$stat; r$lrt_elev_given_model_p <- l$p; r$delta_aic_vs_without <- r$aic - aic_ml(m0)
      rows[[length(rows) + 1]] <- r
      # full coefficient tables of the with-elevation fits, for the appendix
      ct <- coef_table(fit_lmm(g, paste(a[2], "+ elev_m_z"), reml = TRUE), sds) %>% mutate(scope = scope, model = paste(a[1], "+ elevation"))
      coefs[[length(coefs) + 1]] <- ct
    }
  }
  t <- bind_rows(rows) %>% mutate(delta_aic_vs_null = aic - ave(aic, scope, FUN = function(a) a[1]))
  write_table(t, "elevation_models"); write_table(bind_rows(coefs), "elevation_model_coefficients")
  E$models <<- pred_rows(t)
  t
}

# ---- step 5: district x elevation interaction ------------------------------------------------------
step5_interaction <- function(d, sds) {
  base <- rhs_full("region")
  m0 <- fit_lmm(d, paste(base, "+ elev_m_z"), reml = FALSE)
  m1 <- fit_lmm(d, paste(base, "+ elev_m_z + region:elev_m_z"), reml = FALSE); l <- lrt(m1, m0)
  mr <- fit_lmm(d, paste(base, "+ elev_m_z + region:elev_m_z"), reml = TRUE); fe <- fixef(mr)
  f <- UNIT$elev_m[[1]] / sds["elev_m"]
  E$interaction <<- list(lrt_chi2 = l$stat, df = l$df, p = l$p,
                         ewa_slope_per_10m = unname(fe["elev_m_z"] * f),
                         hnl_slope_per_10m = unname((fe["elev_m_z"] + fe["regionHonolulu:elev_m_z"]) * f))
}

# ---- step 6: best subsets with elevation in the pool --------------------------------------------------
step6_subsets <- function(d, res) {
  out <- list(); tabs <- list()
  for (reg in REGIONS) {
    g <- d %>% filter(region == reg)
    cands <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(PRED9, k, simplify = FALSE)), recursive = FALSE))
    rr <- bind_rows(lapply(cands, function(cc) {
      rhs <- if (length(cc)) paste(paste0(cc, "_z"), collapse = " + ") else "1"
      mm <- fit_lmm(g, rhs, reml = FALSE)
      tibble(region = reg, predictors = if (length(cc)) paste(cc, collapse = " + ") else "(null)", k = length(cc),
             aic = aic_ml(mm), r2m = r2_nakagawa(mm)$r2m)
    })) %>% mutate(delta_aic = aic - min(aic), weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))) %>% arrange(aic) %>%
      mutate(rank = row_number())
    tabs[[reg]] <- rr
    orig <- res$regional[[reg]]$best_subset$predictors
    out[[reg]] <- list(n_models = nrow(rr), best = rr$predictors[1], best_aic = rr$aic[1], best_r2m = rr$r2m[1],
                       original_best = orig, original_rank = rr$rank[rr$predictors == orig], original_delta_aic = rr$delta_aic[rr$predictors == orig],
                       elevation_importance = sum(rr$weight[grepl("elev_m", rr$predictors, fixed = TRUE)]),
                       coast_importance = sum(rr$weight[grepl("coast_km", rr$predictors, fixed = TRUE)]),
                       best_with_elevation = (rr %>% filter(grepl("elev_m", predictors, fixed = TRUE)))$predictors[1],
                       best_with_elevation_rank = (rr %>% filter(grepl("elev_m", predictors, fixed = TRUE)))$rank[1],
                       best_with_elevation_delta_aic = (rr %>% filter(grepl("elev_m", predictors, fixed = TRUE)))$delta_aic[1],
                       n_within_2 = sum(rr$delta_aic <= 2), top5 = pred_rows(head(rr, 5)))
  }
  write_table(bind_rows(tabs), "elevation_best_subsets")
  E$subsets <<- out
}

# ---- step 7: leave-one-site-out with and without elevation --------------------------------------------
step7_loso <- function(d, res) {
  site <- d %>% group_by(sensor_id) %>% summarise(dT_site = mean(dT_night), region = first(region)) %>% arrange(sensor_id)
  predict_loso <- function(dd, rhs, ids) vapply(ids, function(sid) {
    m <- fit_lmm(dd %>% filter(sensor_id != sid), rhs, reml = TRUE)
    unname(predict(m, newdata = dd %>% filter(sensor_id == sid) %>% slice(1), re.form = NA)) }, numeric(1))
  cv_row <- function(scope, model, pred, obs) { err <- pred - obs
    tibble(scope = scope, model = model, rmse = sqrt(mean(err^2)), mae = mean(abs(err)), r = cor(pred, obs),
           r2_cv = 1 - sum(err^2) / sum((obs - mean(obs))^2)) }
  rows <- list()
  for (s in list(c("pooled", "parsimonious", "region + imperv_z + height_z"), c("pooled", "full", rhs_full("region")))) {
    for (add in c("", " + elev_m_z")) {
      pr <- predict_loso(d, paste0(s[3], add), site$sensor_id)
      rows[[length(rows) + 1]] <- cv_row(s[1], paste0(s[2], if (nzchar(add)) " + elevation" else ""), pr, site$dT_site)
    }
  }
  for (reg in REGIONS) {
    g <- d %>% filter(region == reg); sr <- site %>% filter(region == reg)
    bs <- res$regional[[reg]]$best_subset$predictors
    f_bs <- paste(paste0(strsplit(bs, " \\+ ")[[1]], "_z"), collapse = " + ")
    for (s in list(c(sprintf("best subset (%s)", gsub("coast_km", "coast", bs)), f_bs), c("full", rhs_full(NULL)))) {
      for (add in c("", " + elev_m_z")) {
        pr <- predict_loso(g, paste0(s[2], add), sr$sensor_id)
        rows[[length(rows) + 1]] <- cv_row(reg, paste0(s[1], if (nzchar(add)) " + elevation" else ""), pr, sr$dT_site)
      }
    }
    # elevation in place of coast (Honolulu) / added to the null (both): the terrain-only model
    pr <- predict_loso(g, "elev_m_z", sr$sensor_id)
    rows[[length(rows) + 1]] <- cv_row(reg, "elevation only", pr, sr$dT_site)
  }
  t <- bind_rows(rows); write_table(t, "elevation_loso")
  E$loso <<- pred_rows(t)
}

# ---- figure ---------------------------------------------------------------------------------------------
fig_elevation <- function(sites, t) {
  # Purpose: show the raw relationship of site-mean dT with elevation, with the slope of the
  # elevation-only mixed model (solid where p < 0.05, dotted otherwise), and how elevation
  # and distance to the coast are related within each district.
  source(file.path(.here, "maplib.R"), encoding = "UTF-8")      # theme_thesis(), save_fig(), REG_COL, REG_LABEL
  reg_f <- function(x) factor(x, levels = REGIONS, labels = REG_LABEL[REGIONS])
  sites <- sites %>% mutate(regionf = reg_f(as.character(region)))
  slope <- t %>% filter(scope %in% REGIONS, model == "elevation only") %>% transmute(region = scope, b = elev_per_10m / 10, p = elev_p)
  lines <- sites %>% group_by(region = as.character(region), regionf) %>%
    summarise(xmin = min(elev_m), xmax = max(elev_m), xm = mean(elev_m), ym = mean(dT_site), .groups = "drop") %>%
    inner_join(slope, by = "region") %>%
    mutate(y0 = ym + b * (xmin - xm), y1 = ym + b * (xmax - xm), sig = ifelse(p < 0.05, "LMM slope, p < 0.05", "p ≥ 0.05"))
  cols <- setNames(REG_COL[REGIONS], REG_LABEL[REGIONS])
  pa <- ggplot(sites, aes(elev_m, dT_site, colour = regionf)) + theme_thesis() +
    geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.25) +
    geom_point(size = 1.2, alpha = 0.75) +
    geom_segment(data = lines, aes(x = xmin, xend = xmax, y = y0, yend = y1, colour = regionf, linetype = sig), linewidth = 0.7) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_linetype_manual(values = c("LMM slope, p < 0.05" = "solid", "p ≥ 0.05" = "dotted"), name = NULL, drop = FALSE) +
    labs(x = "Site elevation (m above sea level)", y = "Site-mean ΔT (°C)", title = "(a) Site-mean ΔT against elevation") +
    theme(legend.position = "bottom", legend.box = "horizontal")
  pb <- ggplot(sites, aes(coast_km, elev_m, colour = regionf)) + theme_thesis() +
    geom_point(size = 1.2, alpha = 0.75) +
    scale_colour_manual(values = cols, name = NULL, guide = "none") +
    labs(x = "Distance to coast (km)", y = "Site elevation (m)", title = "(b) Elevation against distance to the coast")
  p <- (pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  save_fig(p, "T20_elevation", width = 7.2, height = 3.3)
}

# ---- main ---------------------------------------------------------------------------------------------
main <- function() {
  x <- load_elev(); d <- x$d; sds <- x$sds; res <- x$res
  E$source <<- list(dem = "USGS 3DEP 1 m lidar DEM, bare earth, via the USGS Elevation Point Query Service (epqs.nationalmap.gov, v1)",
                    retrieved = "2026-09-19", radius_m = res$adopted_radius, n_sites = n_distinct(d$sensor_id), n_obs = nrow(d),
                    z_sd_elev = unname(sds["elev_m"]), z_mean_elev = unname(x$means["elev_m"]))
  message("Step 1: descriptives"); sites <- step1_descr(d)
  message("Steps 2-4: models"); t <- step2to4_models(d, sds, res)
  message("Step 5: interaction"); step5_interaction(d, sds)
  message("Step 6: best subsets"); step6_subsets(d, res)
  message("Step 7: LOSO"); step7_loso(d, res)
  message("Figure"); fig_elevation(sites, t)
  write_json(E, file.path(OUT, "elevation_check.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
  message("done")
}

if (sys.nframe() == 0) main()
