#!/usr/bin/env Rscript
# 04_replication_check.R -------------------------------------------------------
# Compares the R results (results/, primary) with the independent Python
# replication of the same pipeline (replication/results: pandas / statsmodels
# MixedLM / rasterio / shapely, kept in replication/python/).  Writes
#   replication/replication_check.csv   (quantity, R, Python, difference)
# and prints a summary.  Differences are expected only where the two
# implementations legitimately differ:
#   * coefficient p-values: lmerTest uses Satterthwaite t tests, statsmodels
#     uses Wald z tests (the estimates and SEs are the same);
#   * OLS AIC: R counts the residual variance as a parameter (+2);
#   * raster fractions: the two libraries differ in which boundary pixels of a
#     circle they count (< 0.2 pp at 50 m, < 0.04 pp at 100 m, < 0.01 pp at 200 m);
#   * the night fixed-effect test: lme4 drops the night dummy that is collinear
#     with district (9 df); statsmodels keeps it (10 df).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")

if (!file.exists(file.path(PY, "thesis_results.json"))) stop("no Python results in ", PY)
R_ <- fromJSON(file.path(OUT, "thesis_results.json"), simplifyVector = FALSE)
P_ <- fromJSON(gsub(":\\s*NaN", ": null", paste(readLines(file.path(PY, "thesis_results.json"), warn = FALSE), collapse = "\n")), simplifyVector = FALSE)  # Python writes NaN

rows <- list()
add <- function(quantity, r, p, unit = "") rows[[length(rows) + 1]] <<- tibble(quantity = quantity, R = as.numeric(r), Python = as.numeric(p), unit = unit)
coef_of <- function(coefs, term) { for (cc in coefs) if (!is.null(cc$term) && cc$term %in% term) return(cc); NULL }
pterm <- function(t) { m <- c("(Intercept)" = "Intercept"); m[REGION_TERM] <- paste0("region[T.", REGIONS[1], "]"); if (t %in% names(m)) m[[t]] else t }   # statsmodels' names

# design and data preparation
for (reg in REGIONS) {
  add(sprintf("%s: sensor-nights", reg), R_$design[[reg]]$n_sensor_nights, P_$design[[reg]]$n_sensor_nights)
  add(sprintf("%s: site-mean dT range (min)", reg), R_$design[[reg]]$dT_site_range[[1]], P_$design[[reg]]$dT_site_range[[1]], "°C")
  add(sprintf("%s: site-mean dT range (max)", reg), R_$design[[reg]]$dT_site_range[[2]], P_$design[[reg]]$dT_site_range[[2]], "°C")
}
add("Adopted radius", R_$adopted_radius, P_$adopted_radius, "m")
for (k in names(R_$scale_rule$max_delta_aic)) add(sprintf("Scale rule: max ΔAIC at %s m", k), R_$scale_rule$max_delta_aic[[k]], P_$scale_rule$max_delta_aic[[k]])
# null models
for (sc in c("pooled", REGIONS)) add(sprintf("Null model ICC (%s)", sc), R_$null[[sc]]$icc, P_$null[[sc]]$icc)
# full pooled model
for (t in c("(Intercept)", REGION_TERM, paste0(PRED, "_z"))) {
  r <- coef_of(R_$full_pooled$coefs, t); p <- coef_of(P_$full_pooled$coefs, pterm(t))
  add(sprintf("Full pooled: %s estimate", t), r$estimate, p$estimate, "°C per SD")
  add(sprintf("Full pooled: %s SE", t), r$se, p$se, "°C per SD")
  add(sprintf("Full pooled: %s p (R: Satterthwaite t; Python: Wald z)", t), r$p, p$p)
}
add("Full pooled: R2m", R_$full_pooled$r2m, P_$full_pooled$r2m); add("Full pooled: R2c", R_$full_pooled$r2c, P_$full_pooled$r2c)
add("Full pooled: AIC (ML)", R_$full_pooled$aic, P_$full_pooled$aic)
add("Full pooled: LRT vs null, chi2", R_$full_pooled$lrt_vs_null[[1]], P_$full_pooled$lrt_vs_null[[1]])
add("Full pooled: LRT morphology | district, chi2", R_$full_pooled$lrt_morphology_given_region[[1]], P_$full_pooled$lrt_morphology_given_region[[1]])
add("Full pooled: site SD (random intercept)", sqrt(R_$full_pooled$var_site), sqrt(P_$full_pooled$var_site), "°C")
# parsimonious
for (t in paste0(PARSIMONIOUS, "_z")) { r <- coef_of(R_$parsimonious_pooled$coefs, t); p <- coef_of(P_$parsimonious_pooled$coefs, t)
  add(sprintf("Parsimonious pooled: %s (natural units)", t), r$estimate_natural, p$estimate_natural, r$natural_unit) }
# district models
for (reg in REGIONS) {
  add(sprintf("%s full model: R2m", reg), R_$regional[[reg]]$full$r2m, P_$regional[[reg]]$full$r2m)
  add(sprintf("%s full model: AIC", reg), R_$regional[[reg]]$full$aic, P_$regional[[reg]]$full$aic)
  add(sprintf("%s best subset: R2m (%s)", reg, R_$regional[[reg]]$best_subset$predictors), R_$regional[[reg]]$best_subset$r2m, P_$regional[[reg]]$best_subset$r2m)
  add(sprintf("%s best subset: Akaike weight of best", reg), R_$regional[[reg]]$best_subset$top5[[1]]$weight, P_$regional[[reg]]$best_subset$top5[[1]]$weight)
  for (cc in R_$regional[[reg]]$best_subset$coefs) if (cc$term != "(Intercept)") {
    p <- coef_of(P_$regional[[reg]]$best_subset$coefs, cc$term)
    add(sprintf("%s best subset: %s (natural units)", reg, cc$term), cc$estimate_natural, p$estimate_natural, cc$natural_unit) }
}
# interactions, night effect
add("Interactions: joint LRT chi2", R_$interactions$joint$lrt_chi2, P_$interactions$joint$lrt_chi2)
for (s in R_$interactions$single) { p <- Filter(function(x) x$predictor == s$predictor, P_$interactions$single)[[1]]
  add(sprintf("Interaction %s: LRT p", s$predictor), s$p, p$p) }
add("Night effect: LRT chi2", R_$night_effect$lrt_chi2, P_$night_effect$lrt_chi2)
add("Night effect: df (lme4 drops the collinear night dummy)", R_$night_effect$df, P_$night_effect$df)
for (reg in REGIONS) add(sprintf("%s: mean between-night correlation of site dT", reg), R_$night_effect$between_night_site_correlation[[reg]]$mean_r, P_$night_effect$between_night_site_correlation[[reg]]$mean_r)
# thresholds
tr <- bind_rows(lapply(R_$threshold_sensitivity, as_tibble)); tp <- bind_rows(lapply(P_$threshold_sensitivity, as_tibble))
for (i in seq_len(nrow(tr))) {
  j <- which(tp$wind_max == tr$wind_max[i] & tp$cloud_max == tr$cloud_max[i] & tp$min_frac == tr$min_frac[i])
  if (length(j) == 1) add(sprintf("Thresholds (wind %g, cloud %g, frac %g): R2m full", tr$wind_max[i], tr$cloud_max[i], tr$min_frac[i]), tr$r2m_full[i], tp$r2m_full[j])
}
# BLUPs, LOSO, calm core, shared nights
add("BLUP SD", R_$blups$sd, P_$blups$sd, "°C")
for (m in R_$loso_cv) { p <- Filter(function(x) x$model == m$model, P_$loso_cv)[[1]]
  add(sprintf("LOSO %s: r", m$model), m$r, p$r); add(sprintf("LOSO %s: RMSE", m$model), m$rmse, p$rmse, "°C") }
for (reg in REGIONS) add(sprintf("Calm core %s: R2m", reg), R_$calm_core[[reg]]$r2m, P_$calm_core[[reg]]$r2m)
for (i in seq_along(R_$shared_nights$median_gap)) add(sprintf("Shared night %s: Honolulu − ʻEwa median", R_$shared_nights$nights[[i]]), R_$shared_nights$median_gap[[i]], P_$shared_nights$median_gap[[i]], "°C")
# predictors (site level)
pr_r <- read_csv(file.path(DER, "site_predictors_multiscale.csv"), show_col_types = FALSE)
pr_p <- read_csv(file.path(PY, "data", "site_predictors_multiscale.csv"), show_col_types = FALSE)
m <- inner_join(pr_r, pr_p, by = c("sensor_id", "radius_m"), suffix = c("_r", "_p"))
for (v in c("imperv", "tree", "water", "bldg", "height", "aspect", "svf_point", "coast_km")) for (R in RADII) {
  g <- m %>% filter(radius_m == R)
  rows[[length(rows) + 1]] <- tibble(quantity = sprintf("Predictor %s at %d m: max |R − Python| over 74 sites", v, R),
                                     R = max(abs(g[[paste0(v, "_r")]] - g[[paste0(v, "_p")]])), Python = NA_real_, unit = "")
}
# sensor-nights
sn_r <- read_csv(file.path(PROC, "sensor_nights.csv"), show_col_types = FALSE); sn_p <- read_csv(file.path(PY, "data", "sensor_nights.csv"), show_col_types = FALSE)
mm <- inner_join(sn_r, sn_p, by = c("sensor_id", "night_date"), suffix = c("_r", "_p"))
rows[[length(rows) + 1]] <- tibble(quantity = "Sensor-night dT: max |R − Python| over 400 rows", R = max(abs(mm$dT_night_r - mm$dT_night_p)), Python = NA_real_, unit = "°C")

t <- bind_rows(rows) %>% mutate(difference = R - Python)
write_csv(t, file.path(ROOT, "replication", "replication_check.csv"), na = "")
cat(sprintf("%d quantities compared; %d differ by more than 0.01 (listed below)\n", nrow(t), sum(abs(t$difference) > 0.01, na.rm = TRUE)))
print(as.data.frame(t %>% filter(abs(difference) > 0.01) %>% mutate(across(where(is.numeric), ~ signif(.x, 4)))), row.names = FALSE)
