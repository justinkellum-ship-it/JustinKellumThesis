#!/usr/bin/env Rscript
# 05_thesis_tables.R -----------------------------------------------------------
# Markdown tables for the thesis document, generated from
# results/thesis_results.json and results/tables/*.csv so that
# every number printed in the thesis comes from the analysis files.
# Output: results/thesis_tables/<name>.md (pandoc pipe tables with captions, as they appear in the thesis).
#
# Packages: dplyr/tidyr (shaping), jsonlite (results file), readr (CSV).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")

DST <- file.path(OUT, "thesis_tables"); dir.create(DST, recursive = TRUE, showWarnings = FALSE)
R <- fromJSON(file.path(OUT, "thesis_results.json"), simplifyVector = FALSE)
NAME <- c(imperv = "Impervious surface fraction (%)", tree = "Tree canopy fraction (%)", bldg = "Building footprint fraction (%)",
          water = "Water surface fraction (%)", height = "Mean building height (m)", svf_point = "Sky view factor (0–1)",
          aspect = "Canyon aspect ratio (H/W)", coast_km = "Distance to coast (km)")
REG <- c(Honolulu = "Honolulu", Ewa = "ʻEwa")
RAD <- R$adopted_radius

p_fmt <- function(p) ifelse(p < 0.001, "< 0.001", ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p)))
f2 <- function(x, d = 2) sprintf(paste0("%.", d, "f"), x)
s2 <- function(x, d = 2) sprintf(paste0("%+.", d, "f"), x)
coefs_df <- function(cl) bind_rows(lapply(cl, function(cc) as_tibble(lapply(cc, function(v) if (is.null(v)) NA else v))))
bold <- function(s, yes) ifelse(yes, paste0("**", s, "**"), s)

md_table <- function(df, caption, name, align = NULL) {
  # pandoc pipe table; the separator dashes carry relative column widths
  df <- as.data.frame(df, stringsAsFactors = FALSE); for (j in seq_along(df)) df[[j]] <- ifelse(is.na(df[[j]]), "", as.character(df[[j]]))
  cols <- names(df)
  # a "|" inside a header or cell would split the pipe-table column and shift the columns after it
  bad <- c(cols[grepl("|", cols, fixed = TRUE)], unlist(lapply(df, function(v) v[grepl("|", v, fixed = TRUE)])))
  if (length(bad)) stop("table ", name, ": '|' inside a cell: ", paste(bad, collapse = "; "))
  if (is.null(align)) align <- c("---", rep("---:", length(cols) - 1))
  # relative column widths from the longest entry, bounded so that short columns (e.g. the study area)
  # are not squeezed by one long text column
  widths <- pmax(12, pmin(40, vapply(seq_along(cols), function(j) max(nchar(c(cols[j], df[[j]]))), numeric(1))))
  sep <- vapply(seq_along(cols), function(j) { a <- align[j]; paste0(if (startsWith(a, ":")) ":" else "", strrep("-", widths[j]), if (endsWith(a, ":")) ":" else "") }, character(1))
  lines <- c(paste0("Table: ", caption, " {#tbl:", name, "}"), "",
             paste0("| ", paste(cols, collapse = " | "), " |"), paste0("|", paste(sep, collapse = "|"), "|"),
             vapply(seq_len(nrow(df)), function(i) paste0("| ", paste(df[i, ], collapse = " | "), " |"), character(1)))
  writeLines(lines, file.path(DST, paste0(name, ".md")), useBytes = TRUE)
  message("wrote ", name)
}

# T1 night inventory ---------------------------------------------------------------
inv <- read_csv(file.path(PROC, "night_inventory.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
sel <- inv %>% filter(selected | frac_calm_clear >= 0.75) %>%
  transmute(`Study area` = REG[region], `Night (18:00–06:00 HST)` = trimws(format(as.Date(night_date), "%e %b %Y")),
            `Sensors reporting` = n_sensors, `Night-mean wind (km h⁻¹)` = f2(wind_night, 1), `Night-mean cloud (%)` = round(cloud_night),
            `Share of calm & clear bins` = paste0(round(frac_calm_clear * 100), "%"),
            Decision = recode(decision, "selected" = "analysed", "partial night (deployment/retrieval)" = "excluded: partial night (loggers still being deployed)",
                              "fewer than 20 sensors reporting" = "excluded: one logger still recording"))
md_table(sel, "Nights on which at least 75% of the available night bins were calm and clear (ERA5 10m wind below 10km h⁻¹ and cloud cover below 25%), with the completeness decision (all 24 bins present; at least 20 loggers reporting). Diurnal dates start at 06:00 HST.",
         "nights", align = c("---", "---", "---:", "---:", "---:", "---:", "---"))

# T2 deployment summary --------------------------------------------------------------
hh <- read_csv(file.path(RAW, "logger_readings.csv"), col_types = cols(datetime = col_character(), .default = col_guess())) %>% mutate(datetime = read_stamps(datetime))
des <- lapply(REGIONS, function(reg) { d <- R$design[[reg]]; h <- hh %>% filter(region == reg)
  c(`Loggers deployed / recovered` = sprintf("50 / %d", d$n_sites),
    Record = sprintf("%s – %s", trimws(format(min(h$datetime), "%e %b %Y")), trimws(format(max(h$datetime), "%e %b %Y"))),
    `Half-hourly readings` = format(nrow(h), big.mark = ","), `Calm, clear nights analysed` = as.character(d$n_nights),
    `Sensor-nights` = as.character(d$n_sensor_nights), `Night-mean air temperature (°C)` = f2(d$temp_night_mean, 1),
    `Night-mean wind / cloud on analysed nights` = sprintf("%.1fkm h⁻¹ / %.0f%%", d$wind_mean, d$cloud_mean)) })
t <- tibble(` ` = names(des[[1]]), Honolulu = unname(des[[1]]), `ʻEwa` = unname(des[[2]]))
md_table(t, "Summary of the two deployments and of the analysed record.", "design", align = c("---", "---", "---"))

# T4 predictor descriptives --------------------------------------------------------
pb <- read_csv(file.path(TAB, "predictors_by_region.csv"), show_col_types = FALSE) %>% filter(radius_m == RAD)
rows <- bind_rows(lapply(PRED, function(p) {
  h <- pb %>% filter(region == "Honolulu", predictor == p); e <- pb %>% filter(region == "Ewa", predictor == p)
  d <- if (p %in% c("svf_point", "aspect", "coast_km")) 2 else 1
  tibble(Predictor = NAME[[p]], `Honolulu mean ± SD` = sprintf("%s ± %s", f2(h$mean, d), f2(h$sd, d)), `Honolulu range` = sprintf("%s – %s", f2(h$min, d), f2(h$max, d)),
         `ʻEwa mean ± SD` = sprintf("%s ± %s", f2(e$mean, d), f2(e$sd, d)), `ʻEwa range` = sprintf("%s – %s", f2(e$min, d), f2(e$max, d)),
         `Honolulu vs ʻEwa, p` = p_fmt(R$region_difference_100m[[p]]$p)) }))
md_table(rows, sprintf("Site predictors at the adopted %dm source-area radius, by study area (sky view factor and coastal distance are point properties of the site). The last column tests whether a descriptor differs between the two districts' sites: it is the p-value of a Mann–Whitney (Wilcoxon rank-sum) test, and small values mean the districts differ.", RAD),
         "predictors", align = c("---", rep("---:", 5)))

# T5 scale selection ---------------------------------------------------------------
sc <- coefs_df(R$scale_selection)
st <- read_csv(file.path(TAB, "stepwise_ols_site_means.csv"), show_col_types = FALSE) %>% group_by(region, radius_m) %>% slice_max(step, n = 1) %>% ungroup()
rows <- bind_rows(lapply(list(c("pooled", "Pooled + study area"), c("Honolulu", "Honolulu"), c("Ewa", "ʻEwa")), function(s) bind_rows(lapply(RADII, function(Rr) {
  g <- sc %>% filter(scope == s[1], radius_m == Rr)
  tibble(Scope = if (Rr == 50) s[2] else "", `Radius (m)` = Rr, `AIC (ML)` = f2(g$aic, 1), `ΔAIC` = f2(g$delta_aic, 1), `R²m` = f2(g$r2m), `R²c` = f2(g$r2c),
         `Stepwise OLS R² (site means)` = if (s[1] != "pooled") f2(st$r2[st$region == s[1] & st$radius_m == Rr]) else "") }))))
md_table(rows, "Source-area scale selection. AIC and ΔAIC (difference from the best radius within each scope) of the full mixed model (all eight predictors, random intercept per site) fitted by maximum likelihood; R²m and R²c are its marginal and conditional R². The last column gives the R² of a forward-stepwise ordinary least squares regression on site means, fitted as a cross-check.",
         "scale", align = c("---", rep("---:", 6)))

# T6 collinearity -------------------------------------------------------------------
rows <- bind_rows(lapply(PRED, function(p) tibble(Predictor = NAME[[p]],
  `VIF Honolulu` = f2(R$collinearity$Honolulu$vif[[p]], 1), `Largest correlation, Honolulu` = f2(R$collinearity$Honolulu$max_abs_corr[[p]]),
  `VIF ʻEwa` = f2(R$collinearity$Ewa$vif[[p]], 1), `Largest correlation, ʻEwa` = f2(R$collinearity$Ewa$max_abs_corr[[p]]), `VIF pooled` = f2(R$collinearity$pooled$vif[[p]], 1))))
md_table(rows, "Variance inflation factors (VIF) of the eight predictors, for each district and pooled, and the largest correlation of each predictor with any of the other seven in each district (absolute value of Pearson's r), at the site level.",
         "vif", align = c("---", rep("---:", 5)))

# T7 null models --------------------------------------------------------------------
rows <- bind_rows(lapply(list(c("pooled", "Pooled (74 sites, 400 sensor-nights)"), c("Honolulu", "Honolulu (36 sites, 248)"), c("Ewa", "ʻEwa (38 sites, 152)")), function(s) {
  n <- R$null[[s[1]]]
  tibble(Model = s[2], `Between-site variance σ²ᵤ (°C²)` = f2(n$var_site, 3), `Between-site SD (°C)` = f2(n$site_sd),
         `Residual variance σ²ε (°C²)` = f2(n$var_resid, 3), `Residual SD (°C)` = f2(n$resid_sd), ICC = f2(n$icc)) }))
md_table(rows, "Variance components of the intercept-only (null) mixed models (REML) and the intraclass correlation coefficient ICC = σ²ᵤ / (σ²ᵤ + σ²ε).",
         "null", align = c("---", rep("---:", 5)))

# T8 single predictors ---------------------------------------------------------------
sp <- coefs_df(R$single_predictor)
rows <- bind_rows(lapply(list(c("pooled (region-adjusted)", "Pooled (study area as covariate)"), c("Honolulu", "Honolulu"), c("Ewa", "ʻEwa")), function(s) bind_rows(lapply(seq_along(PRED), function(k) {
  g <- sp %>% filter(scope == s[1], predictor == PRED[k])
  tibble(Model = if (k == 1) s[2] else "", Predictor = NAME[[PRED[k]]],
         `β per SD ± SE (°C)` = bold(sprintf("%s ± %s", s2(g$beta_z), f2(g$se_z)), g$lrt_p < 0.05),
         `Effect in natural units (°C)` = sprintf("%s ± %s %s", s2(g$beta_natural, 3), f2(g$se_natural, 3), g$natural_unit),
         `p (LRT)` = p_fmt(g$lrt_p), `R²m` = f2(g$r2m)) }))))
md_table(rows, "Single-predictor mixed models (random intercept per site). β is the change in ΔT per one standard deviation of the predictor (± standard error, REML); p is from the likelihood-ratio test against the null model (ML). Bold: p < 0.05.",
         "single", align = c("---", "---", rep("---:", 4)))

# T9 full pooled ---------------------------------------------------------------------
coef_rows <- function(cl) bind_rows(lapply(cl, function(cc) {
  name <- if (cc$term == "(Intercept)") "Intercept (ʻEwa, all predictors at their mean)" else if (cc$term == "regionHonolulu") "Study area: Honolulu (vs ʻEwa)" else NAME[[sub("_z$", "", cc$term)]]
  tibble(Term = name, `Estimate per SD (°C)` = bold(s2(cc$estimate, 3), cc$p < 0.05), SE = f2(cc$se, 3),
         `95% CI` = sprintf("[%s, %s]", s2(cc$lo), s2(cc$hi)), p = p_fmt(cc$p),
         `Effect in natural units` = if (grepl("_z$", cc$term)) sprintf("%s ± %s°C %s", s2(cc$estimate_natural, 3), f2(cc$se_natural, 3), cc$natural_unit) else "") }))
f <- R$full_pooled
md_table(coef_rows(f$coefs), sprintf("Full pooled mixed model: ΔT ~ study area + eight standardised predictors, random intercept per site (REML estimates; n = %d sensor-nights, %d sites; p-values from Satterthwaite t tests). Marginal R² = %s, conditional R² = %s; between-site variance %s °C², residual variance %s °C². Likelihood-ratio test of the eight predictors given study area: χ²(%d) = %.1f, p %s. Bold: p < 0.05.",
  f$n_obs, f$n_sites, f2(f$r2m), f2(f$r2c), f2(f$var_site, 3), f2(f$var_resid, 3), f$lrt_morphology_given_region[[2]], f$lrt_morphology_given_region[[1]], p_fmt(f$lrt_morphology_given_region[[3]])),
  "full_pooled", align = c("---", "---:", "---:", "---:", "---:", "---"))

# T10 district full models ----------------------------------------------------------------
rows <- bind_rows(lapply(c("(Intercept)", PRED), function(p) {
  r <- tibble(Term = if (p == "(Intercept)") "Intercept" else NAME[[p]])
  for (reg in REGIONS) { lab <- REG[[reg]]
    cs <- coefs_df(R$regional[[reg]]$full$coefs); cc <- cs[cs$term == (if (p == "(Intercept)") "(Intercept)" else paste0(p, "_z")), ]
    r[[sprintf("%s: β per SD (°C)", lab)]] <- bold(sprintf("%s ± %s", s2(cc$estimate, 3), f2(cc$se, 3)), cc$p < 0.05)
    r[[sprintf("%s: p", lab)]] <- p_fmt(cc$p)
    r[[sprintf("%s: natural units", lab)]] <- if (p == "(Intercept)") "" else sprintf("%s %s", s2(cc$estimate_natural, 3), cc$natural_unit) }
  r }))
h <- R$regional$Honolulu$full; e <- R$regional$Ewa$full
md_table(rows, sprintf("District full mixed models (eight standardised predictors, random intercept per site, REML; Satterthwaite p-values). Honolulu: R²m = %s, R²c = %s, LRT vs null χ²(8) = %.1f, p %s. ʻEwa: R²m = %s, R²c = %s, χ²(8) = %.1f, p %s. Standardisation uses the pooled standard deviations so that coefficients are comparable between areas. Bold: p < 0.05.",
  f2(h$r2m), f2(h$r2c), h$lrt_vs_null[[1]], p_fmt(h$lrt_vs_null[[3]]), f2(e$r2m), f2(e$r2c), e$lrt_vs_null[[1]], p_fmt(e$lrt_vs_null[[3]])),
  "regional", align = c("---", "---:", "---:", "---", "---:", "---:", "---"))

# T11 best subsets --------------------------------------------------------------------------
pretty_pred <- function(s) { s <- gsub("imperv", "impervious", s); s <- gsub("coast_km", "coast distance", s); s <- gsub("bldg", "footprint", s); gsub("svf_point", "sky view", s) }
rows <- bind_rows(lapply(REGIONS, function(reg) { b <- R$regional[[reg]]$best_subset
  bind_rows(lapply(seq_along(b$top5), function(k) { m <- b$top5[[k]]
    tibble(`Study area` = if (k == 1) REG[[reg]] else "", Rank = k, Predictors = pretty_pred(m$predictors), k = as.integer(m$k),
           `ΔAIC` = f2(m$delta_aic, 1), `Akaike weight` = f2(m$weight), `R²m` = f2(m$r2m)) })) }))
md_table(rows, "The five best subsets of at most three predictors in each study area, ranked by maximum-likelihood AIC among all 93 candidate models (null model included). Akaike weights are relative to the whole candidate set.",
         "subsets", align = c("---", "---:", "---", "---:", "---:", "---:", "---:"))
rows <- bind_rows(lapply(PRED, function(p) tibble(Predictor = NAME[[p]], `Importance Honolulu` = f2(R$regional$Honolulu$best_subset$importance[[p]]),
                                                  `Importance ʻEwa` = f2(R$regional$Ewa$best_subset$importance[[p]]))))
md_table(rows, "Relative importance of each predictor: sum of Akaike weights of all candidate models (≤ 3 predictors) that contain it, by study area.", "importance", align = c("---", "---:", "---:"))

# T12 interactions --------------------------------------------------------------------------
rows <- bind_rows(lapply(R$interactions$single, function(cc) tibble(Predictor = NAME[[cc$predictor]], `Slope ʻEwa (per SD)` = s2(cc$ewa_slope_z), `Slope Honolulu (per SD)` = s2(cc$hnl_slope_z),
  `Difference ± SE` = sprintf("%s ± %s", s2(cc$interaction_coef), f2(cc$se)), `LRT χ²(1)` = f2(cc$lrt_chi2, 1), p = bold(p_fmt(cc$p), cc$p < 0.05))))
j <- R$interactions$joint
md_table(rows, sprintf("Study-area × predictor interactions, each added on its own to the full pooled model and tested by likelihood-ratio test (ML). Slopes are from the model with that interaction (REML). Joint test of all eight interactions: χ²(%d) = %.1f, p %s; the interaction model has R²m = %s.", j$df, j$lrt_chi2, p_fmt(j$p), f2(j$r2m)),
         "interactions", align = c("---", rep("---:", 5)))

# T13 model comparison ------------------------------------------------------------------------
mc <- coefs_df(R$model_comparison) %>%
  mutate(Model = gsub("region", "study area", gsub("region x 8 predictors", "study area × 8 predictors", model, fixed = TRUE), fixed = TRUE),
         Parameters = k, AIC = f2(aic, 1), `ΔAIC` = f2(delta_aic, 1)) %>% select(Model, Parameters, AIC, `ΔAIC`)
md_table(mc, "Maximum-likelihood AIC of the candidate model structures (all with the same 400 sensor-nights; the parameter count includes the variance components).", "comparison", align = c("---", "---:", "---:", "---:"))

# T14 thresholds ------------------------------------------------------------------------------
th <- coefs_df(R$threshold_sensitivity) %>%
  filter((cloud_max == 25 & min_frac == 0.75) | (wind_max == 10 & min_frac == 0.75) | wind_max <= 0 | wind_max == 999 | min_frac != 0.75) %>%
  distinct(wind_max, cloud_max, min_frac, .keep_all = TRUE) %>%
  mutate(`Selection rule` = case_when(wind_max == 999 ~ "no weather filter (all complete nights)", wind_max == -1 ~ "windy or cloudy nights only (complement)",
                                      TRUE ~ sprintf("wind < %dkm h⁻¹, cloud < %d%%, ≥ %d%% of bins", as.integer(wind_max), as.integer(cloud_max), as.integer(round(min_frac * 100)))),
         `Nights Hon/ʻEwa` = sprintf("%d / %d", nights_honolulu, nights_ewa), `Sensor-nights` = sensor_nights,
         `Imperv. β` = sprintf("%s ± %s", s2(beta_imperv_z), f2(se_imperv)), `Height β` = sprintf("%s ± %s", s2(beta_height_z), f2(se_height)),
         `R²m (full)` = f2(r2m_full), `SD sites (°C)` = f2(site_sd), `r vs base` = f2(corr_with_baseline),
         ord = ifelse(wind_max == 999, 900, ifelse(wind_max == -1, 950, wind_max)) + cloud_max / 100 + min_frac / 1000) %>% arrange(ord) %>%
  select(`Selection rule`, `Nights Hon/ʻEwa`, `Sensor-nights`, `Imperv. β`, `Height β`, `R²m (full)`, `SD sites (°C)`, `r vs base`)
md_table(th, "Sensitivity of the results to the night-selection thresholds. Nights were re-selected from the complete record under each rule and the pooled models refitted. β: effect per pooled SD (± SE) in the parsimonious model ΔT ~ study area + impervious + height; R²m (full): marginal R² of the full pooled model; SD sites: standard deviation of site-mean ΔT; r vs base: correlation of site means with those of the adopted rule (wind < 10km h⁻¹, cloud < 25%, ≥ 75% of bins).",
         "thresholds", align = c("---", rep("---:", 7)))

# T15 LOSO --------------------------------------------------------------------------------------
cv <- coefs_df(R$loso_cv) %>%
  mutate(Model = recode(model, "full" = "Pooled: study area + 8 predictors", "parsimonious" = "Pooled: study area + impervious + height", "region only" = "Pooled: study area only",
                        "Honolulu full" = "Honolulu: 8 predictors", "Ewa full" = "ʻEwa: 8 predictors",
                        "Honolulu best subset (imperv + coast_km)" = "Honolulu: impervious + coast distance", "Ewa best subset (imperv + height + aspect)" = "ʻEwa: impervious + height + aspect ratio"),
         `RMSE (°C)` = f2(rmse), `MAE (°C)` = f2(mae), r = f2(r), `R²cv` = f2(r2_cv)) %>% select(Model, `RMSE (°C)`, `MAE (°C)`, r, `R²cv`)
md_table(cv, "Leave-one-site-out cross-validation of site-mean ΔT (each site predicted from the fixed effects of a model fitted without it). R²cv is 1 − SSE/SST over the held-out sites (negative values mean worse than the mean).", "loso", align = c("---", rep("---:", 4)))

# T16 shared nights ---------------------------------------------------------------------------
sh <- coefs_df(R$shared_nights$rows) %>%
  transmute(Night = trimws(format(as.Date(night_date), "%e %b %Y")), `Study area` = REG[region], Sensors = n_sensors, `Network median (°C)` = f2(T_median_night, 1),
            `10th–90th percentile of sites (°C)` = sprintf("%s – %s", f2(T_p10, 1), f2(T_p90, 1)), `Coolest – warmest site (°C)` = sprintf("%s – %s", f2(T_min_site, 1), f2(T_max_site, 1)),
            `Wind (km h⁻¹) / cloud (%)` = sprintf("%.1f / %.0f", wind, cloud))
md_table(sh, "Absolute night-mean (18:00–06:00 HST) air temperatures of the two networks on the two calm, clear nights they share.", "shared", align = c("---", "---", rep("---:", 5)))

# Appendix: site tables ----------------------------------------------------------------------
site <- read_csv(file.path(TAB, "site_summary.csv"), show_col_types = FALSE) %>%
  inner_join(read_csv(file.path(TAB, "blups.csv"), show_col_types = FALSE) %>% select(sensor_id, blup, latitude, longitude), by = "sensor_id") %>% arrange(region, sensor_id)
ta <- site %>% transmute(Area = REG[region], Site = sensor_id, `Lat (°N)` = f2(latitude, 4), `Lon (°E)` = f2(longitude, 4), Nights = n_nights, `ΔT (°C)` = s2(dT_site), `SD (°C)` = s2(dT_sd), `BLUP (°C)` = s2(blup))
md_table(ta, "All 74 sites: coordinates (WGS84), number of analysed nights, site-mean ΔT with its night-to-night standard deviation, and the random intercept (BLUP) of the full pooled model.", "sites", align = c("---", "---", rep("---:", 6)))
tb <- site %>% transmute(Area = REG[region], Site = sensor_id, `Imperv. %` = f2(imperv, 1), `Tree %` = f2(tree, 1), `Bldg %` = f2(bldg, 1), `Water %` = f2(water, 1), `Height m` = f2(height, 1), SVF = f2(svf_point), `H/W` = f2(aspect), `Coast km` = f2(coast_km))
md_table(tb, sprintf("All 74 sites: the eight descriptors at the %d m source-area radius (sky view factor and coastal distance are point properties).", RAD), "sites_predictors", align = c("---", "---", rep("---:", 8)))

# Software table (Appendix C) --------------------------------------------------------------------------
ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "")
sw <- tribble(~Package, ~Version, ~Role, ~`Functions used`, ~Reference,
  "R", sub("R version ([0-9.]+).*", "\\1", R.version.string), "the language: arithmetic, median(), cor(), lm(), wilcox.test(), cor.test()", "", "[@rcore2024]",
  "dplyr, tidyr, tibble, purrr, readr", paste(ver("dplyr"), ver("tidyr"), ver("readr"), sep = " / "), "grouped summaries, joins and reshaping of the record; typed CSV input/output", "group_by(), summarise(), inner_join(), pivot_wider(), read_csv()", "[@wickham2019]",
  "lubridate", ver("lubridate"), "rounding of the loggers' time stamps to the nearest 30 minutes", "round_date(), ymd_hms()", "[@grolemund2011]",
  "jsonlite", ver("jsonlite"), "the results file read by the tables and figures", "write_json(), fromJSON()", "[@ooms2014]",
  "sf", ver("sf"), "vector geometry: projection, buffers, footprint clipping, sight lines for the sky view factor, shoreline distance", "st_transform(), st_buffer(), st_intersects(), st_intersection(), st_distance()", "[@pebesma2018]",
  "terra", ver("terra"), "the 1 m C-CAP masks as virtual mosaics; pixel extraction inside circles; aggregation for base maps", "vrt(), extract(), crop(), aggregate()", "[@hijmans2024]",
  "lme4", ver("lme4"), "the mixed-effects models (random intercept per site; REML and ML fits; BLUPs)", "lmer(), fixef(), ranef(), VarCorr(), logLik(), AIC(), predict()", "[@bates2015]",
  "lmerTest", ver("lmerTest"), "Satterthwaite degrees of freedom and t tests for the coefficients", "lmer(), summary()", "[@kuznetsova2017]",
  "performance", ver("performance"), "check of the marginal and conditional R² and ICC computed from the variance components", "r2_nakagawa(), icc()", "[@ludecke2021]",
  "car", ver("car"), "variance inflation factors", "vif()", "[@fox2019]",
  "ggplot2, patchwork, scales, cowplot", paste(ver("ggplot2"), ver("patchwork"), sep = " / "), "all figures; multi-panel layout; color scales", "ggplot(), geom_sf(), geom_raster(), coord_sf(), wrap_plots()", "[@wickham2016; @pedersen2024]",
  "rmarkdown, knitr", paste(ver("rmarkdown"), ver("knitr"), sep = " / "), "the analysis walkthrough document", "render()", "",
  "Python replication: pandas, NumPy, SciPy", "3.0.2 / 2.4.4 / 1.17.1", "record processing and statistics in the independent replication", "", "[@mckinney2010; @harris2020; @virtanen2020]",
  "Python replication: statsmodels", "0.15.0", "mixed models (MixedLM) in the replication", "", "[@seabold2010]",
  "Python replication: rasterio, Shapely, GeoPandas", "1.4.4 / 2.1.2 / 1.1.4", "rasters and vector geometry in the replication", "", "[@gillies2013; @gillies2007]")
md_table(sw, "Software used for the analysis, with the role of each package and the functions called. Versions are those of the run that produced this document (`results/sessionInfo.txt`).", "software", align = c("---", "---", "---", "---", "---"))

# Replication check (Appendix C) ------------------------------------------------------------------
rcf <- file.path(ROOT, "replication", "replication_check.csv")
if (file.exists(rcf)) {
  rp <- read_csv(rcf, show_col_types = FALSE) %>% filter(!is.na(Python))
  pick <- rp %>% filter(grepl("^Full pooled: .*_z estimate|^Full pooled: R2|^Null model ICC|best subset: R2m|Interactions: joint|^BLUP SD|^LOSO .*: r$|Scale rule|Full pooled: LRT vs null", quantity)) %>%
    transmute(Quantity = quantity, R = sprintf("%.3f", R), Python = sprintf("%.3f", Python), Difference = sprintf("%+.3f", difference))
  md_table(pick, sprintf("Agreement between the R analysis (lme4 / lmerTest) and the independent Python replication (statsmodels MixedLM) on the same inputs: a selection of the %d paired quantities in `replication_check.csv` (which also holds 25 maximum-difference checks over the 74 sites). Differences larger than 0.01 are optimizer tolerance (ΔAIC, joint χ²) or a documented implementation choice (one unstable leave-one-site-out fold of the ʻEwa full model in statsmodels).", nrow(rp)),
           "replication", align = c("---", "---:", "---:", "---:"))
}

# Terrain check (Section 4.7; 06_elevation_check.R) ----------------------------------------------------
ecf <- file.path(TAB, "elevation_models.csv")
if (file.exists(ecf)) {
  em <- read_csv(ecf, show_col_types = FALSE)
  E <- fromJSON(file.path(OUT, "elevation_check.json"), simplifyVector = FALSE)
  est <- function(b, se, p) ifelse(is.na(b), "", sprintf("%s ± %s (p %s)", s2(b), f2(se), ifelse(p < 0.001, "< 0.001", paste0("= ", p_fmt(p)))))
  keep <- list(Honolulu = c("coast only", "elevation only", "coast + elevation", "best subset (imperv + coast)", "best subset (imperv + coast) + elevation",
                            "full (eight descriptors)", "full (eight descriptors) + elevation"),
               Ewa = c("coast only", "elevation only", "best subset (imperv + height + aspect)", "best subset (imperv + height + aspect) + elevation",
                       "full (eight descriptors)", "full (eight descriptors) + elevation"),
               pooled = c("elevation only", "parsimonious (district + impervious + height)", "parsimonious (district + impervious + height) + elevation",
                          "full (district + eight descriptors)", "full (district + eight descriptors) + elevation"))
  rows <- bind_rows(lapply(names(keep), function(sc) {
    x <- em %>% filter(scope == sc, model %in% keep[[sc]]) %>% mutate(model = factor(model, levels = keep[[sc]])) %>% arrange(model)
    x <- x %>% mutate(d_aic = case_when(model == "elevation only" ~ delta_aic_vs_null, model == "coast + elevation" ~ aic - em$aic[em$scope == sc & em$model == "coast only"],
                                        !is.na(delta_aic_vs_without) ~ delta_aic_vs_without, TRUE ~ NA_real_),
                      p_el = case_when(model == "elevation only" ~ lrt_p, model == "coast + elevation" ~ lrt_elev_given_coast_p,
                                       !is.na(lrt_elev_given_model_p) ~ lrt_elev_given_model_p, TRUE ~ NA_real_))
    tibble(`Study area` = c(if (sc == "pooled") "Pooled (with study area)" else REG[sc], rep("", nrow(x) - 1)),
           Model = gsub("\\bimperv\\b", "impervious", gsub("district", "study area", as.character(x$model))),
           `Coast (°C per km)` = est(x$coast_per_km, x$coast_se, x$coast_p),
           `Elevation (°C per 10m)` = est(x$elev_per_10m, x$elev_se, x$elev_p),
           `Marginal R²` = f2(x$r2m), `ΔAIC (adding elevation)` = ifelse(is.na(x$d_aic), "", s2(x$d_aic, 1)),
           `LRT p (adding elevation)` = ifelse(is.na(x$p_el), "", p_fmt(x$p_el)))
  }))
  md_table(rows, sprintf("Terrain check: the models of Sections 4.5–4.7 refitted with the bare-earth elevation of each site (USGS 3DEP 1m lidar DEM) as an additional standardised term. Coefficients are REML estimates in natural units (± SE; Satterthwaite p-values); ΔAIC and the likelihood-ratio p-value refer to adding elevation to the same model without it (ML fits; for the elevation-only rows, to the null model). R²m is that of the REML fit, so the best-subset rows read 0.34 and 0.30 where Table {{tbl:subsets}}, which ranks ML fits, gives 0.36 and 0.31. With elevation in the candidate pool (%d models of ≤ 3 terms per study area), the best subsets of Table {{tbl:subsets}} remain the best; the summed Akaike weight of models containing elevation is %s in Honolulu and %s in ʻEwa, against %s and %s for distance to the coast.",
                         E$subsets$Honolulu$n_models, f2(E$subsets$Honolulu$elevation_importance), f2(E$subsets$Ewa$elevation_importance),
                         f2(E$subsets$Honolulu$coast_importance), f2(E$subsets$Ewa$coast_importance)),
           "elevation", align = c("---", "---", "---:", "---:", "---:", "---:", "---:"))
}
message("done")
