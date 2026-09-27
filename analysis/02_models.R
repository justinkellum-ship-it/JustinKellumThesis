#!/usr/bin/env Rscript
# 02_models.R --------------------------------------------------------------------
# Mixed-effects analysis of nocturnal canopy-layer temperature departures
# (dT = logger - district network median, 18:00-06:00 HST means on calm, clear
# nights) against two- and three-dimensional urban-form predictors at 50, 100 and
# 200 m.
#
# Unit of analysis: the sensor-night (400 rows: 36 Honolulu sites x 7 nights and
# 38 ʻEwa sites x 4 nights).  Every model has a random intercept per site,
# (1 | sensor_id).  Fixed effects are tested with likelihood-ratio tests on ML
# fits; estimates, standard errors and variance components are reported from
# REML fits (Zuur et al. 2009); p-values of individual coefficients use
# Satterthwaite degrees of freedom (lmerTest).
#
# Steps (mirroring the Methods chapter)
#   1  descriptive statistics of the predictors by district and radius
#   2  scale selection: pooled and district full models per radius (ML AIC, R2m);
#      a forward-stepwise OLS on site means (AIC) as a secondary check
#   3  correlation and variance-inflation screen at the adopted radius
#   4  null models and intraclass correlation
#   5  single-predictor models (district-adjusted, and per district)
#   6  full eight-predictor model with district (pooled)
#   7  district models, best-subset (<= 3 predictors) AIC ranking, Akaike weights
#   8  district x predictor interaction tests (transferability)
#   9  night fixed-effect check (the reference removes the night mean)
#  10  calm/clear threshold sensitivity (night selection re-run from the record)
#  11  model comparison table (ML AIC)
#  12  BLUPs, leave-one-site-out cross-validation, calm-core refit
#  13  absolute temperature comparison on the two shared nights
#
# Outputs: results/thesis_results.json and results/tables/*.csv
#
# Packages: lme4 (lmer: the mixed models), lmerTest (Satterthwaite df and p-values
#           for the coefficients), performance (Nakagawa R2 and ICC), car (VIF),
#           dplyr/tidyr/purrr (data handling), jsonlite (results file).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")
# the district term of the pooled models: the first district against the second (reference level)
stopifnot(length(REGIONS) == 2)     # the district comparisons below are written for two districts
REGION_TERM <- paste0("region", REGIONS[1])
source(file.path(.here, "00_prepare_inputs.R"), encoding = "UTF-8")          # load_binned(), night_inventory(), sensor_nights()
suppressPackageStartupMessages({ library(lme4); library(lmerTest); library(performance); library(car) })

RES <- list()

# ---- model helpers ---------------------------------------------------------------
fit_lmm <- function(df, rhs, reml = TRUE) {
  # random-intercept model  dT_night ~ rhs + (1 | sensor_id)
  # bobyqa with a generous evaluation budget; singular fits (site variance -> 0)
  # are legitimate outcomes for the sensitivity models and are kept
  f <- as.formula(paste0("dT_night ~ ", rhs, " + (1 | sensor_id)"))
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5),
                      check.conv.singular = "ignore")
  m <- suppressWarnings(suppressMessages(lmerTest::lmer(f, data = df, REML = reml, control = ctrl)))
  if (length(m@optinfo$conv$lme4$messages)) {          # fall back to Nelder-Mead if bobyqa complained
    m2 <- suppressWarnings(suppressMessages(lmerTest::lmer(f, data = df, REML = reml,
            control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1e5), check.conv.singular = "ignore"))))
    if (!length(m2@optinfo$conv$lme4$messages)) m <- m2
  }
  m
}
n_fixed <- function(m) length(fixef(m))
n_params <- function(m) n_fixed(m) + 2L                     # + site variance + residual variance
aic_ml <- function(m) as.numeric(AIC(m))                    # ML fits only: -2 logLik + 2 (k + 2)
llf <- function(m) as.numeric(logLik(m))
lrt <- function(m_full, m_null) {                           # both ML fits
  stat <- 2 * (llf(m_full) - llf(m_null)); df <- n_fixed(m_full) - n_fixed(m_null)
  list(stat = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE))
}
var_comp <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  list(var_site = vc$vcov[vc$grp == "sensor_id"], var_resid = vc$vcov[vc$grp == "Residual"])
}
r2_nakagawa <- function(m) {
  # Nakagawa & Schielzeth (2013): R2m = var_fixed / (var_fixed + var_site + var_resid),
  # R2c = (var_fixed + var_site) / total, with var_fixed the variance of the
  # fixed-effect predictions (performance::r2_nakagawa gives the same numbers)
  vc <- var_comp(m)
  var_f <- var(as.numeric(model.matrix(m) %*% fixef(m)))
  tot <- var_f + vc$var_site + vc$var_resid
  list(r2m = var_f / tot, r2c = (var_f + vc$var_site) / tot, var_fixed = var_f, var_site = vc$var_site,
       var_resid = vc$var_resid, icc_cond = vc$var_site / (vc$var_site + vc$var_resid))
}
coef_table <- function(m, sds = NULL) {
  # fixed effects with Satterthwaite df, t and p (lmerTest) and t-based 95 % CI;
  # z-scored slopes are also expressed in natural units (UNIT in helpers.R)
  s <- summary(m)$coefficients
  t <- tibble(term = rownames(s), estimate = s[, "Estimate"], se = s[, "Std. Error"], df = s[, "df"],
              t = s[, "t value"], p = s[, "Pr(>|t|)"]) %>%
    mutate(lo = estimate - qt(0.975, df) * se, hi = estimate + qt(0.975, df) * se)
  if (!is.null(sds)) {
    t <- t %>% mutate(natural_unit = "", estimate_natural = NA_real_, se_natural = NA_real_)
    for (i in seq_len(nrow(t))) {
      base <- sub("_z$", "", t$term[i])
      if (grepl("_z$", t$term[i]) && base %in% names(sds)) {
        f <- UNIT[[base]][[1]] / sds[[base]]
        t$natural_unit[i] <- UNIT[[base]][[2]]
        t$estimate_natural[i] <- t$estimate[i] * f
        t$se_natural[i] <- t$se[i] * abs(f)
      }
    }
  }
  t
}
rhs_full <- function(extra = "region") paste(c(extra, paste0(PRED, "_z")), collapse = " + ")
pred_rows <- function(t) lapply(seq_len(nrow(t)), function(i) as.list(t[i, ]))

load_radius <- function(radius) {
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
  pr <- read_csv(file.path(DER, "site_predictors_multiscale.csv"), show_col_types = FALSE)
  p <- pr %>% filter(radius_m == radius) %>% select(-radius_m, -region)
  d <- sn %>% inner_join(p, by = "sensor_id") %>%
    mutate(region = factor(region, levels = rev(REGIONS)))     # the second district is the reference level
  list(d = d, pr = pr)
}
site_means <- function(d) d %>% group_by(region, sensor_id) %>%
  summarise(dT_site = mean(dT_night), dT_sd = sd(dT_night), n_nights = n(),
            across(all_of(PRED), first), .groups = "drop")

# ---- step 1: descriptives ------------------------------------------------------------
step1_descriptives <- function(pr) {
  t <- pr %>% pivot_longer(all_of(PRED), names_to = "predictor") %>%
    group_by(region, radius_m, predictor) %>%
    summarise(mean = mean(value), sd = sd(value), min = min(value), median = median(value), max = max(value), .groups = "drop") %>%
    mutate(predictor = factor(predictor, levels = PRED)) %>% arrange(region, radius_m, predictor) %>%
    mutate(predictor = as.character(predictor))
  write_table(t, "predictors_by_region")
  RES$predictors_by_region <<- setNames(lapply(c(outer(REGIONS, RADII, paste, sep = "_")), function(k) {
    r <- strsplit(k, "_")[[1]]
    setNames(lapply(PRED, function(p) { x <- t %>% filter(region == r[1], radius_m == as.numeric(r[2]), predictor == p)
                                        list(mean = x$mean, sd = x$sd) }), PRED)
  }), c(outer(REGIONS, RADII, paste, sep = "_")))
  # district difference at 100 m: Mann-Whitney (Wilcoxon rank-sum, normal approximation)
  p100 <- pr %>% filter(radius_m == 100)
  RES$region_difference_100m <<- setNames(lapply(PRED, function(p) {
    a <- p100[[p]][p100$region == REGIONS[1]]; b <- p100[[p]][p100$region == REGIONS[2]]
    list(honolulu_median = median(a), ewa_median = median(b),
         p = suppressWarnings(wilcox.test(a, b, exact = FALSE, correct = TRUE))$p.value)
  }), PRED)
}

# ---- step 2: scale selection -----------------------------------------------------------
step2_scale <- function() {
  rows <- list(); single <- list()
  for (R in RADII) {
    d <- load_radius(R)$d
    z <- zscore_frame(d, PRED); dz <- z$df
    m <- fit_lmm(dz, rhs_full("region"), reml = FALSE); r2 <- r2_nakagawa(m)
    rows[[length(rows) + 1]] <- tibble(scope = "pooled", radius_m = R, model = "region + 8 predictors", aic = aic_ml(m),
                                       llf = llf(m), r2m = r2$r2m, r2c = r2$r2c, var_site = r2$var_site)
    for (reg in REGIONS) {
      dr <- zscore_frame(d %>% filter(region == reg), PRED)$df
      mr <- fit_lmm(dr, rhs_full(NULL), reml = FALSE); r2r <- r2_nakagawa(mr)
      rows[[length(rows) + 1]] <- tibble(scope = reg, radius_m = R, model = "8 predictors", aic = aic_ml(mr),
                                         llf = llf(mr), r2m = r2r$r2m, r2c = r2r$r2c, var_site = r2r$var_site)
      null <- fit_lmm(dr, "1", reml = FALSE)
      for (p in PRED) {                     # at which scale does each predictor work best?
        ms <- fit_lmm(dr, paste0(p, "_z"), reml = FALSE); l <- lrt(ms, null)
        single[[length(single) + 1]] <- tibble(region = reg, radius_m = R, predictor = p,
          beta_z = unname(fixef(ms)[paste0(p, "_z")]), se = unname(sqrt(diag(vcov(ms)))[2]),
          lrt_p = l$p, r2m = r2_nakagawa(ms)$r2m, aic = aic_ml(ms))
      }
    }
  }
  t <- bind_rows(rows) %>% group_by(scope) %>% mutate(delta_aic = aic - min(aic)) %>% ungroup()
  write_table(t, "scale_selection_lmm"); write_table(bind_rows(single), "scale_single_predictor")
  RES$scale_selection <<- pred_rows(t)

  # secondary check: forward stepwise OLS on site means (AIC), per district and radius
  steps <- list()
  for (R in RADII) {
    sm_ <- site_means(load_radius(R)$d)
    for (reg in REGIONS) {
      g <- zscore_frame(sm_ %>% filter(region == reg), PRED)$df
      remaining <- paste0(PRED, "_z"); chosen <- character(0)
      cur_aic <- AIC(lm(dT_site ~ 1, data = g)); k <- 0
      steps[[length(steps) + 1]] <- tibble(region = reg, radius_m = R, step = 0, added = "(null)", aic = cur_aic, r2 = 0, adj_r2 = 0)
      while (length(remaining)) {
        cand <- lapply(remaining, function(p) { f <- lm(reformulate(c(chosen, p), "dT_site"), data = g)
                                                list(p = p, aic = AIC(f), r2 = summary(f)$r.squared, adj = summary(f)$adj.r.squared) })
        best <- cand[[which.min(sapply(cand, `[[`, "aic"))]]
        if (best$aic < cur_aic - 1e-9) {
          chosen <- c(chosen, best$p); remaining <- setdiff(remaining, best$p); cur_aic <- best$aic; k <- k + 1
          steps[[length(steps) + 1]] <- tibble(region = reg, radius_m = R, step = k, added = sub("_z$", "", best$p),
                                               aic = best$aic, r2 = best$r2, adj_r2 = best$adj)
        } else break
      }
    }
  }
  st <- bind_rows(steps); write_table(st, "stepwise_ols_site_means")
  RES$stepwise_final <<- pred_rows(st %>% group_by(region, radius_m) %>% slice_max(step, n = 1) %>% ungroup() %>%
                                     select(region, radius_m, step, aic, r2, adj_r2))
}

# ---- step 3: collinearity ----------------------------------------------------------------
step3_collinearity <- function(d) {
  out <- list()
  for (scope in c("pooled", REGIONS)) {
    g <- if (scope == "pooled") d else d %>% filter(region == scope)
    s <- g %>% distinct(sensor_id, .keep_all = TRUE) %>% select(all_of(PRED))
    corr <- cor(s); write.csv(round(corr, 6), file.path(TAB, paste0("correlation_", scope, ".csv")))
    sz <- as.data.frame(scale(s)); sz$y <- rnorm(nrow(sz))          # VIF depends on the predictors only
    v <- car::vif(lm(y ~ ., data = sz))
    out[[scope]] <- list(vif = as.list(v[PRED]),
                         max_abs_corr = setNames(lapply(PRED, function(p) max(abs(corr[p, setdiff(PRED, p)]))), PRED))
  }
  vif_t <- as.data.frame(sapply(out, function(o) unlist(o$vif))); write.csv(vif_t, file.path(TAB, "vif.csv"))
  RES$collinearity <<- out
}

# ---- step 4: null models and ICC -----------------------------------------------------------
step4_null <- function(d) {
  out <- list()
  for (scope in c("pooled", REGIONS)) {
    g <- if (scope == "pooled") d else d %>% filter(region == scope)
    m <- fit_lmm(g, "1", reml = TRUE); vc <- var_comp(m)
    out[[scope]] <- list(intercept = unname(fixef(m)[1]), var_site = vc$var_site, var_resid = vc$var_resid,
                         icc = vc$var_site / (vc$var_site + vc$var_resid), n_sites = n_distinct(g$sensor_id), n_obs = nrow(g),
                         site_sd = sqrt(vc$var_site), resid_sd = sqrt(vc$var_resid))
  }
  m <- fit_lmm(d, "region", reml = TRUE); vc <- var_comp(m); ct <- coef_table(m)
  out$pooled_region <- list(var_site = vc$var_site, var_resid = vc$var_resid, icc = vc$var_site / (vc$var_site + vc$var_resid),
                            region_coef = ct$estimate[ct$term == REGION_TERM], region_p = ct$p[ct$term == REGION_TERM])
  nt <- bind_rows(lapply(names(out), function(k) as_tibble(out[[k]]) %>% mutate(scope = k))) %>% select(scope, everything())
  write.csv(nt, file.path(TAB, "null_models.csv"), row.names = FALSE)
  RES$null <<- out
}

# ---- step 5: single predictors ---------------------------------------------------------------
step5_single <- function(d, sds) {
  one <- function(g, scope, base_rhs, null_ml) {
    bind_rows(lapply(PRED, function(p) {
      rhs <- if (nzchar(base_rhs)) paste(base_rhs, paste0(p, "_z"), sep = " + ") else paste0(p, "_z")
      m_ml <- fit_lmm(g, rhs, reml = FALSE); m <- fit_lmm(g, rhs, reml = TRUE); l <- lrt(m_ml, null_ml)
      ct <- coef_table(m, sds); r <- ct[ct$term == paste0(p, "_z"), ]
      tibble(scope = scope, predictor = p, beta_z = r$estimate, se_z = r$se, lrt_chi2 = l$stat, lrt_p = l$p,
             r2m = r2_nakagawa(m)$r2m, beta_natural = r$estimate_natural, se_natural = r$se_natural,
             natural_unit = r$natural_unit, aic = aic_ml(m_ml))
    }))
  }
  t <- bind_rows(one(d, "pooled (region-adjusted)", "region", fit_lmm(d, "region", reml = FALSE)),
                 bind_rows(lapply(REGIONS, function(reg) { g <- d %>% filter(region == reg)
                                                          one(g, reg, "", fit_lmm(g, "1", reml = FALSE)) })))
  write_table(t, "single_predictor_models")
  RES$single_predictor <<- pred_rows(t)
}

# ---- step 6: full pooled model ---------------------------------------------------------------
step6_full <- function(d, sds) {
  f <- rhs_full("region")
  m <- fit_lmm(d, f, reml = TRUE); m_ml <- fit_lmm(d, f, reml = FALSE)
  null_ml <- fit_lmm(d, "1", reml = FALSE); null_reg_ml <- fit_lmm(d, "region", reml = FALSE)
  t <- coef_table(m, sds); write_table(t, "full_model_pooled")
  r2 <- r2_nakagawa(m); l_all <- lrt(m_ml, null_ml); l_m <- lrt(m_ml, null_reg_ml)
  fe <- fixef(m)
  RES$full_pooled <<- c(list(coefs = pred_rows(t)), r2,
                        list(aic = aic_ml(m_ml), lrt_vs_null = unname(unlist(l_all)), lrt_morphology_given_region = unname(unlist(l_m)),
                             converged = !length(m@optinfo$conv$lme4$messages), n_obs = nrow(d), n_sites = n_distinct(d$sensor_id)))
  # district gap: the region coefficient is the gap at the pooled mean morphology; what is
  # predicted for each district's own average site, and for an ʻEwa-average site placed in Honolulu?
  prof <- d %>% distinct(sensor_id, .keep_all = TRUE) %>% group_by(region) %>% summarise(across(paste0(PRED, "_z"), mean))
  pr_of <- function(reg) setNames(as.numeric(prof[prof$region == reg, paste0(PRED, "_z")]), PRED)
  pred <- function(region, pf) unname(fe["(Intercept)"] + (if (region == REGIONS[1]) fe[REGION_TERM] else 0) +
                                        sum(fe[paste0(PRED, "_z")] * pf[PRED]))
  hnl <- pr_of(REGIONS[1]); ewa <- pr_of(REGIONS[2])
  RES$full_pooled$region_gap_at_pooled_mean <<- unname(fe[REGION_TERM])
  RES$full_pooled$pred_hnl_profile_in_hnl <<- pred(REGIONS[1], hnl)
  RES$full_pooled$pred_ewa_profile_in_ewa <<- pred(REGIONS[2], ewa)
  RES$full_pooled$pred_ewa_profile_in_hnl <<- pred(REGIONS[1], ewa)
  RES$full_pooled$morphology_contribution_hnl_minus_ewa <<- unname(sum(fe[paste0(PRED, "_z")] * (hnl - ewa)))
  # parsimonious pooled model (district + impervious + height) used by the sensitivity analyses
  mp <- fit_lmm(d, "region + imperv_z + height_z", reml = TRUE); mp_ml <- fit_lmm(d, "region + imperv_z + height_z", reml = FALSE)
  tp <- coef_table(mp, sds); write_table(tp, "parsimonious_model_pooled")
  RES$parsimonious_pooled <<- c(list(coefs = pred_rows(tp)), r2_nakagawa(mp),
                                list(aic = aic_ml(mp_ml), lrt_full_vs_parsimonious = unname(unlist(lrt(m_ml, mp_ml)))))
  list(m = m, m_ml = m_ml)
}

# ---- step 7: district models and best subsets --------------------------------------------------
step7_regional <- function(d, sds) {
  out <- list(); subsets <- list()
  for (reg in REGIONS) {
    g <- d %>% filter(region == reg)
    f <- rhs_full(NULL)
    m <- fit_lmm(g, f, reml = TRUE); m_ml <- fit_lmm(g, f, reml = FALSE); null_ml <- fit_lmm(g, "1", reml = FALSE)
    t <- coef_table(m, sds) %>% mutate(region = reg); write_table(t, paste0("full_model_", reg))
    l <- lrt(m_ml, null_ml)
    out[[reg]] <- list(full = c(list(coefs = pred_rows(t)), r2_nakagawa(m),
                                list(aic = aic_ml(m_ml), lrt_vs_null = unname(unlist(l)),
                                     converged = !length(m@optinfo$conv$lme4$messages),
                                     n_sites = n_distinct(g$sensor_id), n_obs = nrow(g))))
    # all subsets of <= 3 predictors, ML AIC and Akaike weights (what MuMIn::dredge would do)
    cands <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(PRED, k, simplify = FALSE)), recursive = FALSE))
    rr <- bind_rows(lapply(cands, function(cc) {
      rhs <- if (length(cc)) paste(paste0(cc, "_z"), collapse = " + ") else "1"
      mm <- fit_lmm(g, rhs, reml = FALSE)
      tibble(region = reg, predictors = if (length(cc)) paste(cc, collapse = " + ") else "(null)", k = length(cc),
             aic = aic_ml(mm), r2m = r2_nakagawa(mm)$r2m, llf = llf(mm))
    })) %>% mutate(delta_aic = aic - min(aic), weight = exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))) %>% arrange(aic)
    subsets[[reg]] <- rr
    importance <- setNames(lapply(PRED, function(p) sum(rr$weight[grepl(p, rr$predictors, fixed = TRUE)])), PRED)
    best <- rr[1, ]
    fb <- if (best$k > 0) paste(paste0(strsplit(best$predictors, " \\+ ")[[1]], "_z"), collapse = " + ") else "1"
    mb <- fit_lmm(g, fb, reml = TRUE)
    tb <- coef_table(mb, sds) %>% mutate(region = reg); write_table(tb, paste0("best_subset_model_", reg))
    r2b <- r2_nakagawa(mb); r2b$r2m <- NULL
    out[[reg]]$best_subset <- c(list(predictors = best$predictors, aic = best$aic, r2m = best$r2m, coefs = pred_rows(tb)),
                                r2b, list(importance = importance, top5 = pred_rows(head(rr, 5)), n_within_2 = sum(rr$delta_aic <= 2)))
  }
  write_table(bind_rows(subsets), "best_subsets")
  RES$regional <<- out
}

# ---- step 8: interactions ------------------------------------------------------------------------
step8_interactions <- function(d, m_ml_base) {
  base <- rhs_full("region")
  t <- bind_rows(lapply(PRED, function(pred) {
    pz <- paste0(pred, "_z"); key <- paste0(REGION_TERM, ":", pz)
    f <- paste0(base, " + region:", pz)
    m <- fit_lmm(d, f, reml = FALSE); l <- lrt(m, m_ml_base); mr <- fit_lmm(d, f, reml = TRUE)
    fe <- fixef(mr); se_key <- unname(sqrt(diag(vcov(mr)))[names(fe) == key])
    tibble(predictor = pred, lrt_chi2 = l$stat, df = l$df, p = l$p, interaction_coef = unname(fe[key]),
           se = se_key, ewa_slope_z = unname(fe[pz]), hnl_slope_z = unname(fe[pz] + fe[key]))
  }))
  f_all <- paste0(base, " + ", paste0("region:", PRED, "_z", collapse = " + "))
  m_all <- fit_lmm(d, f_all, reml = FALSE); l <- lrt(m_all, m_ml_base)
  write_table(t, "interaction_tests")
  RES$interactions <<- list(single = pred_rows(t), joint = list(lrt_chi2 = l$stat, df = l$df, p = l$p, aic = aic_ml(m_all),
                                                                r2m = r2_nakagawa(m_all)$r2m))
  m_all
}

# ---- step 9: night effect and night-to-night stability ---------------------------------------------
step9_night <- function(d, m_ml_base) {
  m <- fit_lmm(d %>% mutate(night_id = factor(night_id)), paste0(rhs_full("region"), " + night_id"), reml = FALSE)
  l <- lrt(m, m_ml_base)
  nm <- d %>% group_by(region, night_id) %>% summarise(mean = mean(dT_night), std = sd(dT_night), size = n(), .groups = "drop")
  write_table(nm, "night_means")
  RES$night_effect <<- list(lrt_chi2 = l$stat, df = l$df, p = l$p, max_abs_night_mean = max(abs(nm$mean)),
                            night_sd_range = c(min(nm$std), max(nm$std)))
  out <- list()
  for (reg in REGIONS) {
    w <- d %>% filter(region == reg) %>% select(sensor_id, night_id, dT_night) %>%
      pivot_wider(names_from = night_id, values_from = dT_night) %>% select(-sensor_id)
    cc <- cor(w, use = "pairwise.complete.obs"); v <- cc[upper.tri(cc)]
    out[[reg]] <- list(mean_r = mean(v, na.rm = TRUE), min_r = min(v, na.rm = TRUE), max_r = max(v, na.rm = TRUE))
  }
  RES$night_effect$between_night_site_correlation <<- out
}

# ---- step 10: threshold sensitivity ------------------------------------------------------------------
step10_thresholds <- function(sds) {
  b <- load_binned(); hh <- b$hh; era <- b$era
  base_inv <- night_inventory(hh, era)
  pr <- read_csv(file.path(DER, "site_predictors_multiscale.csv"), show_col_types = FALSE)
  p100 <- pr %>% filter(radius_m == RES$adopted_radius) %>% select(-radius_m, -region)
  base_sn <- sensor_nights(hh, base_inv)
  base_site <- base_sn %>% group_by(sensor_id) %>% summarise(dT = mean(dT_night))
  prep_d <- function(sn) {
    d <- sn %>% inner_join(p100, by = "sensor_id") %>% mutate(region = factor(region, levels = rev(REGIONS)))
    for (p in PRED) d[[paste0(p, "_z")]] <- (d[[p]] - RES$z_means[[p]]) / sds[[p]]
    d
  }
  one <- function(d, inv, w, cc, f) {
    m <- fit_lmm(d, "region + imperv_z + height_z", reml = TRUE); mf <- fit_lmm(d, rhs_full("region"), reml = TRUE)
    site <- d %>% group_by(sensor_id) %>% summarise(dT = mean(dT_night)) %>% inner_join(base_site, by = "sensor_id", suffix = c("", "_base"))
    nn <- inv %>% filter(selected) %>% count(region)
    fe <- fixef(m); se <- sqrt(diag(vcov(m)))
    tibble(wind_max = w, cloud_max = cc, min_frac = f,
           nights_honolulu = sum(nn$n[nn$region == REGIONS[1]]), nights_ewa = sum(nn$n[nn$region == REGIONS[2]]), sensor_nights = nrow(d),
           beta_imperv_z = unname(fe["imperv_z"]), se_imperv = unname(se[names(fe) == "imperv_z"]),
           beta_height_z = unname(fe["height_z"]), se_height = unname(se[names(fe) == "height_z"]),
           r2m_parsimonious = r2_nakagawa(m)$r2m, r2m_full = r2_nakagawa(mf)$r2m,
           site_sd = sd(site$dT), corr_with_baseline = cor(site$dT, site$dT_base),
           mean_wind = mean(d$wind_night), mean_cloud = mean(d$cloud_night))
  }
  grid <- expand.grid(w = c(8, 10, 12, 15, 20), cc = c(15, 25, 35, 50), f = 0.75)
  grid <- rbind(grid, data.frame(w = c(10, 10, 10, 999), cc = c(25, 25, 25, 999), f = c(0.6, 0.9, 1.0, 0)))   # last = no weather filter
  rows <- list()
  for (i in seq_len(nrow(grid))) {
    w <- grid$w[i]; cc <- grid$cc[i]; f <- grid$f[i]
    inv <- night_inventory(hh, era, wind_max = w, cloud_max = cc, min_frac = f)
    if (w == 999) inv <- inv %>% mutate(meets_weather = coalesce(n_era5_bins == NIGHT_BINS, FALSE), selected = meets_weather & meets_network)
    sn <- sensor_nights(hh, inv)
    if (!nrow(sn) || n_distinct(sn$region) < 2) next
    rows[[length(rows) + 1]] <- one(prep_d(sn), inv, w, cc, f)
  }
  # complement: the nights that fail the weather rule (windy and/or cloudy), >= 20 loggers
  inv <- base_inv %>% mutate(selected = !meets_weather & meets_network & coalesce(n_era5_bins == NIGHT_BINS, FALSE))
  rows[[length(rows) + 1]] <- one(prep_d(sensor_nights(hh, inv)), inv, -1, -1, -1)
  t <- bind_rows(rows); write_table(t, "threshold_sensitivity")
  RES$threshold_sensitivity <<- pred_rows(t)
  # night-level relation between the spread of dT across the network and the ERA5 weather
  inv_all <- base_inv %>% filter(coalesce(n_era5_bins == NIGHT_BINS, FALSE), meets_network) %>% mutate(selected = TRUE)
  sn_all <- sensor_nights(hh, inv_all)
  spread <- sn_all %>% group_by(region, night_date) %>%
    summarise(sd_dT = sd(dT_night), wind = mean(wind_night), cloud = mean(cloud_night), n = n(), .groups = "drop") %>%
    inner_join(base_inv %>% select(region, night_date, selected, frac_calm_clear), by = c("region", "night_date"))
  write_table(spread, "night_spread_vs_weather")
  RES$night_spread <<- setNames(lapply(REGIONS, function(reg) {
    g <- spread %>% filter(region == reg)
    cw <- suppressWarnings(cor.test(g$wind, g$sd_dT, method = "spearman", exact = FALSE))
    cc <- suppressWarnings(cor.test(g$cloud, g$sd_dT, method = "spearman", exact = FALSE))
    list(rho_wind = unname(cw$estimate), p_wind = cw$p.value, rho_cloud = unname(cc$estimate), p_cloud = cc$p.value,
         n_nights = nrow(g), sd_selected = mean(g$sd_dT[g$selected]), sd_other = mean(g$sd_dT[!g$selected]))
  }), REGIONS)
  RES$night_inventory <<- setNames(lapply(REGIONS, function(reg) {
    g <- base_inv %>% filter(region == reg); comp <- g %>% filter(coalesce(n_era5_bins == NIGHT_BINS, FALSE))
    list(n_complete = sum(coalesce(g$n_era5_bins == NIGHT_BINS, FALSE) & g$meets_network), n_selected = sum(g$selected),
         n_meets_weather = sum(g$meets_weather),
         excluded_network = I(g$night_date[g$meets_weather & !g$meets_network]),
         excluded_partial = I(g$night_date[coalesce(g$n_era5_bins < NIGHT_BINS, FALSE) & coalesce(g$frac_calm_clear >= MIN_FRAC_STEPS, FALSE)]),
         frac_calm_only = mean(comp$frac_calm >= 0.75), frac_clear_only = mean(comp$frac_clear >= 0.75))
  }), REGIONS)
}

# ---- step 11: model comparison -----------------------------------------------------------------------
step11_comparison <- function(d, m_int_ml) {
  ols <- lm(reformulate(c("region", paste0(PRED, "_z")), "dT_night"), data = d)
  rows <- list(tibble(model = "OLS: region + 8 predictors (no random effect)", k = length(coef(ols)) + 1, aic = AIC(ols), llf = as.numeric(logLik(ols))))
  specs <- list(c("Null: random intercept only", "1"), c("Random intercept + region", "region"),
                c("Random intercept + 8 predictors", rhs_full(NULL)),
                c("Random intercept + region + impervious + height", "region + imperv_z + height_z"),
                c("Random intercept + region + 8 predictors", rhs_full("region")))
  for (s in specs) { m <- fit_lmm(d, s[2], reml = FALSE); rows[[length(rows) + 1]] <- tibble(model = s[1], k = n_params(m), aic = aic_ml(m), llf = llf(m)) }
  rows[[length(rows) + 1]] <- tibble(model = "Random intercept + region x 8 predictors (interactions)", k = n_params(m_int_ml), aic = aic_ml(m_int_ml), llf = llf(m_int_ml))
  t <- bind_rows(rows) %>% mutate(delta_aic = aic - min(aic)); write_table(t, "model_comparison")
  RES$model_comparison <<- pred_rows(t)
}

# ---- step 12: BLUPs, cross-validation, calm core ------------------------------------------------------
step12_blups_cv <- function(d, m_full, sds) {
  re <- ranef(m_full, condVar = TRUE)$sensor_id
  pv <- attr(re, "postVar")[1, 1, ]
  fe <- fixef(m_full)
  sites <- d %>% distinct(sensor_id, .keep_all = TRUE)
  b <- tibble(sensor_id = rownames(re), blup = re[, 1], blup_sd = sqrt(pv)) %>%
    inner_join(sites %>% select(sensor_id, region, latitude, longitude, all_of(paste0(PRED, "_z"))), by = "sensor_id") %>%
    mutate(fixed_pred = unname(fe["(Intercept)"]) + ifelse(region == REGIONS[1], unname(fe[REGION_TERM]), 0) +
             as.numeric(as.matrix(across(all_of(paste0(PRED, "_z")))) %*% fe[paste0(PRED, "_z")])) %>%
    inner_join(d %>% group_by(sensor_id) %>% summarise(dT_site = mean(dT_night)), by = "sensor_id") %>%
    mutate(region = as.character(region)) %>%
    select(sensor_id, region, blup, fixed_pred, dT_site, latitude, longitude, blup_sd) %>% arrange(blup)
  write_table(b, "blups")
  RES$blups <<- list(range = c(min(b$blup), max(b$blup)), sd = sd(b$blup),
                     coolest = pred_rows(head(b, 3) %>% select(sensor_id, region, blup)),
                     warmest = pred_rows(tail(b, 3) %>% select(sensor_id, region, blup)),
                     by_region = setNames(lapply(REGIONS, function(reg) list(mean = mean(b$blup[b$region == reg]), sd = sd(b$blup[b$region == reg]))), REGIONS))
  # leave-one-site-out cross-validation of the site-mean dT (fixed effects only for the held-out site)
  site <- d %>% group_by(sensor_id) %>% summarise(dT_site = mean(dT_night), region = first(region)) %>% arrange(sensor_id)
  predict_loso <- function(dd, rhs, ids) {
    vapply(ids, function(sid) {
      m <- fit_lmm(dd %>% filter(sensor_id != sid), rhs, reml = TRUE)
      te <- dd %>% filter(sensor_id == sid) %>% slice(1)
      unname(predict(m, newdata = te, re.form = NA))
    }, numeric(1))
  }
  cv_row <- function(model, pred, obs, reg_vec, by_region = TRUE) {
    err <- pred - obs
    r <- tibble(model = model, rmse = sqrt(mean(err^2)), mae = mean(abs(err)), r = cor(pred, obs),
                r2_cv = 1 - sum(err^2) / sum((obs - mean(obs))^2))
    for (reg in REGIONS) r[[paste0("rmse_", reg)]] <- if (by_region) sqrt(mean(err[reg_vec == reg]^2)) else NA_real_
    r
  }
  cv <- list()
  for (s in list(c(rhs_full("region"), "full"), c("region + imperv_z + height_z", "parsimonious"), c("region", "region only"))) {
    site[[paste0("pred_", s[2])]] <- predict_loso(d, s[1], site$sensor_id)
    cv[[length(cv) + 1]] <- cv_row(s[2], site[[paste0("pred_", s[2])]], site$dT_site, site$region)
  }
  site$pred_regional_full <- NA_real_; site$pred_regional_best <- NA_real_
  for (reg in REGIONS) {
    g <- d %>% filter(region == reg); sr <- site %>% filter(region == reg)
    bs <- RES$regional[[reg]]$best_subset$predictors
    f_bs <- paste(paste0(strsplit(bs, " \\+ ")[[1]], "_z"), collapse = " + ")
    for (s in list(c(rhs_full(NULL), paste(reg, "full"), "pred_regional_full"), c(f_bs, sprintf("%s best subset (%s)", reg, bs), "pred_regional_best"))) {
      pr <- predict_loso(g, s[1], sr$sensor_id)
      site[[s[3]]][site$region == reg] <- pr
      r <- cv_row(s[2], pr, sr$dT_site, sr$region, by_region = FALSE); r[[paste0("rmse_", reg)]] <- r$rmse
      cv[[length(cv) + 1]] <- r
    }
  }
  write_table(site %>% mutate(region = as.character(region)), "loso_predictions")
  cvt <- bind_rows(cv); write_table(cvt, "loso_cv")
  RES$loso_cv <<- pred_rows(cvt)
  # night-to-night correlation matrices of the site pattern
  for (reg in REGIONS) {
    w <- d %>% filter(region == reg) %>% select(sensor_id, night_id, dT_night) %>%
      pivot_wider(names_from = night_id, values_from = dT_night) %>% select(-sensor_id)
    cc <- round(cor(w, use = "pairwise.complete.obs"), 3)
    write.csv(cbind(night_id = rownames(cc), as.data.frame(cc)), file.path(TAB, paste0("night_correlation_", reg, ".csv")), row.names = FALSE)
  }
  # calm-core check: district best-subset models refitted on nights with mean ERA5 wind < 5 km/h
  RES$calm_core <<- setNames(lapply(REGIONS, function(reg) {
    g <- d %>% filter(region == reg)
    keep <- g %>% group_by(night_id) %>% summarise(w = mean(wind_night)) %>% filter(w < 5) %>% pull(night_id)
    gg <- g %>% filter(night_id %in% keep)
    bs <- RES$regional[[reg]]$best_subset$predictors
    m <- fit_lmm(gg, paste(paste0(strsplit(bs, " \\+ ")[[1]], "_z"), collapse = " + "), reml = TRUE)
    list(nights_kept = I(sort(keep)), n_obs = nrow(gg), coefs = pred_rows(coef_table(m, sds)),
         r2m = r2_nakagawa(m)$r2m, icc_cond = r2_nakagawa(m)$icc_cond)
  }), REGIONS)
}

# ---- step 13: shared nights (absolute temperatures) ----------------------------------------------------
step13_shared_nights <- function() {
  hh <- read_csv(file.path(PROC, "halfhourly_calm_clear.csv"), col_types = cols(datetime = col_character(), time_bin = col_character(), night_date = col_character(), .default = col_guess()))
  shared <- sort(intersect(unique(hh$night_date[hh$region == REGIONS[1]]), unique(hh$night_date[hh$region == REGIONS[2]])))
  t <- bind_rows(lapply(shared, function(nd) bind_rows(lapply(REGIONS, function(reg) {
    g <- hh %>% filter(region == reg, night_date == nd, is_night)
    med <- g %>% group_by(time_bin) %>% summarise(m = first(network_median_temp))
    sm <- g %>% group_by(sensor_id) %>% summarise(T = mean(temp_c))
    tibble(night_date = nd, region = reg, n_sensors = n_distinct(g$sensor_id), T_median_night = mean(med$m), T_mean_all = mean(g$temp_c),
           T_p10 = unname(quantile(sm$T, 0.1)), T_p90 = unname(quantile(sm$T, 0.9)), T_min_site = min(sm$T), T_max_site = max(sm$T),
           sd_sites = sd(sm$T), wind = mean(g$era5_wind_speed_kmh), cloud = mean(g$era5_cloud_cover_pct))
  }))))
  write_table(t, "shared_nights")
  RES$shared_nights <<- list(nights = I(shared), rows = pred_rows(t),
                             median_gap = sapply(shared, function(nd) t$T_median_night[t$night_date == nd & t$region == REGIONS[1]] -
                                                                       t$T_median_night[t$night_date == nd & t$region == REGIONS[2]]))
  prof <- bind_rows(lapply(shared, function(nd) bind_rows(lapply(REGIONS, function(reg)
    hh %>% filter(region == reg, night_date == nd) %>% group_by(hour_dec) %>%
      summarise(network_median_temp = first(network_median_temp), .groups = "drop") %>% mutate(region = reg, night_date = nd)))))
  write_table(prof, "shared_nights_profiles")
}

# ---- main ------------------------------------------------------------------------------------------------
main <- function() {
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
  pr <- read_csv(file.path(DER, "site_predictors_multiscale.csv"), show_col_types = FALSE)
  RES$design <<- setNames(lapply(REGIONS, function(reg) {
    g <- sn %>% filter(region == reg); sm <- g %>% group_by(sensor_id) %>% summarise(dT = mean(dT_night))
    list(n_sites = n_distinct(g$sensor_id), n_nights = n_distinct(g$night_date), n_sensor_nights = nrow(g),
         nights = sort(unique(g$night_date)), dT_site_range = c(min(sm$dT), max(sm$dT)), dT_site_sd = sd(sm$dT),
         temp_night_mean = mean(g$temp_night), wind_mean = mean(g$wind_night), cloud_mean = mean(g$cloud_night))
  }), REGIONS)
  message("Step 1: descriptives"); step1_descriptives(pr)
  message("Step 2: scale selection"); step2_scale()
  # adopt the radius that is best, or equivalent to the best (delta AIC <= 2, Burnham &
  # Anderson 2002), in all three comparisons: the smallest maximum delta AIC across the
  # pooled, Honolulu and ʻEwa full models
  sc <- bind_rows(RES$scale_selection)
  worst <- sc %>% group_by(radius_m) %>% summarise(d = max(delta_aic)) %>% arrange(d)
  adopted <- worst$radius_m[1]
  RES$adopted_radius <<- adopted
  RES$scale_rule <<- list(rule = "minimax delta AIC over pooled / Honolulu / Ewa full models",
                          max_delta_aic = setNames(as.list(worst$d), worst$radius_m))
  message("   adopted radius: ", adopted)
  d <- load_radius(adopted)$d
  z <- zscore_frame(d, PRED); d <- z$df; sds <- z$sds; means <- z$means
  RES$z_sds <<- as.list(sds); RES$z_means <<- as.list(means)
  write_table(d %>% mutate(region = as.character(region)), "analysis_dataset")
  write_table(site_means(d) %>% mutate(region = as.character(region)), "site_summary")
  message("Step 3: collinearity"); step3_collinearity(d)
  message("Step 4: null models"); step4_null(d)
  message("Step 5: single predictors"); step5_single(d, sds)
  message("Step 6: full pooled"); mm <- step6_full(d, sds)
  message("Step 7: district models"); step7_regional(d, sds)
  message("Step 8: interactions"); m_int_ml <- step8_interactions(d, mm$m_ml)
  message("Step 9: night effect"); step9_night(d, mm$m_ml)
  message("Step 10: thresholds"); step10_thresholds(sds)
  message("Step 11: comparison"); step11_comparison(d, m_int_ml)
  message("Step 12: BLUPs, LOSO, calm core"); step12_blups_cv(d, mm$m, sds)
  message("Step 13: shared nights"); step13_shared_nights()
  RES$software <<- list(R = R.version.string, lme4 = as.character(packageVersion("lme4")),
                        lmerTest = as.character(packageVersion("lmerTest")), performance = as.character(packageVersion("performance")))
  write_json(RES, file.path(OUT, "thesis_results.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  message("done")
}

if (sys.nframe() == 0) main()
