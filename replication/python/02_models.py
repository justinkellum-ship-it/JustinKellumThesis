#!/usr/bin/env python3
"""
02_models.py
------------
Mixed-effects analysis of nocturnal canopy-layer temperature departures
(dT = sensor - network median, 18:00-06:00 HST means on calm, clear nights)
against two- and three-dimensional urban-form predictors at 50, 100 and
200 m: the Python counterpart of analysis/02_models.R, used to check its results.

Unit of analysis: sensor-night.  Random intercept per sensor site.  Fixed
effects are tested with likelihood-ratio tests on ML fits; estimates and
variance components are reported from REML fits (Zuur et al. 2009).

Steps
  1  descriptive statistics of the predictors by region and radius
  2  scale selection (pooled and regional LMM AIC / R2m per radius, and a
     forward-stepwise OLS on site means as a secondary check)
  3  correlation and variance-inflation screen at the adopted radius
  4  null models and intraclass correlation
  5  single-predictor models (region-adjusted, and per region)
  6  full eight-predictor model with region (pooled)
  7  regional models, best-subset (<= 3 predictors) AIC ranking and
     Akaike-weight predictor importance
  8  region x predictor interaction tests
  9  night fixed-effect check (the reference removes the night mean)
 10  calm/clear threshold sensitivity (night selection re-run from the raw
     record) and the no-filter comparison
 11  model comparison table (ML AIC)
 12  BLUPs and leave-one-site-out cross-validation
 13  absolute temperature comparison on the two shared nights

Outputs: replication/results/thesis_results.json and replication/results/tables/*.csv
"""
import importlib.util
import itertools
import json
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "replication", "results")
PDATA = os.path.join(OUT, "data")   # the processed tables written by 00_prepare_inputs.py and 01_site_predictors.py
TAB = os.path.join(OUT, "tables")
os.makedirs(TAB, exist_ok=True)

spec = importlib.util.spec_from_file_location("prep", os.path.join(HERE, "00_prepare_inputs.py"))
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)

PRED = ["imperv", "tree", "bldg", "water", "height", "svf_point", "aspect", "coast_km"]
LABEL = {"imperv": "Impervious surface fraction (%)", "tree": "Tree canopy fraction (%)",
         "bldg": "Building footprint fraction (%)", "water": "Water surface fraction (%)",
         "height": "Mean building height (m)", "svf_point": "Sky view factor (building, at sensor)",
         "aspect": "Canyon aspect ratio (H/W)", "coast_km": "Distance to coast (km)"}
UNIT = {"imperv": (10, "per +10 pp"), "tree": (10, "per +10 pp"), "bldg": (10, "per +10 pp"),
        "water": (10, "per +10 pp"), "height": (1, "per +1 m"), "svf_point": (-0.1, "per -0.1"),
        "aspect": (0.1, "per +0.1"), "coast_km": (1, "per +1 km")}
RADII = [50, 100, 200]
REGIONS = ["Honolulu", "Ewa"]
RES = {}


# ----------------------------------------------------------------------
# helpers
def _sane(r):
    return bool(np.all(np.isfinite(r.fe_params.values)) and np.abs(r.fe_params.values).max() < 50)


def fit_lmm(df, formula, reml=True):
    md = smf.mixedlm(formula, df, groups=df["sensor_id"])
    r = None
    try:
        r = md.fit(reml=reml, maxiter=3000)
    except Exception:
        pass
    if r is None or not r.converged or not _sane(r):
        for meth in ("powell", "nm", "cg"):
            try:
                r2 = md.fit(reml=reml, method=meth, maxiter=8000)
                if _sane(r2) and (r is None or not _sane(r) or r2.converged):
                    r = r2
                    if r2.converged:
                        break
            except Exception:
                continue
    return r


def n_params(m):
    return int(len(m.fe_params) + 2)      # + random-intercept variance + residual variance


def aic_ml(m):
    return float(-2 * m.llf + 2 * n_params(m))


def lrt(m_full, m_null):
    stat = 2 * (m_full.llf - m_null.llf)
    df = len(m_full.fe_params) - len(m_null.fe_params)
    return float(stat), int(df), float(stats.chi2.sf(stat, df))


def r2_nakagawa(m):
    X = m.model.exog
    var_f = float(np.var(X @ m.fe_params.values, ddof=0))
    var_s = float(m.cov_re.iloc[0, 0])
    var_e = float(m.scale)
    tot = var_f + var_s + var_e
    return dict(r2m=var_f / tot, r2c=(var_f + var_s) / tot, var_fixed=var_f, var_site=var_s, var_resid=var_e,
                icc_cond=var_s / (var_s + var_e))


def coef_table(m, sds=None):
    fe = m.fe_params
    se = m.bse_fe
    rows = []
    for k in fe.index:
        est, s = float(fe[k]), float(se[k])
        rows.append(dict(term=k, estimate=est, se=s, z=est / s, p=float(2 * stats.norm.sf(abs(est / s))),
                         lo=est - 1.96 * s, hi=est + 1.96 * s))
    t = pd.DataFrame(rows)
    if sds is not None:
        # natural-unit slopes for z-scored predictors
        t["natural_unit"] = ""
        t["estimate_natural"] = np.nan
        t["se_natural"] = np.nan
        for i, k in enumerate(t.term):
            base = k.replace("_z", "")
            if k.endswith("_z") and base in sds:
                u, lab = UNIT[base]
                f = u / sds[base]
                t.loc[i, "natural_unit"] = lab
                t.loc[i, "estimate_natural"] = t.loc[i, "estimate"] * f
                t.loc[i, "se_natural"] = t.loc[i, "se"] * abs(f)
    return t


def zscore_frame(df, cols, ref=None):
    """z-scores; scaling parameters taken from `ref` (defaults to df) so that
    coefficients stay comparable across subsets."""
    ref = df if ref is None else ref
    out = df.copy()
    sds, means = {}, {}
    for c in cols:
        mu, sd = ref[c].mean(), ref[c].std(ddof=1)
        out[c + "_z"] = (out[c] - mu) / sd
        sds[c], means[c] = float(sd), float(mu)
    return out, sds, means


def formula_full(extra="region"):
    terms = ([extra] if extra else []) + [p + "_z" for p in PRED]
    return "dT_night ~ " + " + ".join(terms)


# ----------------------------------------------------------------------
def load(radius):
    sn = pd.read_csv(os.path.join(PDATA, "sensor_nights.csv"))
    pr = pd.read_csv(os.path.join(PDATA, "site_predictors_multiscale.csv"))
    p = pr[pr.radius_m == radius].drop(columns=["radius_m", "region"])
    d = sn.merge(p, on="sensor_id")
    d["region"] = pd.Categorical(d.region, categories=["Ewa", "Honolulu"])   # Ewa = reference
    return d, pr


def site_means(d):
    return (d.groupby(["region", "sensor_id"], observed=True)
             .agg(dT_site=("dT_night", "mean"), dT_sd=("dT_night", "std"), n_nights=("dT_night", "size"),
                  **{p: (p, "first") for p in PRED}).reset_index())


# ----------------------------------------------------------------------
def step1_descriptives(pr):
    rows = []
    for (reg, R), g in pr.groupby(["region", "radius_m"]):
        for p in PRED:
            rows.append(dict(region=reg, radius_m=R, predictor=p, mean=g[p].mean(), sd=g[p].std(),
                             min=g[p].min(), median=g[p].median(), max=g[p].max()))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "predictors_by_region.csv"), index=False)
    RES["predictors_by_region"] = {f"{r}_{R}": {p: dict(mean=float(t[(t.region == r) & (t.radius_m == R) & (t.predictor == p)]["mean"].iloc[0]),
                                                         sd=float(t[(t.region == r) & (t.radius_m == R) & (t.predictor == p)]["sd"].iloc[0]))
                                                for p in PRED} for r in REGIONS for R in RADII}


def step1b_region_difference(pr, radius):
    # district difference tests at the adopted radius (site level): Mann-Whitney
    pa = pr[pr.radius_m == radius]
    RES["region_difference"] = dict(radius_m=int(radius))
    for p in PRED:
        a, b = pa[pa.region == REGIONS[0]][p], pa[pa.region == REGIONS[1]][p]
        u, pv = stats.mannwhitneyu(a, b)
        RES["region_difference"][p] = dict(honolulu_median=float(a.median()), ewa_median=float(b.median()), p=float(pv))


def step2_scale(pr, sn):
    rows = []
    single = []
    for R in RADII:
        d, _ = load(R)
        dz, sds, _ = zscore_frame(d, PRED)
        m = fit_lmm(dz, formula_full("region"), reml=False)
        r2 = r2_nakagawa(m)
        rows.append(dict(scope="pooled", radius_m=R, model="region + 8 predictors", aic=aic_ml(m), llf=float(m.llf),
                         r2m=r2["r2m"], r2c=r2["r2c"], var_site=r2["var_site"]))
        for reg in REGIONS:
            dr, sds_r, _ = zscore_frame(d[d.region == reg].copy(), PRED)
            mr = fit_lmm(dr, formula_full(None), reml=False)
            r2r = r2_nakagawa(mr)
            rows.append(dict(scope=reg, radius_m=R, model="8 predictors", aic=aic_ml(mr), llf=float(mr.llf),
                             r2m=r2r["r2m"], r2c=r2r["r2c"], var_site=r2r["var_site"]))
            # single predictors per region (which scale each predictor works best at)
            null = fit_lmm(dr, "dT_night ~ 1", reml=False)
            for p in PRED:
                ms = fit_lmm(dr, f"dT_night ~ {p}_z", reml=False)
                s, df_, pv = lrt(ms, null)
                single.append(dict(region=reg, radius_m=R, predictor=p, beta_z=float(ms.fe_params[f"{p}_z"]),
                                   se=float(ms.bse_fe[f"{p}_z"]), lrt_p=pv, r2m=r2_nakagawa(ms)["r2m"], aic=aic_ml(ms)))
    t = pd.DataFrame(rows)
    t["delta_aic"] = t.groupby("scope").aic.transform(lambda x: x - x.min())
    t.to_csv(os.path.join(TAB, "scale_selection_lmm.csv"), index=False)
    s = pd.DataFrame(single)
    s.to_csv(os.path.join(TAB, "scale_single_predictor.csv"), index=False)
    RES["scale_selection"] = t.to_dict("records")

    # thesis replication: forward stepwise OLS on site means, per region and radius
    steps = []
    for R in RADII:
        d, _ = load(R)
        sm_ = site_means(d)
        for reg in REGIONS:
            g = sm_[sm_.region == reg].copy()
            g, _, _ = zscore_frame(g, PRED)
            remaining = [p + "_z" for p in PRED]
            chosen = []
            y = g.dT_site
            cur = sm.OLS(y, sm.add_constant(np.ones(len(g)) * 0 + 1, has_constant="add")).fit()
            cur_aic = sm.OLS(y, np.ones((len(g), 1))).fit().aic
            steps.append(dict(region=reg, radius_m=R, step=0, added="(null)", aic=cur_aic, r2=0.0, adj_r2=0.0))
            k = 0
            while remaining:
                best = None
                for p in remaining:
                    X = sm.add_constant(g[chosen + [p]])
                    f = sm.OLS(y, X).fit()
                    if best is None or f.aic < best[1]:
                        best = (p, f.aic, f.rsquared, f.rsquared_adj)
                if best[1] < cur_aic - 1e-9:
                    chosen.append(best[0]); remaining.remove(best[0]); cur_aic = best[1]; k += 1
                    steps.append(dict(region=reg, radius_m=R, step=k, added=best[0].replace("_z", ""), aic=best[1],
                                      r2=best[2], adj_r2=best[3]))
                else:
                    break
    st = pd.DataFrame(steps)
    st.to_csv(os.path.join(TAB, "stepwise_ols_site_means.csv"), index=False)
    RES["stepwise_final"] = (st.sort_values("step").groupby(["region", "radius_m"]).tail(1)
                             [["region", "radius_m", "step", "aic", "r2", "adj_r2"]].to_dict("records"))


def step3_collinearity(d):
    out = {}
    for scope, g in [("pooled", d), ("Honolulu", d[d.region == "Honolulu"]), ("Ewa", d[d.region == "Ewa"])]:
        s = g.drop_duplicates("sensor_id")[PRED]
        corr = s.corr(method="pearson")
        corr.to_csv(os.path.join(TAB, f"correlation_{scope}.csv"))
        X = sm.add_constant((s - s.mean()) / s.std())
        vif = {p: float(variance_inflation_factor(X.values, i + 1)) for i, p in enumerate(PRED)}
        out[scope] = dict(vif=vif, max_abs_corr={p: float(corr[p].drop(p).abs().max()) for p in PRED})
    pd.DataFrame({k: v["vif"] for k, v in out.items()}).to_csv(os.path.join(TAB, "vif.csv"))
    RES["collinearity"] = out


def step4_null(d):
    out = {}
    for scope, g in [("pooled", d), ("Honolulu", d[d.region == "Honolulu"]), ("Ewa", d[d.region == "Ewa"])]:
        m = fit_lmm(g, "dT_night ~ 1", reml=True)
        vs, ve = float(m.cov_re.iloc[0, 0]), float(m.scale)
        out[scope] = dict(intercept=float(m.fe_params["Intercept"]), var_site=vs, var_resid=ve, icc=vs / (vs + ve),
                          n_sites=int(g.sensor_id.nunique()), n_obs=int(len(g)),
                          site_sd=float(np.sqrt(vs)), resid_sd=float(np.sqrt(ve)))
    # pooled null with region
    m = fit_lmm(d, "dT_night ~ region", reml=True)
    vs, ve = float(m.cov_re.iloc[0, 0]), float(m.scale)
    out["pooled_region"] = dict(var_site=vs, var_resid=ve, icc=vs / (vs + ve),
                                region_coef=float(m.fe_params["region[T.Honolulu]"]),
                                region_p=float(m.pvalues["region[T.Honolulu]"]))
    pd.DataFrame(out).T.to_csv(os.path.join(TAB, "null_models.csv"))
    RES["null"] = out


def step5_single(d, sds):
    rows = []
    null_p = fit_lmm(d, "dT_night ~ region", reml=False)
    for p in PRED:
        m_ml = fit_lmm(d, f"dT_night ~ region + {p}_z", reml=False)
        m = fit_lmm(d, f"dT_night ~ region + {p}_z", reml=True)
        s, df_, pv = lrt(m_ml, null_p)
        u, lab = UNIT[p]
        b = float(m.fe_params[f"{p}_z"]); se = float(m.bse_fe[f"{p}_z"])
        rows.append(dict(scope="pooled (region-adjusted)", predictor=p, beta_z=b, se_z=se, lrt_chi2=s, lrt_p=pv,
                         r2m=r2_nakagawa(m)["r2m"], beta_natural=b * u / sds[p], se_natural=se * abs(u) / sds[p],
                         natural_unit=lab, aic=aic_ml(m_ml)))
    for reg in REGIONS:
        g = d[d.region == reg]
        null_r = fit_lmm(g, "dT_night ~ 1", reml=False)
        for p in PRED:
            m_ml = fit_lmm(g, f"dT_night ~ {p}_z", reml=False)
            m = fit_lmm(g, f"dT_night ~ {p}_z", reml=True)
            s, df_, pv = lrt(m_ml, null_r)
            u, lab = UNIT[p]
            b = float(m.fe_params[f"{p}_z"]); se = float(m.bse_fe[f"{p}_z"])
            rows.append(dict(scope=reg, predictor=p, beta_z=b, se_z=se, lrt_chi2=s, lrt_p=pv,
                             r2m=r2_nakagawa(m)["r2m"], beta_natural=b * u / sds[p], se_natural=se * abs(u) / sds[p],
                             natural_unit=lab, aic=aic_ml(m_ml)))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "single_predictor_models.csv"), index=False)
    RES["single_predictor"] = t.to_dict("records")


def step6_full(d, sds):
    f = formula_full("region")
    m = fit_lmm(d, f, reml=True)
    m_ml = fit_lmm(d, f, reml=False)
    null_ml = fit_lmm(d, "dT_night ~ 1", reml=False)
    null_reg_ml = fit_lmm(d, "dT_night ~ region", reml=False)
    t = coef_table(m, sds)
    t.to_csv(os.path.join(TAB, "full_model_pooled.csv"), index=False)
    r2 = r2_nakagawa(m)
    s_all, df_all, p_all = lrt(m_ml, null_ml)
    s_m, df_m, p_m = lrt(m_ml, null_reg_ml)
    # region term when predictors are dropped (raw), and the mean-predictor explanation
    RES["full_pooled"] = dict(coefs=t.to_dict("records"), **r2, aic=aic_ml(m_ml), lrt_vs_null=[s_all, df_all, p_all],
                              lrt_morphology_given_region=[s_m, df_m, p_m], converged=bool(m.converged), n_obs=int(len(d)),
                              n_sites=int(d.sensor_id.nunique()))
    # predicted region gap at each region's own mean morphology and at the pooled mean
    Xmean = {p: d.drop_duplicates("sensor_id").groupby("region", observed=True)[p + "_z"].mean() for p in PRED}
    gap_pooled_mean = float(m.fe_params["region[T.Honolulu]"])
    # dT predicted for a site with Honolulu-average morphology placed in each region
    def pred(region, prof):
        v = float(m.fe_params["Intercept"]) + (float(m.fe_params["region[T.Honolulu]"]) if region == "Honolulu" else 0.0)
        for p in PRED:
            v += float(m.fe_params[p + "_z"]) * float(prof[p])
        return v
    hnl_prof = {p: Xmean[p]["Honolulu"] for p in PRED}
    ewa_prof = {p: Xmean[p]["Ewa"] for p in PRED}
    RES["full_pooled"]["region_gap_at_pooled_mean"] = gap_pooled_mean
    RES["full_pooled"]["pred_hnl_profile_in_hnl"] = pred("Honolulu", hnl_prof)
    RES["full_pooled"]["pred_ewa_profile_in_ewa"] = pred("Ewa", ewa_prof)
    RES["full_pooled"]["pred_ewa_profile_in_hnl"] = pred("Honolulu", ewa_prof)
    RES["full_pooled"]["morphology_contribution_hnl_minus_ewa"] = float(sum(
        float(m.fe_params[p + "_z"]) * (hnl_prof[p] - ewa_prof[p]) for p in PRED))
    # parsimonious pooled model (impervious + height + region) for the sensitivity analyses
    mp = fit_lmm(d, "dT_night ~ region + imperv_z + height_z", reml=True)
    mp_ml = fit_lmm(d, "dT_night ~ region + imperv_z + height_z", reml=False)
    tp = coef_table(mp, sds)
    tp.to_csv(os.path.join(TAB, "parsimonious_model_pooled.csv"), index=False)
    RES["parsimonious_pooled"] = dict(coefs=tp.to_dict("records"), **r2_nakagawa(mp), aic=aic_ml(mp_ml),
                                      lrt_full_vs_parsimonious=list(lrt(m_ml, mp_ml)))
    return m, m_ml


def step7_regional(d, sds):
    out = {}
    subsets_rows = []
    for reg in REGIONS:
        g = d[d.region == reg].copy()
        f = formula_full(None)
        m = fit_lmm(g, f, reml=True)
        m_ml = fit_lmm(g, f, reml=False)
        null_ml = fit_lmm(g, "dT_night ~ 1", reml=False)
        t = coef_table(m, sds)
        t["region"] = reg
        t.to_csv(os.path.join(TAB, f"full_model_{reg}.csv"), index=False)
        s, df_, pv = lrt(m_ml, null_ml)
        out[reg] = dict(full=dict(coefs=t.to_dict("records"), **r2_nakagawa(m), aic=aic_ml(m_ml),
                                  lrt_vs_null=[s, df_, pv], converged=bool(m.converged), n_sites=int(g.sensor_id.nunique()),
                                  n_obs=int(len(g))))
        # best subsets (<= 3 predictors), ML AIC, Akaike weights
        cands = [()] + [c for k in (1, 2, 3) for c in itertools.combinations(PRED, k)]
        rows = []
        for c in cands:
            ff = "dT_night ~ 1" if not c else "dT_night ~ " + " + ".join(p + "_z" for p in c)
            mm = fit_lmm(g, ff, reml=False)
            rows.append(dict(region=reg, predictors=" + ".join(c) if c else "(null)", k=len(c), aic=aic_ml(mm),
                             r2m=r2_nakagawa(mm)["r2m"], llf=float(mm.llf)))
        rr = pd.DataFrame(rows)
        rr["delta_aic"] = rr.aic - rr.aic.min()
        rr["weight"] = np.exp(-0.5 * rr.delta_aic) / np.exp(-0.5 * rr.delta_aic).sum()
        rr = rr.sort_values("aic").reset_index(drop=True)
        subsets_rows.append(rr)
        importance = {p: float(rr[rr.predictors.str.contains(p)].weight.sum()) for p in PRED}
        best = rr.iloc[0]
        fb = "dT_night ~ " + " + ".join(p + "_z" for p in best.predictors.split(" + ")) if best.k > 0 else "dT_night ~ 1"
        mb = fit_lmm(g, fb, reml=True)
        tb = coef_table(mb, sds)
        tb["region"] = reg
        tb.to_csv(os.path.join(TAB, f"best_subset_model_{reg}.csv"), index=False)
        out[reg]["best_subset"] = dict(predictors=best.predictors, aic=float(best.aic), r2m=float(best.r2m),
                                       coefs=tb.to_dict("records"), **{k: v for k, v in r2_nakagawa(mb).items() if k != "r2m"},
                                       importance=importance, top5=rr.head(5).to_dict("records"),
                                       n_within_2=int((rr.delta_aic <= 2).sum()))
    pd.concat(subsets_rows).to_csv(os.path.join(TAB, "best_subsets.csv"), index=False)
    RES["regional"] = out


def step8_interactions(d, m_ml_base):
    rows = []
    base_f = formula_full("region")
    for p in PRED:
        f = base_f + f" + region:{p}_z"
        m = fit_lmm(d, f, reml=False)
        s, df_, pv = lrt(m, m_ml_base)
        mr = fit_lmm(d, f, reml=True)
        key = f"region[T.Honolulu]:{p}_z"
        rows.append(dict(predictor=p, lrt_chi2=s, df=df_, p=pv, interaction_coef=float(mr.fe_params[key]),
                         se=float(mr.bse_fe[key]), ewa_slope_z=float(mr.fe_params[p + "_z"]),
                         hnl_slope_z=float(mr.fe_params[p + "_z"] + mr.fe_params[key])))
    f_all = base_f + " + " + " + ".join(f"region:{p}_z" for p in PRED)
    m_all = fit_lmm(d, f_all, reml=False)
    s, df_, pv = lrt(m_all, m_ml_base)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "interaction_tests.csv"), index=False)
    RES["interactions"] = dict(single=t.to_dict("records"), joint=dict(lrt_chi2=s, df=df_, p=pv, aic=aic_ml(m_all),
                                                                        r2m=r2_nakagawa(m_all)["r2m"]))
    return m_all


def step9_night(d, m_ml_base):
    f = formula_full("region") + " + C(night_id)"
    m = fit_lmm(d, f, reml=False)
    s, df_, pv = lrt(m, m_ml_base)
    night_means = d.groupby(["region", "night_id"], observed=True).dT_night.agg(["mean", "std", "size"]).reset_index()
    night_means.to_csv(os.path.join(TAB, "night_means.csv"), index=False)
    RES["night_effect"] = dict(lrt_chi2=s, df=df_, p=pv, max_abs_night_mean=float(night_means["mean"].abs().max()),
                               night_sd_range=[float(night_means["std"].min()), float(night_means["std"].max())])
    # night-to-night stability of the spatial pattern: correlation of site dT between nights
    out = {}
    for reg in REGIONS:
        w = d[d.region == reg].pivot(index="sensor_id", columns="night_id", values="dT_night")
        c = w.corr()
        vals = c.values[np.triu_indices_from(c.values, 1)]
        out[reg] = dict(mean_r=float(np.nanmean(vals)), min_r=float(np.nanmin(vals)), max_r=float(np.nanmax(vals)))
    RES["night_effect"]["between_night_site_correlation"] = out


def step10_thresholds(sds):
    hh, wx = prep.load_binned()
    base_inv = prep.night_inventory(hh, wx)
    pr = pd.read_csv(os.path.join(PDATA, "site_predictors_multiscale.csv"))
    p100 = pr[pr.radius_m == RES["adopted_radius"]].drop(columns=["radius_m", "region"])
    coords = None
    base_sn = prep.sensor_nights(hh, base_inv)
    base_site = base_sn.groupby("sensor_id").dT_night.mean()
    rows = []
    grid = [(w, c, f) for w in (8, 10, 12, 15, 20) for c in (15, 25, 35, 50) for f in (0.75,)]
    grid += [(10, 25, 0.6), (10, 25, 0.9), (10, 25, 1.0), (999, 999, 0.0)]      # last = no weather filter
    for w, c, f in grid:
        inv = prep.night_inventory(hh, wx, wind_max=w, cloud_max=c, min_frac=f)
        if w == 999:
            inv["meets_weather"] = inv.n_weather_bins == 24
            inv["selected"] = inv.meets_weather & inv.meets_network
        sn = prep.sensor_nights(hh, inv)
        if sn.empty or sn.region.nunique() < 2:
            continue
        d = sn.merge(p100, on="sensor_id")
        d["region"] = pd.Categorical(d.region, categories=["Ewa", "Honolulu"])
        for p in PRED:
            d[p + "_z"] = (d[p] - RES["z_means"][p]) / sds[p]
        m = fit_lmm(d, "dT_night ~ region + imperv_z + height_z", reml=True)
        mf = fit_lmm(d, formula_full("region"), reml=True)
        site = d.groupby("sensor_id").dT_night.mean()
        common = site.index.intersection(base_site.index)
        r = float(np.corrcoef(site[common], base_site[common])[0, 1])
        n_nights = inv[inv.selected].groupby("region").size()
        rows.append(dict(wind_max=w, cloud_max=c, min_frac=f, nights_honolulu=int(n_nights.get("Honolulu", 0)),
                         nights_ewa=int(n_nights.get("Ewa", 0)), sensor_nights=int(len(d)),
                         beta_imperv_z=float(m.fe_params["imperv_z"]), se_imperv=float(m.bse_fe["imperv_z"]),
                         beta_height_z=float(m.fe_params["height_z"]), se_height=float(m.bse_fe["height_z"]),
                         r2m_parsimonious=r2_nakagawa(m)["r2m"], r2m_full=r2_nakagawa(mf)["r2m"],
                         site_sd=float(site.std()), corr_with_baseline=r,
                         mean_wind=float(d.wind_night.mean()), mean_cloud=float(d.cloud_night.mean())))
    # complement: nights that fail the baseline weather rule (windy/cloudy), >=20 sensors
    inv = base_inv.copy()
    inv["selected"] = (~inv.meets_weather) & inv.meets_network & (inv.n_weather_bins == 24)
    sn = prep.sensor_nights(hh, inv)
    d = sn.merge(p100, on="sensor_id")
    d["region"] = pd.Categorical(d.region, categories=["Ewa", "Honolulu"])
    for p in PRED:
        d[p + "_z"] = (d[p] - RES["z_means"][p]) / sds[p]
    m = fit_lmm(d, "dT_night ~ region + imperv_z + height_z", reml=True)
    mf = fit_lmm(d, formula_full("region"), reml=True)
    site = d.groupby("sensor_id").dT_night.mean()
    common = site.index.intersection(base_site.index)
    n_nights = inv[inv.selected].groupby("region").size()
    rows.append(dict(wind_max=-1, cloud_max=-1, min_frac=-1, nights_honolulu=int(n_nights.get("Honolulu", 0)),
                     nights_ewa=int(n_nights.get("Ewa", 0)), sensor_nights=int(len(d)),
                     beta_imperv_z=float(m.fe_params["imperv_z"]), se_imperv=float(m.bse_fe["imperv_z"]),
                     beta_height_z=float(m.fe_params["height_z"]), se_height=float(m.bse_fe["height_z"]),
                     r2m_parsimonious=r2_nakagawa(m)["r2m"], r2m_full=r2_nakagawa(mf)["r2m"],
                     site_sd=float(site.std()), corr_with_baseline=float(np.corrcoef(site[common], base_site[common])[0, 1]),
                     mean_wind=float(d.wind_night.mean()), mean_cloud=float(d.cloud_night.mean())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "threshold_sensitivity.csv"), index=False)
    RES["threshold_sensitivity"] = t.to_dict("records")
    # night-level relation: spread of dT across the network vs wind, all complete nights with >=20 sensors
    inv_all = base_inv[(base_inv.n_weather_bins == 24) & base_inv.meets_network].copy()
    inv_all["selected"] = True
    sn_all = prep.sensor_nights(hh, inv_all)
    spread = (sn_all.groupby(["region", "night_date"]).agg(sd_dT=("dT_night", "std"), wind=("wind_night", "mean"),
                                                            cloud=("cloud_night", "mean"), n=("dT_night", "size"))
              .reset_index())
    spread = spread.merge(base_inv[["region", "night_date", "selected", "frac_calm_clear"]], on=["region", "night_date"])
    spread.to_csv(os.path.join(TAB, "night_spread_vs_weather.csv"), index=False)
    RES["night_spread"] = {reg: dict(rho_wind=float(stats.spearmanr(g.wind, g.sd_dT)[0]), p_wind=float(stats.spearmanr(g.wind, g.sd_dT)[1]),
                                    rho_cloud=float(stats.spearmanr(g.cloud, g.sd_dT)[0]), p_cloud=float(stats.spearmanr(g.cloud, g.sd_dT)[1]),
                                    n_nights=int(len(g)), sd_selected=float(g[g.selected].sd_dT.mean()),
                                    sd_other=float(g[~g.selected].sd_dT.mean()))
                           for reg, g in spread.groupby("region")}
    # inventory summary for the text
    RES["night_inventory"] = {reg: dict(n_complete=int(((g.n_weather_bins == 24) & g.meets_network).sum()),
                                        n_selected=int(g.selected.sum()),
                                        n_meets_weather=int(g.meets_weather.sum()),
                                        excluded_network=g[g.meets_weather & ~g.meets_network].night_date.tolist(),
                                        excluded_partial=g[(g.n_weather_bins < 24) & (g.frac_calm_clear >= 0.75)].night_date.tolist(),
                                        frac_calm_only=float((g[(g.n_weather_bins == 24)].frac_calm >= 0.75).mean()),
                                        frac_clear_only=float((g[(g.n_weather_bins == 24)].frac_clear >= 0.75).mean()))
                              for reg, g in base_inv.groupby("region")}


def step11_comparison(d, m_full_ml, m_int_ml):
    rows = []
    X = sm.add_constant(pd.get_dummies(d[["region"] + [p + "_z" for p in PRED]], drop_first=True, dtype=float))
    ols = sm.OLS(d.dT_night, X).fit()
    rows.append(dict(model="OLS: region + 8 predictors (no random effect)", k=int(X.shape[1] + 1), aic=float(ols.aic), llf=float(ols.llf)))
    for name, f in [("Null: random intercept only", "dT_night ~ 1"),
                    ("Random intercept + region", "dT_night ~ region"),
                    ("Random intercept + 8 predictors", formula_full(None)),
                    ("Random intercept + region + impervious + height", "dT_night ~ region + imperv_z + height_z"),
                    ("Random intercept + region + 8 predictors", formula_full("region"))]:
        m = fit_lmm(d, f, reml=False)
        rows.append(dict(model=name, k=n_params(m), aic=aic_ml(m), llf=float(m.llf)))
    rows.append(dict(model="Random intercept + region x 8 predictors (interactions)", k=n_params(m_int_ml), aic=aic_ml(m_int_ml), llf=float(m_int_ml.llf)))
    t = pd.DataFrame(rows)
    t["delta_aic"] = t.aic - t.aic.min()
    t.to_csv(os.path.join(TAB, "model_comparison.csv"), index=False)
    RES["model_comparison"] = t.to_dict("records")


def step12_blups_cv(d, m_full, sds):
    re = m_full.random_effects
    rows = []
    fe = m_full.fe_params
    for sid, v in re.items():
        g = d[d.sensor_id == sid].iloc[0]
        fixed = float(fe["Intercept"]) + (float(fe["region[T.Honolulu]"]) if g.region == "Honolulu" else 0.0) + \
            sum(float(fe[p + "_z"]) * float(g[p + "_z"]) for p in PRED)
        rows.append(dict(sensor_id=sid, region=str(g.region), blup=float(v.iloc[0]), fixed_pred=fixed,
                         dT_site=float(d[d.sensor_id == sid].dT_night.mean()), latitude=float(g.latitude), longitude=float(g.longitude)))
    b = pd.DataFrame(rows)
    # conditional SD of the BLUPs (approx.: from the model's random_effects_cov)
    try:
        b["blup_sd"] = [float(np.sqrt(m_full.random_effects_cov[s].iloc[0, 0])) for s in b.sensor_id]
    except Exception:
        b["blup_sd"] = np.nan
    b = b.sort_values("blup")
    b.to_csv(os.path.join(TAB, "blups.csv"), index=False)
    RES["blups"] = dict(range=[float(b.blup.min()), float(b.blup.max())], sd=float(b.blup.std()),
                        coolest=b.head(3)[["sensor_id", "region", "blup"]].to_dict("records"),
                        warmest=b.tail(3)[["sensor_id", "region", "blup"]].to_dict("records"),
                        by_region={reg: dict(mean=float(b[b.region == reg].blup.mean()), sd=float(b[b.region == reg].blup.std()))
                                   for reg in REGIONS})
    # leave-one-site-out cross-validation of site-mean dT
    site = d.groupby("sensor_id").agg(dT_site=("dT_night", "mean"), region=("region", "first")).reset_index()
    cv = []
    for f, name in [(formula_full("region"), "full"), ("dT_night ~ region + imperv_z + height_z", "parsimonious"),
                    ("dT_night ~ region", "region only")]:
        preds = []
        for sid in site.sensor_id:
            tr = d[d.sensor_id != sid]
            te = d[d.sensor_id == sid].iloc[[0]]
            m = fit_lmm(tr, f, reml=True)
            Xte = pd.DataFrame({k: [1.0 if k == "Intercept" else
                                    (1.0 if (k == "region[T.Honolulu]" and te.region.iloc[0] == "Honolulu") else
                                     (0.0 if k == "region[T.Honolulu]" else float(te[k].iloc[0])))]
                                for k in m.fe_params.index})
            preds.append(float((Xte.values @ m.fe_params.values)[0]))
        site["pred_" + name] = preds
        err = site["pred_" + name] - site.dT_site
        cv.append(dict(model=name, rmse=float(np.sqrt((err ** 2).mean())), mae=float(err.abs().mean()),
                       r=float(np.corrcoef(site["pred_" + name], site.dT_site)[0, 1]),
                       r2_cv=float(1 - (err ** 2).sum() / ((site.dT_site - site.dT_site.mean()) ** 2).sum()),
                       **{f"rmse_{reg}": float(np.sqrt((err[site.region == reg] ** 2).mean())) for reg in REGIONS}))
    # regional models (full and best subset) cross-validated within their region
    for reg in REGIONS:
        g_all = d[d.region == reg]
        sr = site[site.region == reg].copy()
        bs = RES["regional"][reg]["best_subset"]["predictors"]
        f_bs = "dT_night ~ " + " + ".join(p + "_z" for p in bs.split(" + "))
        for f, name in [(formula_full(None), f"{reg} full"), (f_bs, f"{reg} best subset ({bs})")]:
            preds = []
            for sid in sr.sensor_id:
                tr = g_all[g_all.sensor_id != sid]
                te = g_all[g_all.sensor_id == sid].iloc[[0]]
                m = fit_lmm(tr, f, reml=True)
                x = np.array([1.0 if k == "Intercept" else float(te[k].iloc[0]) for k in m.fe_params.index])
                preds.append(float(x @ m.fe_params.values))
            err = np.array(preds) - sr.dT_site.values
            col = "pred_regional_full" if name.endswith("full") else "pred_regional_best"
            site.loc[site.region == reg, col] = preds
            cv.append(dict(model=name, rmse=float(np.sqrt((err ** 2).mean())), mae=float(np.abs(err).mean()),
                           r=float(np.corrcoef(preds, sr.dT_site)[0, 1]),
                           r2_cv=float(1 - (err ** 2).sum() / ((sr.dT_site - sr.dT_site.mean()) ** 2).sum()),
                           **{f"rmse_{r_}": (float(np.sqrt((err ** 2).mean())) if r_ == reg else np.nan) for r_ in REGIONS}))
    site.to_csv(os.path.join(TAB, "loso_predictions.csv"), index=False)
    pd.DataFrame(cv).to_csv(os.path.join(TAB, "loso_cv.csv"), index=False)
    RES["loso_cv"] = cv
    # night-to-night correlation matrices of the site pattern
    for reg in REGIONS:
        w = d[d.region == reg].pivot(index="sensor_id", columns="night_id", values="dT_night")
        w.corr().round(3).to_csv(os.path.join(TAB, f"night_correlation_{reg}.csv"))
    # calm-core check: regional best-subset models refitted on nights with mean wind < 5 km/h
    core = {}
    for reg in REGIONS:
        g = d[(d.region == reg)]
        nights_ok = g.groupby("night_id").wind_night.mean()
        keep = nights_ok[nights_ok < 5].index
        gg = g[g.night_id.isin(keep)]
        bs = RES["regional"][reg]["best_subset"]["predictors"]
        f_bs = "dT_night ~ " + " + ".join(p + "_z" for p in bs.split(" + "))
        m = fit_lmm(gg, f_bs, reml=True)
        core[reg] = dict(nights_kept=sorted(keep.tolist()), n_obs=int(len(gg)), coefs=coef_table(m, sds).to_dict("records"),
                         r2m=r2_nakagawa(m)["r2m"], icc_cond=r2_nakagawa(m)["icc_cond"])
    RES["calm_core"] = core


def step13_shared_nights():
    hh = pd.read_csv(os.path.join(PDATA, "halfhourly_calm_clear.csv"), parse_dates=["datetime", "time_bin"])
    shared = sorted(set(hh[hh.region == "Honolulu"].night_date) & set(hh[hh.region == "Ewa"].night_date))
    rows = []
    for nd in shared:
        for reg in REGIONS:
            g = hh[(hh.region == reg) & (hh.night_date == nd) & hh.is_night]
            med = g.groupby("time_bin").network_median_temp.first()
            rows.append(dict(night_date=nd, region=reg, n_sensors=int(g.sensor_id.nunique()),
                             T_median_night=float(med.mean()), T_mean_all=float(g.temp_c.mean()),
                             T_p10=float(g.groupby("sensor_id").temp_c.mean().quantile(0.1)),
                             T_p90=float(g.groupby("sensor_id").temp_c.mean().quantile(0.9)),
                             T_min_site=float(g.groupby("sensor_id").temp_c.mean().min()),
                             T_max_site=float(g.groupby("sensor_id").temp_c.mean().max()),
                             sd_sites=float(g.groupby("sensor_id").temp_c.mean().std()),
                             wind=float(g.wind_kmh.mean()), cloud=float(g.cloud_pct.mean())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "shared_nights.csv"), index=False)
    RES["shared_nights"] = dict(nights=shared, rows=t.to_dict("records"),
                                median_gap=[float(t[(t.night_date == nd) & (t.region == "Honolulu")].T_median_night.iloc[0]
                                                  - t[(t.night_date == nd) & (t.region == "Ewa")].T_median_night.iloc[0]) for nd in shared])
    # diurnal median profile difference on shared nights
    prof = []
    for nd in shared:
        for reg in REGIONS:
            g = hh[(hh.region == reg) & (hh.night_date == nd)]
            m = g.groupby("hour_dec").network_median_temp.first().reset_index()
            m["region"] = reg; m["night_date"] = nd
            prof.append(m)
    pd.concat(prof).to_csv(os.path.join(TAB, "shared_nights_profiles.csv"), index=False)


def main():
    sn = pd.read_csv(os.path.join(PDATA, "sensor_nights.csv"))
    pr = pd.read_csv(os.path.join(PDATA, "site_predictors_multiscale.csv"))
    RES["design"] = {reg: dict(n_sites=int(sn[sn.region == reg].sensor_id.nunique()),
                               n_nights=int(sn[sn.region == reg].night_date.nunique()),
                               n_sensor_nights=int((sn.region == reg).sum()),
                               nights=sorted(sn[sn.region == reg].night_date.unique().tolist()),
                               dT_site_range=[float(sn[sn.region == reg].groupby("sensor_id").dT_night.mean().min()),
                                              float(sn[sn.region == reg].groupby("sensor_id").dT_night.mean().max())],
                               dT_site_sd=float(sn[sn.region == reg].groupby("sensor_id").dT_night.mean().std()),
                               temp_night_mean=float(sn[sn.region == reg].temp_night.mean()),
                               wind_mean=float(sn[sn.region == reg].wind_night.mean()),
                               cloud_mean=float(sn[sn.region == reg].cloud_night.mean()))
                     for reg in REGIONS}
    print("Step 1: descriptives"); step1_descriptives(pr)
    print("Step 2: scale selection"); step2_scale(pr, sn)
    # adopt the radius that is best, or equivalent to the best (delta AIC <= 2,
    # Burnham & Anderson 2002), in all three comparisons: the smallest maximum
    # delta AIC across the pooled, Honolulu and Ewa models
    sc = pd.DataFrame(RES["scale_selection"])
    worst = sc.groupby("radius_m").delta_aic.max().sort_values()
    adopted = int(worst.index[0])
    RES["adopted_radius"] = adopted
    RES["scale_rule"] = dict(rule="minimax delta AIC over pooled / Honolulu / Ewa full models",
                             max_delta_aic={int(k): float(v) for k, v in worst.items()})
    print("   adopted radius:", adopted, worst.to_dict())
    step1b_region_difference(pr, adopted)
    d, _ = load(adopted)
    d, sds, means = zscore_frame(d, PRED)
    RES["z_sds"], RES["z_means"] = sds, means
    d.to_csv(os.path.join(TAB, "analysis_dataset.csv"), index=False)
    site_means(d).to_csv(os.path.join(TAB, "site_summary.csv"), index=False)
    print("Step 3: collinearity"); step3_collinearity(d)
    print("Step 4: null models"); step4_null(d)
    print("Step 5: single predictors"); step5_single(d, sds)
    print("Step 6: full pooled"); m_full, m_full_ml = step6_full(d, sds)
    print("Step 7: regional"); step7_regional(d, sds)
    print("Step 8: interactions"); m_int_ml = step8_interactions(d, m_full_ml)
    print("Step 9: night effect"); step9_night(d, m_full_ml)
    print("Step 10: thresholds"); step10_thresholds(sds)
    print("Step 11: comparison"); step11_comparison(d, m_full_ml, m_int_ml)
    print("Step 12: BLUPs and LOSO CV"); step12_blups_cv(d, m_full, sds)
    print("Step 13: shared nights"); step13_shared_nights()
    with open(os.path.join(OUT, "thesis_results.json"), "w") as f:
        json.dump(RES, f, indent=1, default=float)
    print("done")


if __name__ == "__main__":
    main()
