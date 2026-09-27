Table: Agreement between the R analysis (lme4 / lmerTest) and the independent Python replication (statsmodels MixedLM) on the same inputs: a selection of the 117 paired quantities in `replication_check.csv` (which also holds 25 maximum-difference checks: the descriptors over the 74 sites and ΔT over the 400 sensor-nights). Apart from the night-effect degrees of freedom (an implementation choice: lme4 drops the night dummy that is collinear with the district term), the differences larger than 0.01 are optimizer tolerance (ΔAIC, joint χ²); the leave-one-site-out folds of the over-parameterized ʻEwa full model are the quantities most sensitive to the optimizer in statsmodels. {#tbl:replication}

| Quantity | R | Python | Difference |
|----------------------------------------|------------:|------------:|------------:|
| Scale rule: max ΔAIC at 100 m | 1.671 | 1.674 | -0.003 |
| Scale rule: max ΔAIC at 200 m | 10.196 | 10.206 | -0.010 |
| Scale rule: max ΔAIC at 50 m | 12.344 | 12.324 | +0.020 |
| Null model ICC (pooled) | 0.626 | 0.626 | +0.000 |
| Null model ICC (Honolulu) | 0.810 | 0.810 | +0.000 |
| Null model ICC (Ewa) | 0.365 | 0.365 | +0.000 |
| Full pooled: imperv_z estimate | 0.403 | 0.403 | -0.000 |
| Full pooled: tree_z estimate | -0.143 | -0.143 | +0.000 |
| Full pooled: bldg_z estimate | -0.254 | -0.255 | +0.000 |
| Full pooled: water_z estimate | 0.111 | 0.111 | +0.000 |
| Full pooled: height_z estimate | 0.277 | 0.277 | -0.000 |
| Full pooled: svf_point_z estimate | -0.031 | -0.031 | -0.000 |
| Full pooled: aspect_z estimate | -0.081 | -0.081 | +0.000 |
| Full pooled: coast_km_z estimate | -0.196 | -0.196 | +0.000 |
| Full pooled: R2m | 0.254 | 0.254 | +0.000 |
| Full pooled: R2c | 0.652 | 0.652 | +0.000 |
| Full pooled: LRT vs null, chi2 | 35.083 | 35.093 | -0.009 |
| Honolulu best subset: R2m (imperv + coast_km) | 0.356 | 0.355 | +0.001 |
| Ewa best subset: R2m (imperv + height + aspect) | 0.312 | 0.311 | +0.001 |
| Interactions: joint LRT chi2 | 30.429 | 30.419 | +0.010 |
| BLUP SD | 0.453 | 0.453 | +0.000 |
| LOSO full: r | 0.391 | 0.391 | -0.000 |
| LOSO parsimonious: r | 0.356 | 0.357 | -0.000 |
| LOSO region only: r | -0.467 | -0.467 | -0.000 |
| LOSO Honolulu full: r | 0.538 | 0.538 | +0.000 |
| LOSO Honolulu best subset (imperv + coast_km): r | 0.578 | 0.578 | -0.000 |
| LOSO Ewa full: r | 0.649 | 0.658 | -0.009 |
| LOSO Ewa best subset (imperv + height + aspect): r | 0.705 | 0.705 | -0.000 |
