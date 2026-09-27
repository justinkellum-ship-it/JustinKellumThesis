Table: Software used for the analysis, with the role of each package and the functions called. Versions are those of the run that produced this document (`results/sessionInfo.txt`). {#tbl:software}

| Package | Version | Role | Functions used | Reference |
|----------------------------------------|----------------------|----------------------------------------|----------------------------------------|----------------------------------------|
| R | 4.3.3 | the language: arithmetic, median(), cor(), lm(), wilcox.test(), cor.test() |  | [@rcore2024] |
| dplyr, tidyr, tibble, purrr, readr | 1.1.4 / 1.3.1 / 2.1.5 | grouped summaries, joins and reshaping of the record; typed CSV input/output | group_by(), summarise(), inner_join(), pivot_wider(), read_csv() | [@wickham2019] |
| lubridate | 1.9.3 | rounding of the loggers' time stamps to the nearest 30 minutes | round_date(), ymd_hms() | [@grolemund2011] |
| jsonlite | 1.8.8 | the results file read by the tables and figures | write_json(), fromJSON() | [@ooms2014] |
| sf | 1.0.15 | vector geometry: projection, buffers, footprint clipping, sight lines for the sky view factor, shoreline distance | st_transform(), st_buffer(), st_intersects(), st_intersection(), st_distance() | [@pebesma2018] |
| terra | 1.7.65 | the 1 m C-CAP masks as virtual mosaics; pixel extraction inside circles; aggregation for base maps | vrt(), extract(), crop(), aggregate() | [@hijmans2024] |
| lme4 | 1.1.35.1 | the mixed-effects models (random intercept per site; REML and ML fits; BLUPs) | lmer(), fixef(), ranef(), VarCorr(), logLik(), AIC(), predict() | [@bates2015] |
| lmerTest | 3.1.3 | Satterthwaite degrees of freedom and t tests for the coefficients | lmer(), summary() | [@kuznetsova2017] |
| performance | 0.10.8 | check of the marginal and conditional R² and ICC computed from the variance components | r2_nakagawa(), icc() | [@ludecke2021] |
| car | 3.1.2 | variance inflation factors | vif() | [@fox2019] |
| ggplot2, patchwork, scales, cowplot | 3.4.4 / 1.2.0 | all figures; multi-panel layout; color scales | ggplot(), geom_sf(), geom_raster(), coord_sf(), wrap_plots() | [@wickham2016; @pedersen2024] |
| rmarkdown, knitr | 2.25 / 1.45 | the analysis walkthrough document | render() |  |
| Python replication: pandas, NumPy, SciPy | 3.0.2 / 2.4.4 / 1.17.1 | record processing and statistics in the independent replication |  | [@mckinney2010; @harris2020; @virtanen2020] |
| Python replication: statsmodels | 0.15.0 | mixed models (MixedLM) in the replication |  | [@seabold2010] |
| Python replication: rasterio, Shapely, GeoPandas | 1.4.4 / 2.1.2 / 1.1.4 | rasters and vector geometry in the replication |  | [@gillies2013; @gillies2007] |
