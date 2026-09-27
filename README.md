# Nocturnal heat in leeward coastal Oʻahu: sensor network, data and analysis

Data, code and results of the M.A. thesis *Island Urban Climatology: Measuring and Modeling
Nocturnal Heat in Leeward Coastal Oʻahu, Hawaiʻi* (Justin Kellum, Department of Geography and
Environment, University of Hawaiʻi at Mānoa, 2027).

One hundred low-cost iButton stations recorded canopy-layer air temperature every 30 minutes in
two districts of leeward Oʻahu between November 2024 and January 2025: the compact, partly
high-rise core of Honolulu (36 recovered sites) and the low-rise suburban plain of ʻEwa–Kapolei
(38 sites). On the 11 calm, clear nights identified from the ECMWF analysis winds and cloud cover
(400 sensor-nights), each site's
night-mean departure from its district network median (ΔT) is related to eight descriptors of
urban form and setting, computed at 50, 100 and 200 m from 1 m NOAA C-CAP land cover, FEMA
building footprints and the GSHHG shoreline, with linear mixed-effects models (random intercept
per site). Every number, table and figure in the thesis is produced from the files in this
repository by a single command.

**Project page:** https://justinkellum-ship-it.github.io/JustinKellumThesis/ (thesis, figures,
and the step-by-step walkthrough of the analysis)

**Using the pipeline on your own logger network:** `AGENTS.md` is written for an AI coding
agent (Claude Code, Codex, Gemini and the like) and its operator: the questions to settle first
(logging interval, time zone, districts, data sources), the configuration block in
`analysis/helpers.R` that they map to, the run order and the checks at each step.

## Layout

```
data/
  raw/           the study's own measurements
                   exports/              the loggers' export tables as downloaded (one per district: every reading,
                                         plus one metadata row per logger from the download application)
                   logger_readings.csv   the 151,552 readings compiled from them: district, logger, time stamp (HST), °C
                   logger_metadata.csv   one row per logger: model, serial number, sample rate, mission start, sample counts
                   sites.csv             the 74 loggers: district, WGS84 position, setting recorded at deployment
                                         (sensor_id is the label given to the station at deployment; the analysis
                                         uses it only as an identifier)
  external/      data from other providers (small files included; large ones downloaded, see below)
                   weather_hourly_series.csv         10 m wind and cloud cover of the ECMWF analysis at one grid point
                                                     per district, from Open-Meteo: the hourly value at the nearest
                                                     hour, one row per logger time stamp (see Data sources)
                   oahu_gshhs_f.geojson              Oʻahu polygon from the GSHHG full-resolution shoreline
                   ne_10m_land_hawaii.geojson, ne_10m_roads_oahu.geojson   Natural Earth 10 m land and roads, the
                                                     features around the islands (analysis/inputs/extract_basemap_layers.py)
                   oahu_dem_500m_q5.txt, honolulu_dem_150m_q5.txt         USGS 10 m DEM subsamples for the maps of
                                                     Appendix B (the header of each file states the ERDDAP request)
                   noaa_1612340_water_temperature.csv                     Honolulu tide-station water temperature
                                                     (quoted in the text; not used by any script; the header
                                                     states the NOAA API request)
  processed/     written by analysis/00_prepare_inputs.R
                   night_inventory.csv         every district × night with the inclusion decision
                   sensor_nights.csv           the 400 sensor-nights: the unit of analysis
                   halfhourly_calm_clear.csv   the half-hourly record of the analysed diurnal cycles
  descriptors/   written by analysis/01_site_predictors.R and the checks (site_elevation.csv by
                 analysis/inputs/fetch_site_elevation.py)
                   site_predictors_multiscale.csv   the eight descriptors of every site at 50, 100 and 200 m
                   site_elevation.csv, site_default_height_share.csv, island_default_height_share.csv,
                   site_height_fill3m_100m.csv, site_svf_footprint_check.csv
analysis/        the R pipeline, numbered by step; run_all.R runs the scripts in the order of the table below
  inputs/        how the external inputs were made: fetch_weather_series.py (the weather series from Open-Meteo),
                 fetch_site_elevation.py (site elevations from the USGS EPQS), prepare_buildings.R (the Oʻahu
                 footprints from the FEMA geodatabase), extract_island_polygon.R (the island polygon from GSHHG),
                 extract_basemap_layers.py (the Natural Earth layers)
results/         thesis_results.json (every statistic in the thesis), tables/*.csv, figures/ (PNG at 300 dpi, plus
                 PDF for the statistical figures; FIGURES.md maps the files to the thesis figure numbers),
                 thesis_tables/*.md (21 of the thesis's 22 tables as printed; the descriptor-definition table is
                 text), elevation_check.json, sensitivity_checks.json, island_landcover.json (the island and
                 domain land-cover shares quoted in Section 1.2 and Appendix B), sessionInfo.txt (the R and package versions), cache/
                 (4 m land-cover rasters written by 03_figures.R for the maps)
replication/     an independent Python implementation of the same pipeline (python/), its results
                 (results/) and the quantity-by-quantity comparison with R (replication_check.csv)
walkthrough/     Analysis_walkthrough.Rmd and its rendered HTML and PDF: every processing rule, model and
                 check explained and run, with the package that does each step
docs/            the project website
```

## The analysis, script by script

| Script | What it does | Runtime |
|---|---|---|
| `00a_compile_logger_readings.R` | reads the loggers' export tables, separates the metadata rows, writes `logger_readings.csv` and `logger_metadata.csv`, and checks the stored record against the exports reading for reading | seconds |
| `00_prepare_inputs.R` | snaps the loggers' clocks to 30-minute bins, computes the district network median and ΔT, joins the weather series, classifies every night, writes the sensor-night table | seconds |
| `01_site_predictors.R` | the eight descriptors of every site at 50, 100 and 200 m: impervious, tree-canopy, building-footprint and water fractions, mean building height, canyon aspect ratio, ray-cast sky view factor, distance to the coast | ~7 min |
| `02_models.R` | every mixed-effects model: scale selection, collinearity, null models, single descriptors, full pooled and district models, best subsets and Akaike weights, district × descriptor interactions, night effect, threshold sensitivity, model comparison, site effects and leave-one-site-out cross-validation, shared nights | ~1 min |
| `03_figures.R` | the statistical figures and maps of the main text; also `island_landcover.json` (the land-cover shares of the island and the two domains) and the 4 m rasters in `results/cache/` that the maps draw | ~10 min |
| `04_replication_check.R` | compares the R results with the Python implementation | seconds |
| `05_thesis_tables.R` | 21 of the thesis's 22 tables, in Markdown, as printed (the descriptor-definition table is text) | seconds |
| `06_elevation_check.R` | the terrain check: site elevation (USGS 3DEP, `data/descriptors/site_elevation.csv` from `inputs/fetch_site_elevation.py`) added to the models, best subsets and cross-validation | ~1 min |
| `07_sensitivity_checks.R` | three site-level checks: without the warmest site; with 3 m instead of 5 m for the buildings whose height the inventory does not record (and without the sites dominated by them); with the sky view factor of positions inside footprints recomputed | ~7 min |
| `08_appendix_figures.R` | the land-cover and terrain maps and the sky view factor diagram of Appendix B | ~1 min |

`helpers.R` holds the paths, the study configuration (every threshold is a named constant), the
model formulas and the z-scoring; the model-fitting helpers (`fit_lmm()`, `lrt()`, `coef_table()`,
`r2_nakagawa()`) live in `02_models.R`, which the check scripts source; `maplib.R` holds the map
layers. `run_all.R` runs `00a` to `08` and then `05` (whose terrain-check and software tables use
the results of `06` and the versions of the run). The walkthrough shows the code of each step with
its output.

### Processing rules (constants in `analysis/helpers.R`, implemented in `00_prepare_inputs.R`)

1. Readings snapped to the nearest 30-minute bin (the loggers' clocks differ by up to 2 min).
2. Reference = median of all loggers of the district reporting in the bin; ΔT = T − median.
3. Night = 18:00–06:00 HST; a diurnal date starts at 06:00.
4. Calm and clear night: ECMWF analysis 10 m wind < 10 km/h and cloud < 25 % in ≥ 75 % of the 24 night bins; all 24 weather bins present.
5. ≥ 20 loggers reporting; a sensor-night needs ≥ 18 of its 24 bins.

## Reproducing the results

```
Rscript analysis/requirements.R     # installs the packages once (CRAN)
Rscript analysis/run_all.R          # ~30 minutes: data/processed, data/descriptors, results (all scripts, in order)
```

Requires R ≥ 4.3. The versions used for the thesis are in `results/sessionInfo.txt` (R 4.3;
lme4, lmerTest, performance, car, sf, terra, ggplot2, patchwork, dplyr, tidyr, readr, lubridate,
jsonlite). `00_prepare_inputs.R`, `02_models.R` and everything else that reads the stored
descriptors run from the files in the repository; recomputing the descriptors
(`01_site_predictors.R`), the maps and the site checks needs three large inputs placed in
`data/external/`:

| Input | Where | Put it at |
|---|---|---|
| NOAA C-CAP 2021 high-resolution (1 m) land cover for Hawaiʻi: the impervious, canopy and water masks for Oʻahu | https://coast.noaa.gov/digitalcoast/data/ccaphighres.html | `data/external/ccap/*.tif` |
| FEMA / ORNL USA Structures, Hawaiʻi (`HI_Structures.gdb`) | https://gis-fema.hub.arcgis.com/pages/usa-structures | `data/external/buildings_oahu.gpkg`, made from the geodatabase by `analysis/inputs/prepare_buildings.R` (the 207,031 Oʻahu footprints with `height_m`, `height_recorded` and `area_m2`; the 46,926 footprints without a recorded height carry 5 m) |
| GSHHG shoreline, full resolution | https://www.soest.hawaii.edu/pwessel/gshhg/ | `data/external/gshhg/` (the Oʻahu polygon extracted from it by `analysis/inputs/extract_island_polygon.R` is already in `oahu_gshhs_f.geojson`) |

The other external files are small and included; `analysis/inputs/` holds the scripts that
retrieved them (the weather series, the site elevations), so every input can be traced to
its source or rebuilt from it.

### The Python replication

`replication/python/` implements the same pipeline a second time (pandas, rasterio, Shapely,
statsmodels `MixedLM`) and writes to `replication/results/`; `04_replication_check.R` compares
the two results files quantity by quantity (`replication/replication_check.csv`): 117 paired
quantities and 25 maximum-difference checks. The ΔT values agree to machine precision, every
coefficient and variance component to three decimals and every AIC and R² to two; the few
larger differences are documented implementation choices (Satterthwaite versus Wald tests,
boundary pixels of the raster circles, a collinear night dummy) or optimizer tolerance.

```
pip install -r replication/python/requirements.txt
python replication/python/00_prepare_inputs.py
python replication/python/01_site_predictors.py     # needs the large inputs
python replication/python/02_models.py
```

## Data sources

Logger readings: Thermochron iButton DS1921H loggers (0.5 °C resolution, 2,048-reading memory) in
polystyrene-cup radiation shields at about 1.5 m on public sign posts, poles and fences, read out at
recovery with the manufacturer's application; every logger returned exactly 2,048 readings, the first
42.7 days of its mission for 72 loggers and the last 42.7 days for two whose memory rolled over
(`logger_metadata.csv`). The positions in `sites.csv` are those of the posts. Weather: hourly 10 m
wind speed and total cloud cover retrieved through the Open-Meteo historical weather API
(Zippenfenig 2023) for the centroid of each district; Open-Meteo's default historical dataset
combines the ERA5 and ERA5-Land reanalyses with the operational high-resolution analysis of the
ECMWF Integrated Forecasting System (IFS, 9 km grid), and for this period and place the values
are those of the IFS analysis (verified hour by hour against each dataset the API offers), at
the grid points 21.34°N 157.80°W (Honolulu) and 21.34°N 158.07°W (ʻEwa).
`analysis/inputs/fetch_weather_series.py` retrieves the series and reproduces the committed
file byte for byte. Surface: NOAA C-CAP 2021
high-resolution land cover (Hawaiʻi); FEMA/ORNL USA Structures; GSHHG full-resolution shoreline
(Wessel & Smith 1996). Terrain: USGS 3DEP 1 m DEM via the Elevation Point Query Service
(https://epqs.nationalmap.gov, queried 19 September 2026) for the site elevations; the USGS 10 m
DEM of Oʻahu served by PacIOOS ERDDAP (dataset `usgs_dem_10m_oahu`, subsampled, retrieved
26 September 2026) for the terrain maps. Sea: NOAA CO-OPS station 1612340 hourly water
temperature (retrieved 26 September 2026). Base maps: Natural Earth 10 m land and roads.

## Licence and citation

The code (`analysis/`, `replication/python/`) is released under the MIT License (`LICENSE`).
The measurements, processed data, results and figures (`data/`, `results/`,
`replication/results/`) are released under the Creative Commons Attribution 4.0 International
licence (`LICENSE-DATA.md`); the third-party inputs keep their providers' terms. If you use the
data or code, please cite the thesis (`CITATION.cff`):

> Kellum, J. (2027). *Island Urban Climatology: Measuring and Modeling Nocturnal Heat in Leeward
> Coastal Oʻahu, Hawaiʻi*. M.A. thesis, University of Hawaiʻi at Mānoa.
> https://github.com/justinkellum-ship-it/JustinKellumThesis
