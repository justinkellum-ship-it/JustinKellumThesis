# Nocturnal heat in leeward coastal Oʻahu: sensor network, data and analysis

Data, code and results of the M.A. thesis *Island Urban Climatology: Measuring and Modeling
Nocturnal Heat in Leeward Coastal Oʻahu, Hawaiʻi* (Justin Kellum, Department of Geography and
Environment, University of Hawaiʻi at Mānoa, 2027).

One hundred low-cost iButton stations recorded canopy-layer air temperature every 30 minutes in
two districts of leeward Oʻahu between November 2024 and January 2025: the compact, partly
high-rise core of Honolulu (36 recovered sites) and the low-rise suburban plain of ʻEwa–Kapolei
(38 sites). On the 11 calm, clear nights identified from ERA5 (400 sensor-nights), each site's
night-mean departure from its district network median (ΔT) is related to eight descriptors of
urban form and setting, computed at 50, 100 and 200 m from 1 m NOAA C-CAP land cover, FEMA
building footprints and the GSHHG shoreline, with linear mixed-effects models (random intercept
per site). Every number, table and figure in the thesis is produced from the files in this
repository by a single command.

**Project page:** https://justinkellum-ship-it.github.io/JustinKellumThesis/ (thesis, figures,
and the step-by-step walkthrough of the analysis)

## Layout

```
data/
  raw/           the study's own measurements
                   exports/              the loggers' export tables as downloaded (one per district: every reading,
                                         plus one metadata row per logger from the download application)
                   logger_readings.csv   the 151,552 readings compiled from them: district, logger, time stamp (HST), °C
                   logger_metadata.csv   one row per logger: model, serial number, sample rate, mission start, sample counts
                   sites.csv             the 74 loggers: district, WGS84 position, setting recorded at deployment
  external/      data from other providers (small files included; large ones downloaded, see below)
                   era5_hourly_series.csv            ERA5 10 m wind and cloud cover, one grid point per district
                   oahu_gshhs_f.geojson              Oʻahu polygon from the GSHHG full-resolution shoreline
                   ne_10m_land_hawaii.geojson, ne_10m_roads_oahu.geojson   Natural Earth base-map layers
                   oahu_dem_500m_q5.txt, honolulu_dem_150m_q5.txt         USGS 10 m DEM subsamples (Appendix B)
                   noaa_1612340_water_temperature.csv                     Honolulu tide-station water temperature
  processed/     written by analysis/00_prepare_inputs.R
                   night_inventory.csv         every district × night with the inclusion decision
                   sensor_nights.csv           the 400 sensor-nights: the unit of analysis
                   halfhourly_calm_clear.csv   the half-hourly record of the analysed diurnal cycles
  descriptors/   written by analysis/01_site_predictors.R and the checks
                   site_predictors_multiscale.csv   the eight descriptors of every site at 50, 100 and 200 m
                   site_elevation.csv, site_default_height_share.csv, island_default_height_share.csv,
                   site_height_fill3m_100m.csv, site_svf_footprint_check.csv
analysis/        the R pipeline, numbered in the order it runs (below)
results/         thesis_results.json (every statistic in the thesis), tables/*.csv, figures/ (PNG at 300 dpi,
                 plus PDF for the statistical figures),
                 thesis_tables/*.md (the 22 tables as printed), elevation_check.json,
                 sensitivity_checks.json, island_landcover.json, sessionInfo.txt
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
| `00_prepare_inputs.R` | snaps the loggers' clocks to 30-minute bins, computes the district network median and ΔT, joins the ERA5 series, classifies every night, writes the sensor-night table | seconds |
| `01_site_predictors.R` | the eight descriptors of every site at 50, 100 and 200 m: impervious, tree-canopy, building-footprint and water fractions, mean building height, canyon aspect ratio, ray-cast sky view factor, distance to the coast | ~7 min |
| `02_models.R` | every mixed-effects model: scale selection, collinearity, null models, single descriptors, full pooled and district models, best subsets and Akaike weights, district × descriptor interactions, night effect, threshold sensitivity, model comparison, site effects and leave-one-site-out cross-validation, shared nights | ~1 min |
| `03_figures.R` | the statistical figures and maps of the main text | ~10 min |
| `04_replication_check.R` | compares the R results with the Python implementation | seconds |
| `05_thesis_tables.R` | the 22 tables of the thesis, in Markdown, as printed | seconds |
| `06_elevation_check.R` | the terrain check: site elevation (USGS 3DEP) added to the models, best subsets and cross-validation | ~1 min |
| `07_sensitivity_checks.R` | three site-level checks: without the warmest site; with a different height for the buildings whose height the inventory does not record; with the sky view factor of positions inside footprints recomputed | ~7 min |
| `08_appendix_figures.R` | the land-cover and terrain maps and the sky view factor diagram of Appendix B | ~1 min |

`helpers.R` holds the paths, the processing constants (every threshold is a named constant) and
the model helpers; `maplib.R` the map layers. The walkthrough shows the code of each step with
its output.

### Processing rules (constants in `analysis/helpers.R`, implemented in `00_prepare_inputs.R`)

1. Readings snapped to the nearest 30-minute bin (the loggers' clocks differ by up to 2 min).
2. Reference = median of all loggers of the district reporting in the bin; ΔT = T − median.
3. Night = 18:00–06:00 HST; a diurnal date starts at 06:00.
4. Calm and clear night: ERA5 10 m wind < 10 km/h and cloud < 25 % in ≥ 75 % of the 24 night bins; all 24 ERA5 bins present.
5. ≥ 20 loggers reporting; a sensor-night needs ≥ 18 of its 24 bins.

## Reproducing the results

```
Rscript analysis/requirements.R     # installs the packages once (CRAN)
Rscript analysis/run_all.R          # ~30 minutes: data/processed, data/descriptors, results
Rscript analysis/04_replication_check.R   # optional: compare with replication/results
```

Requires R ≥ 4.3. The versions used for the thesis are in `results/sessionInfo.txt` (R 4.3;
lme4, lmerTest, performance, car, sf, terra, ggplot2, patchwork, dplyr, tidyr, readr, lubridate,
jsonlite). `00_prepare_inputs.R` and everything that reads the stored descriptors run from the
files in the repository; recomputing the descriptors, the maps and the site checks needs three
large inputs placed in `data/external/`:

| Input | Where | Put it at |
|---|---|---|
| NOAA C-CAP 2021 high-resolution (1 m) land cover for Hawaiʻi: the impervious, canopy and water masks for Oʻahu | https://coast.noaa.gov/digitalcoast/data/ccaphighres.html | `data/external/ccap/*.tif` |
| FEMA / ORNL USA Structures, Hawaiʻi | https://gis-fema.hub.arcgis.com/pages/usa-structures | `data/external/buildings_oahu.gpkg` (the Oʻahu footprints with the HEIGHT attribute; missing heights filled with 5 m) |
| GSHHG shoreline, full resolution | https://www.soest.hawaii.edu/pwessel/gshhg/ | `data/external/gshhg/` (the Oʻahu polygon extracted from it is already in `oahu_gshhs_f.geojson`) |

### The Python replication

`replication/python/` implements the same pipeline a second time (pandas, rasterio, Shapely,
statsmodels `MixedLM`) and writes to `replication/results/`; `04_replication_check.R` compares
the two results files quantity by quantity (`replication/replication_check.csv`). The ΔT values
agree to machine precision and every coefficient, variance component, AIC and R² to three
decimals; the few larger differences are documented implementation choices (Satterthwaite
versus Wald tests, boundary pixels of the raster circles, a collinear night dummy).

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
(`logger_metadata.csv`). The positions in `sites.csv` are those of the posts. Weather: ERA5 (Hersbach et al. 2020) hourly 10 m wind speed
and total cloud cover via the Open-Meteo archive interface. Surface: NOAA C-CAP 2021
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
