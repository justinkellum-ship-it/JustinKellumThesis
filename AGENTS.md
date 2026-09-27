# Running this analysis on another sensor network

Instructions for an AI coding agent (Claude Code, Codex, Gemini CLI, Cursor, or similar) working
with a person who has their own temperature-logger network, and for the person directing it.
The repository contains everything needed to repeat the thesis analysis on new data: the scripts
in `analysis/`, worked on the Oʻahu data in `data/`, with every processing rule explained in
`walkthrough/Analysis_walkthrough.html`. Read this file first, then the README, then the
walkthrough. Do not start editing before the intake below is complete.

## What the pipeline does

From the readings of a network of fixed air-temperature loggers, it (1) puts every logger on a
common time grid and computes each site's departure from the median of its district in every
bin (ΔT), (2) selects calm, clear nights from an hourly wind and cloud series, (3) describes
each site by eight properties of its surroundings (impervious, tree-canopy, building-footprint
and water fractions, mean building height, canyon aspect ratio, sky view factor, distance to
the coast) at 50, 100 and 200 m from land-cover rasters, building footprints and a coastline,
and (4) explains the site-mean ΔT with linear mixed-effects models (random intercept per
site), including scale selection, single-descriptor and full models, district comparisons,
best subsets, threshold sensitivity, cross-validation and site effects. Every threshold is a
named constant, every step writes a table, and the results file `results/thesis_results.json`
holds every statistic.

## Intake: settle these with the person before touching anything

Record the answers in your report; each one maps to a constant or a file below.

1. **Logging interval.** 15, 30 or 60 minutes (or another value). → `BIN_MINUTES`.
2. **Time stamps.** Local standard time, without daylight-saving shifts, in the form
   `YYYY-MM-DD HH:MM:SS`. If the exports are in UTC or a zone with DST, convert them once,
   before anything else, and say so; the pipeline applies no offsets.
3. **Night window.** Default 18:00–06:00 local. → `NIGHT_START`, `NIGHT_END`.
4. **Districts.** The names of the networks and which is which. The model comparison steps
   (district term, interactions, shared nights) are written for two districts; with one
   district the district term must be dropped from `02_models.R`, with three or more the
   district steps must be adapted. Say which case applies before running `02_models.R`.
5. **Sites.** One row per logger: `sensor_id, region, latitude, longitude` (WGS84 decimal
   degrees), optionally `setting`. → `data/raw/sites.csv`.
6. **Readings.** Either the loggers' export tables, one per district, named
   `<district>_export.csv` with columns `sensor_id, date, time, temp_c` (a metadata row per
   logger, as the iButton application writes it, is allowed) → `data/raw/exports/`; or one
   compiled file `region, sensor_id, datetime, temp_c` → `data/raw/logger_readings.csv`.
   Temperatures in °C (convert from °F once, and say so). Never edit, fill or drop readings.
7. **Weather.** An hourly series per district for the whole deployment: 10 m wind speed
   (km h⁻¹) and total cloud cover (%), columns `datetime, wind, cloud, region` →
   `data/external/era5_hourly_series.csv`. ERA5 through the Open-Meteo historical-weather
   API gives this anywhere (`hourly=wind_speed_10m,cloud_cover`, `wind_speed_unit=kmh`, one
   request per district at its centroid, `timezone` set to the local standard zone). If no
   series can be obtained, stop and ask; do not invent one.
8. **Land cover.** Class rasters (GeoTIFF tiles) for impervious surface, tree canopy and
   water, at 10 m resolution or finer, in a projected CRS in metres, in one folder under
   `data/external/`. Each class is a set of tiles matched by a pattern in the file name and a
   pixel value that marks the class. Sources: NOAA C-CAP high-resolution (US coasts), a
   state or city 1 m product, ESA WorldCover (10 m; tree = 10, water = 80, impervious =
   built-up 50). → `LANDCOVER_DIR`, `LANDCOVER`, `CRS_M`.
9. **Building footprints.** A GeoPackage layer of polygons with a `height_m` column (NA
   allowed). Sources: FEMA/ORNL USA Structures (US, with heights), Microsoft Global
   Building Footprints or Overture Maps (heights partly), OpenStreetMap. Buildings without a
   height receive `HEIGHT_FILL_M`; the thesis used 5 m and checks 3 m. → `BUILDINGS_FILE`.
10. **Coastline.** A polygon of the land (GSHHG full resolution, or OSM). If the network is
    far from any coast, keep the descriptor out: remove `coast_km` from `PRED`, `LABEL` and
    `UNIT` in `helpers.R` and note it. → `COAST_FILE`.
11. **Radii.** Default 50, 100 and 200 m; `02_models.R` adopts one of them by AIC. → `RADII`.
12. **Thresholds.** Calm = wind < 10 km h⁻¹, clear = cloud < 25 %, in at least 75 % of the
    night bins; at least 20 loggers reporting; a sensor-night needs 75 % of its bins. Keep
    them unless the person has a reason to change them; the threshold-sensitivity step in
    `02_models.R` shows what changing them does. Any change goes in `helpers.R` and in the
    report, never in a script.

## Set the configuration

All of it is the "study configuration" block near the top of `analysis/helpers.R`:
`BIN_MINUTES`, `NIGHT_START`, `NIGHT_END`, the thresholds, `RADII`, `CRS_M`,
`LANDCOVER_DIR`, `LANDCOVER`, `BUILDINGS_FILE`, `HEIGHT_FILL_M`, `HEIGHT_RANGE_M`,
`COAST_FILE`. The district names are read from `sites.csv`. Nothing else in the scripts
repeats a value that lives there; if you find one, that is a bug to report, not to work around.

## Run in order, and check each step before the next

```
Rscript analysis/requirements.R                 # installs the packages (R >= 4.3)
Rscript analysis/00a_compile_logger_readings.R  # only with export tables
Rscript analysis/00_prepare_inputs.R
Rscript analysis/01_site_predictors.R           # minutes per hundred sites
Rscript analysis/02_models.R
```

After `00a`: readings per logger, first and last stamp per logger, no duplicated
(logger, stamp) pair; the script stops by itself if a stored record disagrees with the exports.
After `00`: print `data/processed/night_inventory.csv` (every night with its decision) and the
count of sensor-nights per district; if no night is selected, the weather series or the
thresholds are wrong for this climate, so stop and discuss rather than loosening rules
silently. After `01`: descriptors within range (fractions 0–100, sky view factor 0–1, heights
plausible for the place), `n_bldg` above zero at built sites, `coast_km` sensible; read the
walkthrough's Step 2 if a value looks wrong. After `02`: no convergence warnings, the adopted
radius, the marginal R² of the full models, and `results/thesis_results.json` written.

`03_figures.R`, `05_thesis_tables.R`, `06_elevation_check.R`, `07_sensitivity_checks.R` and
`08_appendix_figures.R` are written for the Oʻahu study (the district names and site counts
in their labels, one Oʻahu site in check A, a Hawaiʻi DEM); run them only after adapting
them, and say what was adapted. The CSV tables and the results file from `02_models.R` hold
every statistic they format.
`04_replication_check.R` needs the Python implementation in `replication/` to have been run
on the same inputs; it is optional.

## Rules

Raw files are never edited. No reading is filled, smoothed or dropped; the completeness rules
decide what enters the analysis, and they are reported. No weather value is invented. Every
constant that differs from the thesis is listed in the report with the reason. Results are
reported with `results/sessionInfo.txt` (package versions) and the exact commands run.

## What to hand back

The configuration used (the block from `helpers.R`), the night inventory with the number of
sensor-nights per district, `data/descriptors/site_predictors_multiscale.csv`,
`results/thesis_results.json` with `results/tables/*.csv`,
and a short account of anything that had to be adapted. If you also want the figures of the
thesis, adapt `03_figures.R` next; its maps expect the same land-cover classes and footprints.
