#!/usr/bin/env Rscript
# 00_prepare_inputs.R ----------------------------------------------------------
# From the complete half-hourly iButton record (data/raw/logger_readings.csv:
# every reading of every logger on both deployments), the ERA5 hourly wind and
# cloud series for the two districts (data/external/era5_hourly_series.csv) and the
# logger positions (data/raw/sites.csv), build the three analysis tables:
#
#   data/processed/night_inventory.csv        one row per district x diurnal date (06:00-06:00
#                                   HST): loggers reporting, ERA5 night means, share of
#                                   calm & clear bins, and the inclusion decision
#   data/processed/halfhourly_calm_clear.csv  one row per logger x 30-min bin for the diurnal
#                                   cycles that contain a qualifying night
#   data/processed/sensor_nights.csv          one row per logger x qualifying night
#                                   (18:00-06:00 means) - the unit of analysis
#
# Processing rules (each is a named constant in helpers.R so it can be varied):
#   1. Time binning.  The loggers were started by hand, so their clocks differ by
#      up to +-2 min (stamps :00/:01/:02 and :30/:31/:32).  Every reading is
#      assigned to the nearest 30-min bin so that all loggers of a district are
#      compared at the same time step (a median taken by exact minute stamp would
#      compare each logger only with the loggers that share its cadence).
#   2. Reference.  The district network median air temperature of all loggers
#      reporting in the bin; dT = T_logger - T_median.
#   3. Night.  18:00-06:00 HST; a diurnal date runs 06:00-06:00 and is labelled by
#      the date of its first 06:00.
#   4. Calm & clear night.  ERA5 10 m wind < 10 km/h AND total cloud cover < 25 %
#      in at least 75 % of the 24 night bins; all 24 ERA5 bins must exist.
#   5. Network completeness.  A night is analysed only if >= 20 loggers reported.
#   6. Sensor-night completeness.  A sensor-night needs >= 18 of its 24 bins.
#
# Packages: dplyr/tidyr (grouped summaries and joins), lubridate (round_date),
#           readr (fast CSV input/output).
# ------------------------------------------------------------------------------
# locate helpers.R next to this script (or under THESIS_ROOT/analysis when knitted)
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")

# ---- 1-2. binned record with the network median -------------------------------
load_binned <- function() {
  hh <- read_csv(file.path(RAW, "logger_readings.csv"), col_types = cols(datetime = col_character(), .default = col_guess())) %>%
    mutate(datetime = read_stamps(datetime),
           time_bin = round_date(datetime, "30 minutes"),
           hour_dec = hour(time_bin) + minute(time_bin) / 60,
           night_date = format(time_bin - hours(NIGHT_END), "%Y-%m-%d"),
           is_night = hour_dec >= NIGHT_START | hour_dec < NIGHT_END,
           hour = hour(time_bin)) %>%
    # one reading per logger and bin (duplicates would only arise from clock resets)
    arrange(region, sensor_id, datetime) %>%
    distinct(region, sensor_id, time_bin, .keep_all = TRUE)

  med <- hh %>%
    group_by(region, time_bin) %>%
    summarise(network_median_temp = median(temp_c), n_sensors_active = n_distinct(sensor_id), .groups = "drop")
  hh <- hh %>% inner_join(med, by = c("region", "time_bin")) %>%
    mutate(dT_network = temp_c - network_median_temp)

  # ERA5 hourly values were attached to the logger time stamps by nearest hour,
  # the half-hour stamp being assigned to the following hour (a :31 reading
  # carries the value of the next full hour).  The bin value is the one carried
  # by the :01 / :31 stamps of the main logger group.
  era <- read_csv(file.path(EXT, "era5_hourly_series.csv"), col_types = cols(datetime = col_character(), .default = col_guess())) %>%
    mutate(datetime = read_stamps(datetime),
           time_bin = round_date(datetime, "30 minutes"),
           pref = as.integer(minute(datetime) %% 30 != 1)) %>%          # :01/:31 rows first
    arrange(region, time_bin, pref) %>%
    distinct(region, time_bin, .keep_all = TRUE) %>%
    transmute(region, time_bin, era5_wind_speed_kmh = wind, era5_cloud_cover_pct = cloud)
  hh <- hh %>% left_join(era, by = c("region", "time_bin"))
  list(hh = hh, era = era)
}

# ---- 3-5. classify every diurnal date of each district ------------------------
night_inventory <- function(hh, era, wind_max = WIND_MAX_KMH, cloud_max = CLOUD_MAX_PCT,
                            min_frac = MIN_FRAC_STEPS, min_sensors = MIN_SENSORS_PER_NIGHT) {
  en <- era %>%
    mutate(hour_dec = hour(time_bin) + minute(time_bin) / 60,
           night_date = format(time_bin - hours(NIGHT_END), "%Y-%m-%d"),
           is_night = hour_dec >= NIGHT_START | hour_dec < NIGHT_END) %>%
    filter(is_night) %>%
    mutate(ok = era5_wind_speed_kmh < wind_max & era5_cloud_cover_pct < cloud_max)
  inv <- en %>% group_by(region, night_date) %>%
    summarise(n_era5_bins = n(), frac_calm_clear = mean(ok),
              wind_night = mean(era5_wind_speed_kmh), cloud_night = mean(era5_cloud_cover_pct),
              frac_calm = mean(era5_wind_speed_kmh < wind_max), frac_clear = mean(era5_cloud_cover_pct < cloud_max),
              .groups = "drop")
  sens <- hh %>% filter(is_night) %>% group_by(region, night_date) %>%
    summarise(n_sensors = n_distinct(sensor_id), n_bins = n_distinct(time_bin),
              T_median_night = mean(network_median_temp), .groups = "drop")
  inv <- full_join(inv, sens, by = c("region", "night_date")) %>%
    mutate(n_sensors = coalesce(n_sensors, 0L),
           meets_weather = coalesce(frac_calm_clear >= min_frac & n_era5_bins == 24, FALSE),
           meets_network = n_sensors >= min_sensors,
           selected = meets_weather & meets_network,
           decision = case_when(
             selected ~ "selected",
             !meets_weather & coalesce(n_era5_bins < 24, FALSE) ~ "partial night (deployment/retrieval)",
             !meets_weather ~ "not calm and clear",
             TRUE ~ sprintf("fewer than %d sensors reporting", min_sensors))) %>%
    arrange(region, night_date)
  inv
}

# ---- 6. sensor-night means ----------------------------------------------------
sensor_nights <- function(hh, inv, min_steps = MIN_STEPS_PER_SENSOR_NIGHT) {
  keep <- inv %>% filter(selected) %>% select(region, night_date)
  sn <- hh %>% filter(is_night) %>% inner_join(keep, by = c("region", "night_date")) %>%
    group_by(region, sensor_id, night_date) %>%
    summarise(dT_night = mean(dT_network), dT_night_median = median(dT_network),
              temp_night = mean(temp_c), ref_night = mean(network_median_temp),
              wind_night = mean(era5_wind_speed_kmh, na.rm = TRUE), cloud_night = mean(era5_cloud_cover_pct, na.rm = TRUE),
              n_steps = n(), .groups = "drop") %>%
    filter(n_steps >= min_steps) %>%
    group_by(region) %>% mutate(night_no = dense_rank(night_date)) %>% ungroup() %>%
    mutate(night_id = paste0(substr(region, 1, 3), "_", night_no))
  sn
}

main <- function() {
  b <- load_binned(); hh <- b$hh; era <- b$era
  inv <- night_inventory(hh, era)
  write_csv(inv, file.path(PROC, "night_inventory.csv"), na = "")
  print(as.data.frame(inv %>% filter(selected | frac_calm_clear >= MIN_FRAC_STEPS)))

  keep <- inv %>% filter(selected) %>% select(region, night_date)
  cyc <- hh %>% inner_join(keep, by = c("region", "night_date")) %>%
    mutate(is_calm = era5_wind_speed_kmh < WIND_MAX_KMH, is_clear = era5_cloud_cover_pct < CLOUD_MAX_PCT) %>%
    arrange(region, sensor_id, datetime) %>%
    transmute(region, sensor_id, datetime = fmt_stamp(datetime), time_bin = fmt_stamp(time_bin), night_date, hour,
              hour_dec, is_night, temp_c, network_median_temp, n_sensors_active, dT_network,
              era5_wind_speed_kmh, era5_cloud_cover_pct, is_calm, is_clear)
  write_csv(cyc, file.path(PROC, "halfhourly_calm_clear.csv"), na = "")

  sn <- sensor_nights(hh, inv)
  coords <- read_csv(file.path(RAW, "sites.csv"), show_col_types = FALSE) %>% select(sensor_id, latitude, longitude)
  sn <- sn %>% left_join(coords, by = "sensor_id") %>% arrange(region, sensor_id, night_date)
  write_csv(sn, file.path(PROC, "sensor_nights.csv"), na = "")
  print(sn %>% group_by(region) %>% summarise(sensors = n_distinct(sensor_id), nights = n_distinct(night_date), rows = n()))
  print(sn %>% count(region, night_date), n = 20)
  invisible(sn)
}

if (sys.nframe() == 0) main()
