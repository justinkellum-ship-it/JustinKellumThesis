#!/usr/bin/env Rscript
# 00a_compile_logger_readings.R ------------------------------------------------
# The first link of the chain: from the loggers' export tables to one record.
#
# Input: data/raw/exports/<district>_export.csv, one table per district as downloaded
# from the iButton loggers with the manufacturer's application: one row per reading
# (logger, date, time, temperature in °C) plus one row per logger that carries the
# application's metadata as JSON (model, serial number, sample rate, mission start,
# sample counts).  Nothing in these tables has been edited.
#
# Output:
#   data/raw/logger_readings.csv   every reading: district, logger, time stamp (HST), °C
#   data/raw/logger_metadata.csv   one row per logger: model, serial number, sample rate,
#                                  mission start, sample counts, range of the record
#
# The script also compares the readings with the stored logger_readings.csv when it
# already exists, so that the record used by the analysis can be checked against the
# exports at any time (the two must agree reading for reading).
#
# Packages: readr (CSV), dplyr (grouping), jsonlite (the metadata rows).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8")

EXPORTS <- c(Honolulu = file.path(RAW, "exports", "honolulu_export.csv"),
             Ewa      = file.path(RAW, "exports", "ewa_export.csv"))

read_export <- function(region, path) {
  x <- read_csv(path, col_types = cols(.default = col_character()), na = character(), show_col_types = FALSE)
  stopifnot(identical(names(x), c("sensor_id", "date", "time", "temp_c")))
  is_meta <- startsWith(x$temp_c, "{")                      # the application's metadata row of each logger
  meta <- x[is_meta, ]
  rd <- x[!is_meta, ] %>%
    transmute(region = region, sensor_id, datetime = paste(date, time), temp_c = as.numeric(temp_c))
  stopifnot(!anyNA(rd$temp_c), all(nchar(rd$datetime) == 19))
  md <- bind_rows(lapply(seq_len(nrow(meta)), function(i) {
    app <- fromJSON(meta$date[i]); m <- fromJSON(meta$temp_c[i])
    tibble(region = region, sensor_id = meta$sensor_id[i],
           model = m$model_number, serial_number = m$serial_number,
           sample_rate_s = m$sample_rate, resolution_c = as.numeric(m$deviceMissionParametersObj$missionParameters$device$resolution),
           mission_start = m$start_timestamp, first_sample_utc = m$dataInfoObj$firstSampleTimestamp,
           mission_samples = m$deviceDownloadData$missionSampleCount, device_samples = m$deviceDownloadData$deviceSampleCount,
           rollover = m$deviceDownloadData$rollOver, min_c = m$dataInfoObj$min, max_c = m$dataInfoObj$max,
           app_version = app$version)
  }))
  list(readings = rd, metadata = md)
}

main <- function() {
  parts <- lapply(names(EXPORTS), function(r) read_export(r, EXPORTS[[r]]))
  rd <- bind_rows(lapply(parts, `[[`, "readings")) %>%
    mutate(region = factor(region, levels = names(EXPORTS))) %>%
    arrange(region, sensor_id, datetime) %>% mutate(region = as.character(region))
  md <- bind_rows(lapply(parts, `[[`, "metadata")) %>% arrange(match(region, names(EXPORTS)), sensor_id)
  # every logger has one metadata row and its readings; no reading is duplicated
  stopifnot(!any(duplicated(rd[c("sensor_id", "datetime")])),
            setequal(md$sensor_id, unique(rd$sensor_id)), !any(duplicated(md$sensor_id)))
  counts <- rd %>% group_by(region, sensor_id) %>%
    summarise(readings = n(), first = min(datetime), last = max(datetime), .groups = "drop")
  md <- md %>% left_join(counts, by = c("region", "sensor_id"))
  cat(sprintf("%d readings from %d loggers (%s); %d metadata rows\n", nrow(rd), n_distinct(rd$sensor_id),
              paste(sprintf("%s %d", names(table(rd$region)), table(rd$region)), collapse = ", "), nrow(md)))
  print(md %>% count(model, sample_rate_s, resolution_c, rollover))

  out <- file.path(RAW, "logger_readings.csv")
  if (file.exists(out)) {                                    # the stored record must agree reading for reading
    old <- read_csv(out, col_types = cols(datetime = col_character(), .default = col_guess()), show_col_types = FALSE)
    j <- full_join(rd, old, by = c("region", "sensor_id", "datetime"), suffix = c("", "_stored"))
    stopifnot(nrow(j) == nrow(rd), nrow(j) == nrow(old), all(abs(j$temp_c - j$temp_c_stored) < 1e-9))
    cat("the stored logger_readings.csv agrees with the exports reading for reading\n")
  }
  write_csv(rd, out, na = "")
  write_csv(md, file.path(RAW, "logger_metadata.csv"), na = "")
  invisible(rd)
}

if (sys.nframe() == 0) main()
