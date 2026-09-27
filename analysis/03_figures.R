#!/usr/bin/env Rscript
# 03_figures.R -----------------------------------------------------------------
# The thesis figures, drawn with ggplot2 (statistics) and sf/terra (maps).
# Every figure has one stated purpose (in the comment of its function and in
# the thesis caption).  Maps use the 1 m NOAA C-CAP masks and the FEMA building
# footprints as base map (maplib.R).
#
#   T01 study area, three panels (A Oʻahu, B Honolulu, C ʻEwa)
#   T02 methodological workflow
#   T03 night selection (ERA5 wind / cloud, network median, selected nights)
#   T04 why calm and clear: spread of dT across the network vs wind and cloud
#   T05 how dT is defined (one shared night, both districts)
#   T06 night-time temperature distributions, all vs calm, clear nights
#   T07 predictor maps, Honolulu (8 panels)     T08 predictor maps, ʻEwa
#   T09 scale selection (delta AIC by radius; single-predictor R2m by radius)
#   T10 correlation matrices with VIF
#   T11 bivariate site-mean dT vs predictors
#   T12 fixed-effect coefficients (pooled, Honolulu, ʻEwa)
#   T13 district x predictor interactions (height, coast)
#   T14 threshold sensitivity
#   T15 BLUPs (sorted) and BLUP maps
#   T16 site-mean dT maps
#   T17 shared nights: absolute network medians
#   T18 leave-one-site-out cross-validation
#   T19 night-to-night stability of the site pattern
#
# Usage: Rscript 03_figures.R [01 02 ...]   (no argument = all figures)
#
# Packages: ggplot2 (all plots), patchwork (multi-panel layout), sf + terra
#           (maps), dplyr/tidyr (data shaping), jsonlite (results file),
#           scales (colour scales), lubridate (dates).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8"); source(file.path(.here, "maplib.R"), encoding = "UTF-8")
suppressPackageStartupMessages({ library(scales) })

RES <- fromJSON(file.path(OUT, "thesis_results.json"), simplifyVector = TRUE)
RADIUS <- RES$adopted_radius
PLAB <- c(imperv = "Impervious surface (%)", tree = "Tree canopy (%)", bldg = "Building footprint (%)",
          water = "Water surface (%)", height = "Mean building height (m)", svf_point = "Sky view factor",
          aspect = "Canyon aspect ratio (H/W)", coast_km = "Distance to coast (km)")
PSHORT <- c(imperv = "Impervious", tree = "Tree canopy", bldg = "Bldg footprint", water = "Water",
            height = "Bldg height", svf_point = "Sky view", aspect = "Aspect ratio", coast_km = "Coast dist.")
reg_f <- function(x) factor(x, levels = REGIONS, labels = REG_LABEL[REGIONS])
scale_reg <- function(...) list(scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL, ...),
                                scale_fill_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL, ...))
ptitle <- function(s) ggtitle(s)

island <- load_island()
sites  <- load_sites(RADIUS)
DOM <- list(Honolulu = domain_bounds(sites, "Honolulu", 600), Ewa = domain_bounds(sites, "Ewa", 600))
EQB <- equal_bounds(DOM)
.BLD <- list()
bld <- function(region) { if (is.null(.BLD[[region]])) .BLD[[region]] <<- load_buildings(EQB[[region]]); .BLD[[region]] }

# ==============================================================================
# T01  Study area
figT01_study_area <- function() {
  # Purpose: locate the two study areas on Oʻahu relative to the Koʻolau and
  # Waiʻanae ranges and the trade-wind flow, and show their contrasting urban
  # form at the same map scale.
  ib <- st_bbox(island)
  bounds <- bounds4(ib$xmin - 1500, ib$ymin - 1500, ib$xmax + 9000, ib$ymax + 1500)
  imp_all <- raster_df(tint_raster("impervious", 1, bounds, 20, "island"))
  can_all <- raster_df(tint_raster("canopy", 1, bounds, 20, "island"))
  imp <- imp_all %>% filter(v > 0.01); can <- can_all %>% filter(v > 0.01)     # drop empty cells (sea, bare land) to save memory
  roads <- tryCatch(st_transform(st_read(file.path(EXT, "ne_10m_roads_oahu.geojson"), quiet = TRUE), CRS_M), error = function(e) NULL)
  boxes <- bind_rows(lapply(names(EQB), function(r) tibble(region = r, xmin = EQB[[r]]["xmin"], xmax = EQB[[r]]["xmax"],
                                                           ymin = EQB[[r]]["ymin"], ymax = EQB[[r]]["ymax"],
                                                           lab = c(Honolulu = "B", Ewa = "C")[r])))
  ranges <- tibble(name = c("KOʻOLAU  RANGE", "WAIʻANAE  RANGE"), lon = c(-157.905, -158.165), lat = c(21.475, 21.475), rot = c(-47, -70))
  ranges <- bind_cols(ranges, as_tibble(t(sapply(seq_len(nrow(ranges)), function(i) lonlat_to_xy(ranges$lon[i], ranges$lat[i]))), .name_repair = ~ c("x", "y")))
  places <- tribble(~name, ~lon, ~lat, ~st,
    "Pearl Harbor", -157.965, 21.355, "italic", "Māmala Bay", -157.93, 21.275, "italic", "Honolulu", -157.905, 21.335, "plain",
    "Kāneʻohe", -157.80, 21.415, "plain", "Kapolei", -158.06, 21.337, "plain", "Wahiawā", -158.02, 21.50, "plain",
    "Kailua", -157.74, 21.395, "plain", "Central\nOʻahu plain", -158.03, 21.44, "italic", "ʻEwa plain", -158.02, 21.31, "italic",
    "PACIFIC\nOCEAN", -158.24, 21.30, "plain")
  places <- bind_cols(places, as_tibble(t(sapply(seq_len(nrow(places)), function(i) lonlat_to_xy(places$lon[i], places$lat[i]))), .name_repair = ~ c("x", "y")))
  a0 <- lonlat_to_xy(-157.60, 21.60); a1 <- lonlat_to_xy(-157.70, 21.52); at <- lonlat_to_xy(-157.635, 21.578)
  pA <- ggplot() +
    theme_map() + theme(panel.background = element_rect(fill = C_OCEAN, colour = NA), legend.position = c(0.86, 0.10),
                        legend.background = element_rect(fill = alpha("white", 0.9), colour = NA), legend.key.size = unit(0.3, "cm"),
                        legend.text = element_text(size = 6), plot.margin = margin(1, 1, 1, 1)) +
    geom_sf(data = st_geometry(island), fill = C_LAND, colour = NA) +
    geom_raster(data = imp, aes(x, y, alpha = pmin(v / 0.8, 1)), fill = "#8a8a8a", interpolate = TRUE) +
    geom_raster(data = can %>% mutate(v = v * 0.85), aes(x, y, alpha = v), fill = "#6fae62", interpolate = TRUE) +
    scale_alpha_identity() +
    geom_sf(data = st_boundary(st_geometry(island)), colour = C_COAST, linewidth = 0.3) +
    { if (!is.null(roads)) geom_sf(data = st_geometry(roads), colour = "#a08c78", linewidth = 0.3, alpha = 0.8) } +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = region), fill = NA, linewidth = 0.7) +
    geom_text(data = boxes, aes(x = xmin + 300, y = ymax - 300, label = lab, colour = region), hjust = 0, vjust = 1, size = 3.2, fontface = "bold") +
    geom_sf(data = sites, aes(colour = region), size = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = REG_COL, guide = "none") +
    geom_text(data = ranges, aes(x, y, label = name, angle = rot), size = 2.6, colour = "#2f5d2f", fontface = "bold") +
    geom_label(data = places, aes(x, y, label = name, fontface = st), size = 2.2, colour = "#333333", label.size = 0,
               fill = alpha("white", 0.6), label.padding = unit(0.05, "lines"), lineheight = 0.85) +
    annotate("segment", x = a0[1], y = a0[2], xend = a1[1], yend = a1[2], colour = "#1f3b57", linewidth = 0.7,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    annotate("text", x = at[1], y = at[2], label = "NE trade winds", angle = -38, size = 2.4, colour = "#1f3b57", fontface = "bold", vjust = 0) +
    scale_bar(bounds, 10000, loc = c(0.05, 0.05), label = "10km") + north_arrow(bounds, loc = c(0.965, 0.80), size = 0.07) +
    panel_tag(bounds, "A  Oʻahu", size = 3.6) +
    # legend built from dummy layers
    geom_point(data = tibble(x = -Inf, y = -Inf, k = factor(c("Impervious surface (C-CAP 2021)", "Tree canopy (C-CAP 2021)", "Major road", "Sensor site"),
                                                             levels = c("Impervious surface (C-CAP 2021)", "Tree canopy (C-CAP 2021)", "Major road", "Sensor site"))),
               aes(x, y, fill = k), shape = 22, size = 2.5, colour = NA) +
    scale_fill_manual(values = c("Impervious surface (C-CAP 2021)" = "#8a8a8a", "Tree canopy (C-CAP 2021)" = "#6fae62",
                                 "Major road" = "#a08c78", "Sensor site" = "#444444"), name = NULL) +
    map_coord(bounds, graticule = FALSE)
  # inset: the Hawaiian islands
  landne <- tryCatch(st_read(file.path(EXT, "ne_10m_land_hawaii.geojson"), quiet = TRUE), error = function(e) NULL)
  isl <- tibble(name = c("Kauaʻi", "Oʻahu", "Maui", "Hawaiʻi"), lon = c(-159.5, -157.75, -156.3, -155.5), lat = c(22.25, 21.85, 20.45, 19.2))
  pI <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = C_OCEAN, colour = "#444444", linewidth = 0.3)) +
    { if (!is.null(landne)) geom_sf(data = st_geometry(landne), fill = "#dedbd2", colour = C_COAST, linewidth = 0.2) } +
    geom_sf(data = st_geometry(st_transform(island, 4326)), fill = "#f0b3ab", colour = "#c0392b", linewidth = 0.3) +
    geom_text(data = isl, aes(lon, lat, label = name), size = 1.9, colour = "#333333") +
    coord_sf(xlim = c(-160.6, -154.6), ylim = c(18.7, 22.4), expand = FALSE, datum = NA)
  pA <- pA + inset_element(pI, left = 0.015, bottom = 0.60, right = 0.255, top = 0.925, align_to = "panel")
  # panels B and C at the same scale
  pBC <- lapply(c("Honolulu", "Ewa"), function(region) {
    b <- EQB[[region]]; s <- sites[sites$region == region, ]; lab <- c(Honolulu = "B", Ewa = "C")[region]
    ggplot() + theme_map() + basemap(b, island, buildings = bld(region), key = region) +
      geom_sf(data = s, colour = "white", fill = REG_COL[region], shape = 21, size = 2.2, stroke = 0.4) +
      place_labels(region, b, size = 2.0) + scale_bar(b, 1000) + north_arrow(b, loc = c(0.94, 0.10)) +
      panel_tag(b, sprintf("%s  %s (%d sites)", lab, REG_LABEL[region], nrow(s))) +
      graticule_breaks() + map_coord(b, graticule = TRUE) +
      theme(panel.border = element_rect(fill = NA, colour = REG_COL[region], linewidth = 0.9), axis.text = element_text(size = 5.5),
            axis.text.y = element_text(angle = 90, hjust = 0.5), plot.margin = margin(2, 2, 2, 2))
  })
  leg <- ggplot(tibble(k = factor(c("Buildings (FEMA)", "Paved surfaces", "Tree canopy", "Water", "Sensor site"),
                                  levels = c("Buildings (FEMA)", "Paved surfaces", "Tree canopy", "Water", "Sensor site")), x = 1:5), aes(x, 1, fill = k)) +
    geom_point(shape = 22, size = 3, colour = NA) +
    scale_fill_manual(values = c("Buildings (FEMA)" = C_BLDG, "Paved surfaces" = C_PAVED, "Tree canopy" = C_TREE, "Water" = C_WATER, "Sensor site" = "#555555"), name = NULL) +
    theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.key.size = unit(0.3, "cm")) +
    guides(fill = guide_legend(nrow = 1))
  leg <- cowplot::get_legend(leg)
  hA <- (bounds["ymax"] - bounds["ymin"]) / (bounds["xmax"] - bounds["xmin"]) * 7.2 * 0.96
  hB <- (EQB$Honolulu["ymax"] - EQB$Honolulu["ymin"]) / (EQB$Honolulu["xmax"] - EQB$Honolulu["xmin"]) * 3.5
  p <- (pA / (pBC[[1]] | pBC[[2]]) / patchwork::wrap_elements(leg)) + plot_layout(heights = c(hA, hB + 0.35, 0.25))
  save_fig(p, "T01_study_area", width = 7.2, height = hA + hB + 0.95)
  # island land-cover numbers for the text (same aggregated mosaics)
  land <- st_geometry(island); land_km2 <- as.numeric(st_area(land)) / 1e6; px_km2 <- 20 * 20 / 1e6
  out <- list(oahu_land_km2 = land_km2, impervious_km2 = sum(imp_all$v) * px_km2, impervious_pct = 100 * sum(imp_all$v) * px_km2 / land_km2,
              canopy_km2 = sum(can_all$v) * px_km2, canopy_pct = 100 * sum(can_all$v) * px_km2 / land_km2)
  for (r in REGIONS) {
    b <- EQB[[r]]; sel <- function(d) d %>% filter(x >= b["xmin"], x <= b["xmax"], y >= b["ymin"], y <= b["ymax"])
    dom_land <- as.numeric(st_area(st_intersection(land, st_as_sfc(st_bbox(b, crs = st_crs(CRS_M)))))) / 1e6
    out[[paste0(r, "_domain_land_km2")]] <- dom_land
    out[[paste0(r, "_impervious_km2")]] <- sum(sel(imp_all)$v) * px_km2
    out[[paste0(r, "_impervious_pct_of_land")]] <- 100 * sum(sel(imp_all)$v) * px_km2 / dom_land
    out[[paste0(r, "_canopy_pct_of_land")]] <- 100 * sum(sel(can_all)$v) * px_km2 / dom_land
  }
  write_json(out, file.path(OUT, "island_landcover.json"), auto_unbox = TRUE, digits = NA, pretty = TRUE)
}

# ==============================================================================
# T02  Workflow
figT02_workflow <- function() {
  # Purpose: give the reader the whole method in one view before the details.
  d <- RES$design
  boxes <- tribble(~x, ~y, ~w, ~h, ~title, ~body, ~fc, ~ec,
    1, 70, 30, 28, "1  Question and response",
    "Which surface and form properties\nexplain where night air is warmer or\ncooler than the district network?\n\nResponse: ΔT = T_site − network\nmedian, 18:00–06:00 HST mean per\nsensor-night", "#f4f6f8", "#1f3b57",
    35, 70, 30, 28, "2  Measurement",
    sprintf("iButton loggers in ventilated cup\nshields, ≈1.5m, 30-min interval\nHonolulu %d sites · ʻEwa %d sites\nNov 2024 – Jan 2025\nRandom-point site selection\n(≥ 50m apart)", d$Honolulu$n_sites, d$Ewa$n_sites), "#f4f6f8", "#1f3b57",
    69, 70, 30, 28, "3  Processing",
    "Readings snapped to 30-min bins\nNetwork median per bin and district\nERA5 wind and cloud joined\nNight = 18:00–06:00 (diurnal date\nstarts 06:00)", "#f4f6f8", "#1f3b57",
    69, 37, 30, 28, "4  Night selection",
    sprintf("Calm and clear: wind < 10km h⁻¹\nand cloud < 25%% in ≥ 75%% of night\nbins; ≥ 20 sensors; ≥ 18 of 24 bins\nper sensor-night\n→ %d Honolulu + %d ʻEwa nights,\n%d sensor-nights", d$Honolulu$n_nights, d$Ewa$n_nights, d$Honolulu$n_sensor_nights + d$Ewa$n_sensor_nights), "#f4f6f8", "#1f3b57",
    35, 37, 30, 28, "5  Spatial predictors",
    "C-CAP 1m: impervious, tree, water\nFEMA footprints: building fraction,\nheight, canyon H/W, sky view factor\nDistance to the shoreline (GSHHG\ncoastline database)\nRadii 50 / 100 / 200m", "#eef5ee", "#2f5d2f",
    1, 37, 30, 28, "6  Scale choice",
    sprintf("Full model at each radius (ML AIC)\npooled, Honolulu, ʻEwa\n→ %dm adopted (best or within\n2 AIC in all three)\nCollinearity screen (r, VIF)", RADIUS), "#eef5ee", "#2f5d2f",
    1, 4, 30, 28, "7  Mixed-effects models (lme4)",
    "Random intercept per site\nNull → single predictors → full\n(8 predictors + district)\nDistrict models, best subsets\nDistrict × predictor interactions", "#f9f1e7", "#8a4b08",
    35, 4, 30, 28, "8  Checks",
    "Threshold sensitivity (re-select\nnights, refit); calm-core refit\nTerrain (elevation) check\nSite checks: HNL08, missing\nheights, sky view factor\nLeave-one-site-out CV", "#f9f1e7", "#8a4b08",
    69, 4, 30, 28, "9  Interpretation",
    "Effects in natural units (per 10 pp\nimpervious, per m of height)\nBLUP maps of unexplained heat\nAbsolute comparison on the two\nshared nights", "#f9f1e7", "#8a4b08")
  arrows <- tribble(~x0, ~y0, ~x1, ~y1, 31, 84, 35, 84, 65, 84, 69, 84, 84, 70, 84, 65.5, 69, 51, 65, 51, 35, 51, 31, 51,
                    16, 37, 16, 32.5, 31, 18, 35, 18, 65, 18, 69, 18)
  p <- ggplot() + theme_void() +
    geom_rect(data = boxes, aes(xmin = x, xmax = x + w, ymin = y, ymax = y + h, fill = fc, colour = ec), linewidth = 0.4) +
    scale_fill_identity() + scale_colour_identity() +
    geom_text(data = boxes, aes(x = x + w / 2, y = y + h - 2.2, label = title, colour = ec), vjust = 1, size = 3.0, fontface = "bold") +
    geom_text(data = boxes, aes(x = x + 1.5, y = y + h - 7.5, label = body), hjust = 0, vjust = 1, size = 2.45, colour = "#222222", lineheight = 1.0) +
    geom_segment(data = arrows, aes(x = x0, y = y0, xend = x1, yend = y1), colour = "#1f3b57", linewidth = 0.4,
                 arrow = arrow(length = unit(0.14, "cm"), type = "closed")) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE)
  save_fig(p, "T02_workflow", width = 7.2, height = 5.6)
}

# ==============================================================================
# T03  Night selection
figT03_night_selection <- function() {
  # Purpose: show every night of both deployments, the ERA5 criteria, and which
  # nights were retained and why.  All three panels of a district share one time
  # axis; the stretch at the end of the Honolulu record when only one or two
  # loggers were still recording (no network median) is shaded and labelled.
  inv <- read_csv(file.path(PROC, "night_inventory.csv"), col_types = cols(night_date = col_character(), .default = col_guess())) %>%
    mutate(date = as.Date(night_date), mid = as.POSIXct(date, tz = "UTC") + hours(24),
           cat = case_when(selected ~ "sel", decision != "selected" & frac_calm_clear >= 0.75 ~ "exc", TRUE ~ "oth"))
  hh <- read_csv(file.path(RAW, "logger_readings.csv"), col_types = cols(datetime = col_character(), .default = col_guess())) %>%
    mutate(time_bin = round_date(read_stamps(datetime), "30 minutes"))
  hn_all <- hh %>% group_by(region, time_bin) %>%
    summarise(med = median(temp_c), lo = quantile(temp_c, 0.1), hi = quantile(temp_c, 0.9), n = n_distinct(sensor_id), .groups = "drop")
  hn <- hn_all %>% filter(n >= MIN_SENSORS_PER_NIGHT)
  catcol <- c(sel = "#f2c14e", exc = "#c0392b", oth = "#b8c4cf"); catlab <- c(sel = "Selected calm, clear night", exc = "Calm & clear but excluded (partial night / < 20 sensors)", oth = "Other night")
  panels <- list()
  for (j in seq_along(REGIONS)) {
    region <- REGIONS[j]; e <- inv %>% filter(region == !!region); h <- hn %>% filter(region == !!region); ha <- hn_all %>% filter(region == !!region)
    xl <- c(min(ha$time_bin, e$mid - hours(12)), max(ha$time_bin, e$mid + hours(12)))
    xs <- scale_x_datetime(limits = xl, expand = c(0, 0), date_breaks = "1 week", date_labels = "%e %b")
    # loggers still recording after the network fell below the completeness threshold
    last_ok <- max(h$time_bin); tail <- ha %>% filter(time_bin > last_ok)
    sparse <- if (nrow(tail) && difftime(max(tail$time_bin), last_ok, units = "hours") > 24)
      tibble(x0 = last_ok, x1 = max(tail$time_bin), nmax = max(tail$n[tail$time_bin > last_ok + hours(24)])) else NULL
    spans <- e %>% filter(cat != "oth") %>% transmute(x0 = as.POSIXct(date, tz = "UTC") + hours(18), x1 = as.POSIXct(date, tz = "UTC") + hours(30), cat)
    pt <- ggplot() + theme_thesis() +
      { if (!is.null(sparse)) geom_rect(data = sparse, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf), fill = "#e3e3e3") } +
      geom_rect(data = spans, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf, fill = cat), alpha = 0.45) +
      geom_ribbon(data = h, aes(time_bin, ymin = lo, ymax = hi), fill = "#c9d6e3") +
      geom_line(data = h, aes(time_bin, med), colour = "#1f3b57", linewidth = 0.3) +
      { if (!is.null(sparse)) annotate("text", x = sparse$x0 + (sparse$x1 - sparse$x0) / 2, y = mean(range(c(h$lo, h$hi))),
                                       label = sprintf("only %s loggers\nstill recording:\nno network median", if (sparse$nmax > 1) sprintf("1–%d", sparse$nmax) else "1"),
                                       size = 2.1, colour = "#555555", lineheight = 0.9) } +
      scale_fill_manual(values = catcol, labels = catlab, name = NULL, guide = "none") + xs +
      labs(x = NULL, y = "Air temperature (°C)", title = sprintf("(%s) %s", letters[j], REG_LABEL[region])) +
      theme(axis.text.x = element_blank())
    if (j == 1) pt <- pt + annotate("text", x = min(h$time_bin), y = Inf, hjust = 0, vjust = 1.3, size = 2.1, colour = "#1f3b57",
                                    label = "line: network median   band: 10–90th percentile of sensors")
    bars <- lapply(list(c("wind_night", "Night-mean 10m wind\n(km h⁻¹)", 10), c("cloud_night", "Night-mean cloud\ncover (%)", 25)), function(v) {
      ggplot(e, aes(mid, .data[[v[1]]], fill = cat)) + theme_thesis() + geom_col(width = 0.8 * 86400) +
        geom_hline(yintercept = as.numeric(v[3]), colour = "#c0392b", linewidth = 0.3, linetype = "dashed") +
        scale_fill_manual(values = catcol, labels = catlab, name = NULL) + xs +
        labs(x = NULL, y = v[2]) + theme(axis.title.y = element_text(size = 7), legend.position = "none")
    })
    bars[[2]] <- bars[[2]] + labs(x = "Date (2024–2025)") + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6))
    bars[[1]] <- bars[[1]] + theme(axis.text.x = element_blank())
    panels[[region]] <- list(pt, bars[[1]], bars[[2]])
  }
  leg <- cowplot::get_legend(ggplot(inv, aes(date, wind_night, fill = cat)) + geom_col() +
                               scale_fill_manual(values = catcol, labels = catlab, name = NULL) +
                               geom_hline(aes(yintercept = 10, linetype = "Threshold"), colour = "#c0392b") +
                               scale_linetype_manual(values = c(Threshold = "dashed"), name = NULL) +
                               theme_thesis() + theme(legend.position = "bottom", legend.text = element_text(size = 6.5)) +
                               guides(fill = guide_legend(nrow = 2, order = 1), linetype = guide_legend(order = 2)))
  p <- (panels$Honolulu[[1]] | panels$Ewa[[1]]) / (panels$Honolulu[[2]] | panels$Ewa[[2]]) / (panels$Honolulu[[3]] | panels$Ewa[[3]]) /
    patchwork::wrap_elements(leg) + plot_layout(heights = c(1.15, 1, 1, 0.22))
  save_fig(p, "T03_night_selection", width = 7.2, height = 6.8)
}

# ==============================================================================
# T04  Why calm and clear
figT04_why_calm_clear <- function() {
  # Purpose: show, from this record, that the spatial contrast between sites
  # collapses on windy nights (the reason for the calm filter).
  sp <- read_csv(file.path(TAB, "night_spread_vs_weather.csv"), show_col_types = FALSE) %>%
    mutate(regionf = reg_f(region), sel = ifelse(selected, "Selected calm, clear night", "Other night"))
  rho <- RES$night_spread
  one <- function(var, lab, thr, tag) {
    labs_ <- setNames(sprintf("%s (ρ = %.2f)", REG_LABEL[REGIONS], sapply(REGIONS, function(r) rho[[r]][[paste0("rho_", var)]])), REG_LABEL[REGIONS])
    ggplot(sp, aes(.data[[var]], sd_dT)) + theme_thesis() +
      geom_vline(xintercept = thr, colour = "#c0392b", linewidth = 0.3, linetype = "dashed") +
      geom_point(aes(fill = regionf, size = sel, stroke = ifelse(selected, 0.5, 0)), shape = 21, colour = "black", alpha = 0.85) +
      scale_fill_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), labels = labs_, name = NULL) +
      scale_size_manual(values = c("Selected calm, clear night" = 2.6, "Other night" = 1.6), guide = "none") +
      labs(x = lab, y = "Spread of ΔT across the network, SD (°C)", title = tag) +
      theme(legend.position = c(0.75, 0.88)) + guides(fill = guide_legend(override.aes = list(size = 2.5)))
  }
  p <- one("wind", "Night-mean 10m wind (km h⁻¹)", 10, "(a) Wind") + annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6, size = 2.1, label = "black outline = selected calm, clear night") |
       one("cloud", "Night-mean cloud cover (%)", 25, "(b) Cloud cover")
  save_fig(p, "T04_spread_vs_weather", width = 7.2, height = 3.1)
}

# ==============================================================================
# T05  Definition of dT
figT05_dT_definition <- function() {
  # Purpose: make the response variable concrete: one shared night, every
  # sensor, the network median, and the resulting departures.
  hh <- read_csv(file.path(PROC, "halfhourly_calm_clear.csv"), col_types = cols(datetime = col_character(), time_bin = col_character(), night_date = col_character(), .default = col_guess()))
  night <- "2024-12-21"
  sm <- read_csv(file.path(TAB, "site_summary.csv"), show_col_types = FALSE)
  rows <- list()
  for (j in seq_along(REGIONS)) {
    region <- REGIONS[j]
    g <- hh %>% filter(region == !!region, night_date == night) %>% mutate(h = (hour_dec - 6) %% 24) %>% arrange(sensor_id, h)
    s <- sm %>% filter(region == !!region) %>% arrange(dT_site); cool <- s$sensor_id[1]; warm <- s$sensor_id[nrow(s)]
    med <- g %>% group_by(h) %>% summarise(m = first(network_median_temp))
    g <- g %>% mutate(kind = case_when(sensor_id == warm ~ "warm", sensor_id == cool ~ "cool", TRUE ~ "other"))
    cols <- c(other = "#b8c4cf", warm = "#c0392b", cool = "#2471a3", median = "black")
    labs_ <- c(median = "Network median (reference)", warm = sprintf("%s (warmest site)", warm), cool = sprintf("%s (coolest site)", cool))
    base <- function(yvar, ylab) {
      ggplot() + theme_thesis() + annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "#f2f2f2") +
        geom_line(data = g %>% filter(kind == "other"), aes(h, .data[[yvar]], group = sensor_id), colour = "#b8c4cf", linewidth = 0.25) +
        geom_line(data = g %>% filter(kind != "other"), aes(h, .data[[yvar]], group = sensor_id, colour = kind), linewidth = 0.6) +
        scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("06:00", "12:00", "18:00", "00:00", "06:00")) +
        labs(x = NULL, y = ylab)
    }
    pT <- base("temp_c", "Air temperature (°C)") + geom_line(data = med, aes(h, m, colour = "median"), linewidth = 0.8) +
      scale_colour_manual(values = cols, labels = labs_, breaks = c("median", "warm", "cool"), name = NULL) +
      labs(title = sprintf("(%s) %s, %s", letters[j], REG_LABEL[region], format(as.Date(night), "%e %b %Y") %>% trimws())) +
      theme(legend.position = c(0.72, 0.85), legend.text = element_text(size = 6), axis.text.x = element_blank())
    pD <- base("dT_network", "ΔT = T_site − median (°C)") + geom_hline(yintercept = 0, linewidth = 0.5) +
      scale_colour_manual(values = cols, guide = "none") +
      annotate("text", x = 18, y = Inf, vjust = 1.4, size = 2.2, colour = "#555555", label = "night window\n18:00–06:00") +
      labs(x = "Hours since 06:00 HST", title = sprintf("(%s) departures", letters[j + 2]))
    rows[[region]] <- list(pT, pD)
  }
  p <- (rows$Honolulu[[1]] | rows$Ewa[[1]]) / (rows$Honolulu[[2]] | rows$Ewa[[2]])
  save_fig(p, "T05_dT_definition", width = 7.2, height = 5.6)
}

# ==============================================================================
# T06  Distributions
figT06_distributions <- function() {
  # Purpose: show what the calm/clear filter does to the sample of night-time
  # readings in each district.
  hh <- read_csv(file.path(RAW, "logger_readings.csv"), col_types = cols(datetime = col_character(), .default = col_guess())) %>%
    mutate(time_bin = round_date(read_stamps(datetime), "30 minutes"), hd = hour(time_bin) + minute(time_bin) / 60,
           night_date = format(time_bin - hours(6), "%Y-%m-%d"), is_night = hd >= 18 | hd < 6) %>% filter(is_night)
  inv <- read_csv(file.path(PROC, "night_inventory.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
  panels <- lapply(seq_along(REGIONS), function(j) {
    region <- REGIONS[j]; g <- hh %>% filter(region == !!region)
    sel <- inv %>% filter(region == !!region, selected) %>% pull(night_date)
    a <- g$temp_c; b <- g$temp_c[g$night_date %in% sel]
    dd <- bind_rows(tibble(t = a, k = sprintf("All nights (n = %s)", format(length(a), big.mark = ","))),
                    tibble(t = b, k = sprintf("Calm, clear nights (n = %s)", format(length(b), big.mark = ","))))
    ggplot() + theme_thesis() +
      geom_histogram(data = dd %>% filter(grepl("^All", k)), aes(t, after_stat(density), fill = k), breaks = seq(14, 32.5, 0.5)) +
      geom_step(data = dd %>% filter(grepl("^Calm", k)) %>% count(k, t = cut(t, seq(14, 32.5, 0.5), right = FALSE)) %>%
                  mutate(x = seq(14, 32, 0.5)[as.integer(t)], y = n / sum(n) / 0.5), aes(x, y, colour = k), linewidth = 0.6, direction = "mid") +
      geom_vline(xintercept = median(a), colour = "#7a8792", linewidth = 0.3, linetype = "dotted") +
      geom_vline(xintercept = median(b), colour = REG_COL[region], linewidth = 0.3, linetype = "dotted") +
      scale_fill_manual(values = setNames("#b8c4cf", unique(dd$k)[1]), name = NULL) +
      scale_colour_manual(values = setNames(unname(REG_COL[region]), unique(dd$k)[2]), name = NULL) +
      labs(x = "Nighttime air temperature (°C)", y = if (j == 1) "Density" else NULL, title = sprintf("(%s) %s", letters[j], REG_LABEL[region])) +
      theme(legend.position = c(0.72, 0.88), legend.text = element_text(size = 6), legend.spacing.y = unit(0, "cm"))
  })
  p <- panels[[1]] | panels[[2]]
  save_fig(p, "T06_distributions", width = 7.2, height = 3.0)
}

# ==============================================================================
# T07 / T08  Predictor maps
pred_maps <- function(region, name) {
  # Purpose: show the spatial pattern of each predictor across the network.
  b <- DOM[[region]]; s <- sites[sites$region == region, ]
  pal <- list(imperv = "Greys", tree = "Greens", bldg = "Purples", water = "Blues", height = "Oranges", svf_point = "cividis", aspect = "Reds", coast_km = "YlGnBu")
  base_layers <- basemap(b, island, buildings = NULL, key = paste0(region, "_dom"), canopy = FALSE, paved = TRUE, water = TRUE)
  panels <- lapply(PRED, function(p) {
    sc <- if (p == "svf_point") scale_fill_viridis_c(option = "cividis", name = NULL) else scale_fill_distiller(palette = pal[[p]], direction = 1, name = NULL)
    ggplot() + theme_map() + base_layers +
      geom_sf(data = s, aes(fill = .data[[p]]), shape = 21, size = 1.7, colour = "black", stroke = 0.2) + sc +
      labs(title = if (p == "imperv") sprintf("Impervious surface (%%, r = %dm)", RADIUS) else PLAB[[p]]) +
      map_coord(b, graticule = FALSE) +
      theme(plot.title = element_text(size = 6.5, face = "plain"), legend.position = "bottom", legend.key.height = unit(0.15, "cm"),
            legend.key.width = unit(0.6, "cm"), legend.text = element_text(size = 5.5), legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 3, 1, 3))
  })
  p <- wrap_plots(panels, ncol = 4)
  save_fig(p, name, width = 7.2, height = 2 * (1.75 * map_aspect(b) + 0.55))
}
figT07_pred_maps_hnl <- function() pred_maps("Honolulu", "T07_predictor_maps_Honolulu")
figT08_pred_maps_ewa <- function() pred_maps("Ewa", "T08_predictor_maps_Ewa")

# ==============================================================================
# T09  Scale selection
figT09_scale <- function() {
  # Purpose: justify the source-area radius with the same models used later.
  sc <- as_tibble(RES$scale_selection) %>% mutate(scope = factor(scope, levels = c("pooled", REGIONS), labels = c("Pooled (+ district)", REG_LABEL[REGIONS])))
  ss <- read_csv(file.path(TAB, "scale_single_predictor.csv"), show_col_types = FALSE)
  pa <- ggplot(sc, aes(factor(radius_m, labels = c("50m", "100m", "200m")), delta_aic, fill = scope)) + theme_thesis() +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f", delta_aic), colour = scope), position = position_dodge(0.75), vjust = -0.4, size = 1.9, show.legend = FALSE) +
    geom_hline(yintercept = 2, colour = "#c0392b", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = 0.55, y = 2.3, label = "ΔAIC = 2", size = 2.1, colour = "#c0392b", hjust = 0) +
    scale_fill_manual(values = c("#555555", unname(REG_COL[REGIONS])), name = NULL) + scale_colour_manual(values = c("#555555", unname(REG_COL[REGIONS]))) +
    labs(x = NULL, y = "ΔAIC of full model (ML)", title = "(a) Full model by radius") + theme(legend.position = c(0.7, 0.85), legend.text = element_text(size = 6))
  pb <- lapply(seq_along(REGIONS), function(j) {
    region <- REGIONS[j]
    g <- ss %>% filter(region == !!region, !predictor %in% c("svf_point", "coast_km", "water")) %>% mutate(predictor = factor(PSHORT[predictor], levels = PSHORT))
    ggplot(g, aes(radius_m, r2m, colour = predictor)) + theme_thesis() + geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
      scale_x_continuous(breaks = c(50, 100, 200)) + scale_colour_brewer(palette = "Dark2", name = NULL) +
      labs(x = "Radius (m)", y = "Marginal R² (single predictor)", title = sprintf("(%s) %s", letters[j + 1], REG_LABEL[region])) +
      theme(legend.position = if (region == "Ewa") c(0.7, 0.78) else "none", legend.text = element_text(size = 5.5), legend.key.size = unit(0.25, "cm"))
  })
  p <- (pa | pb[[1]] | pb[[2]]) + plot_layout(widths = c(1.1, 1, 1))
  save_fig(p, "T09_scale_selection", width = 7.2, height = 2.8)
}

# ==============================================================================
# T10  Correlation and VIF
figT10_correlation <- function() {
  # Purpose: show the collinearity structure that limits the full models.
  panels <- lapply(seq_along(REGIONS), function(j) {
    region <- REGIONS[j]
    cc <- read.csv(file.path(TAB, paste0("correlation_", region, ".csv")), row.names = 1)[PRED, PRED]
    vif <- unlist(RES$collinearity[[region]]$vif)
    d <- expand.grid(i = seq_along(PRED), j = seq_along(PRED)) %>%
      mutate(r = mapply(function(a, b) cc[a, b], i, j), lab = ifelse(j < i, sprintf("%.2f", r), ifelse(j == i, sprintf("VIF\n%.1f", vif[PRED[i]]), "")),
             fillv = ifelse(j < i, r, NA), txt = ifelse(j < i & abs(r) > 0.6, "white", "black"))
    ggplot(d, aes(j, i)) + theme_thesis() +
      geom_tile(data = d %>% filter(j < i), aes(fill = fillv), colour = "white", linewidth = 0.3) +
      geom_text(aes(label = lab, colour = txt), size = 1.9, lineheight = 0.8) + scale_colour_identity() +
      scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "Pearson r\n(site level)", na.value = NA) +
      scale_x_continuous(breaks = seq_along(PRED), labels = PSHORT[PRED], expand = c(0, 0)) +
      scale_y_reverse(breaks = seq_along(PRED), labels = PSHORT[PRED], expand = c(0, 0)) + coord_equal() +
      labs(x = NULL, y = NULL, title = sprintf("(%s) %s", letters[j], REG_LABEL[region])) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6), axis.line = element_blank(),
            axis.ticks = element_blank(), legend.position = if (j == 2) "right" else "none", legend.key.height = unit(0.8, "cm"), legend.key.width = unit(0.25, "cm"))
  })
  p <- panels[[1]] | panels[[2]]
  save_fig(p, "T10_correlation", width = 7.2, height = 3.6)
}

# ==============================================================================
# T11  Bivariate
figT11_bivariate <- function() {
  # Purpose: show the raw relationships that the mixed models summarise.
  sm <- read_csv(file.path(TAB, "site_summary.csv"), show_col_types = FALSE)
  single <- as_tibble(RES$single_predictor); sds <- RES$z_sds
  long <- sm %>% pivot_longer(all_of(PRED), names_to = "predictor", values_to = "x") %>% mutate(regionf = reg_f(region))
  lines <- long %>% group_by(region, regionf, predictor) %>%
    summarise(xmin = min(x), xmax = max(x), xm = mean(x), ym = mean(dT_site), .groups = "drop") %>%
    rowwise() %>% mutate(beta = single$beta_z[single$scope == region & single$predictor == predictor] / sds[[predictor]],
                         p = single$lrt_p[single$scope == region & single$predictor == predictor]) %>% ungroup() %>%
    mutate(y0 = ym + beta * (xmin - xm), y1 = ym + beta * (xmax - xm), sig = ifelse(p < 0.05, "LMM slope, p < 0.05", "p ≥ 0.05"))
  long$predictor <- factor(long$predictor, levels = PRED, labels = PLAB[PRED]); lines$predictor <- factor(lines$predictor, levels = PRED, labels = PLAB[PRED])
  p <- ggplot(long, aes(x, dT_site)) + theme_thesis() + geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.25) +
    geom_point(aes(colour = regionf), size = 1.1, alpha = 0.8) +
    geom_segment(data = lines, aes(x = xmin, xend = xmax, y = y0, yend = y1, colour = regionf, linetype = sig), linewidth = 0.5) +
    facet_wrap(~ predictor, scales = "free_x", ncol = 4, strip.position = "bottom") +
    scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL) +
    scale_linetype_manual(values = c("LMM slope, p < 0.05" = "solid", "p ≥ 0.05" = "dotted"), name = NULL) +
    labs(x = NULL, y = "Site-mean ΔT (°C)") +
    theme(strip.placement = "outside", strip.text = element_text(size = 6.5, face = "plain"), legend.position = "bottom", axis.text = element_text(size = 6))
  save_fig(p, "T11_bivariate", width = 7.2, height = 4.6)
}

# ==============================================================================
# T12  Coefficients
figT12_coefficients <- function() {
  # Purpose: compare the standardised effects of all predictors in the pooled
  # and the two district full models, with the best-subset models overlaid.
  sets <- list(list("(a) Pooled + district", RES$full_pooled$coefs, RES$parsimonious_pooled$coefs, "#555555"),
               list("(b) Honolulu", RES$regional$Honolulu$full$coefs, RES$regional$Honolulu$best_subset$coefs, REG_COL[["Honolulu"]]),
               list("(c) ʻEwa", RES$regional$Ewa$full$coefs, RES$regional$Ewa$best_subset$coefs, REG_COL[["Ewa"]]))
  d <- bind_rows(lapply(sets, function(s) bind_rows(
    as_tibble(s[[2]]) %>% filter(term %in% paste0(PRED, "_z")) %>% mutate(model = "full", off = 0),
    as_tibble(s[[3]]) %>% filter(term %in% paste0(PRED, "_z")) %>% mutate(model = "best", off = -0.28)) %>% mutate(panel = s[[1]], col = s[[4]])))
  d <- d %>% mutate(pred = factor(sub("_z$", "", term), levels = rev(PRED), labels = rev(PSHORT[PRED])), y = as.numeric(pred) + off,
                    fillc = ifelse(model == "full" & p >= 0.05, "white", col), shp = ifelse(model == "full", 21, 22), sz = ifelse(model == "full", 1.8, 1.3),
                    al = ifelse(model == "full", 1, 0.6), panel = factor(panel, levels = sapply(sets, `[[`, 1)))
  p <- ggplot(d, aes(estimate, y)) + theme_thesis() + geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = col, alpha = al), height = 0.25, linewidth = 0.4) +
    geom_point(aes(colour = col, fill = fillc, shape = shp, size = sz, alpha = al), stroke = 0.6) +
    scale_colour_identity() + scale_fill_identity() + scale_shape_identity() + scale_size_identity() + scale_alpha_identity() +
    scale_y_continuous(breaks = seq_along(PRED), labels = rev(PSHORT[PRED])) + coord_cartesian(xlim = c(-1.3, 1.3)) +
    facet_wrap(~ panel) + labs(x = "Effect on ΔT per 1 SD (°C), 95% CI", y = NULL) +
    theme(axis.text.y = element_text(size = 7))
  leg <- ggplot(tibble(k = factor(c("Full model (filled: p < 0.05)", "Full model, p ≥ 0.05", "Parsimonious / best-subset model"),
                                  levels = c("Full model (filled: p < 0.05)", "Full model, p ≥ 0.05", "Parsimonious / best-subset model")), x = 1:3), aes(x, 1)) +
    geom_point(aes(shape = k, fill = k), size = 2, colour = "#555555") +
    scale_shape_manual(values = c(21, 21, 22), name = NULL) + scale_fill_manual(values = c("#555555", "white", "#999999"), name = NULL) +
    theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 6.5))
  p <- p / patchwork::wrap_elements(cowplot::get_legend(leg)) + plot_layout(heights = c(1, 0.08))
  save_fig(p, "T12_coefficients", width = 7.2, height = 3.5)
}

# ==============================================================================
# T13  Interactions
figT13_interactions <- function() {
  # Purpose: show the two drivers whose effect differs between districts.
  sm <- read_csv(file.path(TAB, "site_summary.csv"), show_col_types = FALSE) %>% mutate(regionf = reg_f(region))
  inter <- as_tibble(RES$interactions$single); sds <- RES$z_sds
  panels <- lapply(seq_along(c("height", "coast_km")), function(k) {
    p <- c("height", "coast_km")[k]
    lines <- sm %>% group_by(region, regionf) %>% summarise(xmin = min(.data[[p]]), xmax = max(.data[[p]]), xm = mean(.data[[p]]), ym = mean(dT_site), .groups = "drop") %>%
      mutate(beta = sapply(region, function(r) inter[[if (r == "Honolulu") "hnl_slope_z" else "ewa_slope_z"]][inter$predictor == p]) / sds[[p]],
             y0 = ym + beta * (xmin - xm), y1 = ym + beta * (xmax - xm))
    ggplot(sm, aes(.data[[p]], dT_site, colour = regionf)) + theme_thesis() + geom_hline(yintercept = 0, colour = "#999999", linewidth = 0.25) +
      geom_point(size = 1.2, alpha = 0.75) + geom_segment(data = lines, aes(x = xmin, xend = xmax, y = y0, yend = y1), linewidth = 0.7) +
      scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL) +
      labs(x = PLAB[[p]], y = "Site-mean ΔT (°C)", title = sprintf("(%s) %s: interaction p = %.3f", letters[k], PSHORT[[p]], inter$p[inter$predictor == p])) +
      theme(legend.position = c(0.85, 0.9))
  })
  save_fig(panels[[1]] | panels[[2]], "T13_interactions", width = 7.2, height = 3.0)
}

# ==============================================================================
# T14  Threshold sensitivity
figT14_thresholds <- function() {
  # Purpose: show that the conclusions do not hinge on the exact calm/clear
  # thresholds, and that filtering sharpens the morphological signal.
  t <- as_tibble(RES$threshold_sensitivity)
  nofilter <- t %>% filter(wind_max == 999); comp <- t %>% filter(wind_max == -1)
  g <- t %>% filter(min_frac == 0.75, wind_max > 0) %>% mutate(cloud = factor(cloud_max, labels = sprintf("cloud < %d%%", sort(unique(cloud_max)))))
  cloud_cols <- c("#08306b", "#2171b5", "#6baed6", "#bdd7e7")
  one <- function(v, se, ylab, title, legend = FALSE) {
    p <- ggplot(g, aes(wind_max, .data[[v]], colour = cloud)) + theme_thesis() +
      geom_hline(yintercept = nofilter[[v]], colour = "#c0392b", linewidth = 0.4, linetype = "dashed") +
      geom_hline(yintercept = comp[[v]], colour = "#7f8c8d", linewidth = 0.4, linetype = "dotted") +
      geom_vline(xintercept = 10, colour = "#999999", linewidth = 0.25) + geom_line(linewidth = 0.4) + geom_point(size = 1.2)
    if (!is.null(se)) p <- p + geom_errorbar(aes(ymin = .data[[v]] - .data[[se]], ymax = .data[[v]] + .data[[se]]), width = 0.5, linewidth = 0.3)
    p + scale_colour_manual(values = cloud_cols, name = "Cloud threshold") + scale_x_continuous(breaks = c(8, 10, 12, 15, 20)) +
      labs(x = "Wind threshold (km h⁻¹)", y = ylab, title = title) +
      theme(legend.position = if (legend) c(0.78, 0.82) else "none", legend.text = element_text(size = 5.5), legend.title = element_text(size = 5.5), legend.key.size = unit(0.25, "cm"))
  }
  pa <- one("beta_imperv_z", "se_imperv", "Impervious effect per SD (°C)", "(a) Impervious surface")
  pb <- one("beta_height_z", "se_height", "Height effect per SD (°C)", "(b) Building height")
  pc <- one("r2m_full", NULL, "Marginal R², full model", "(c) Explained variance", legend = TRUE)
  leg <- cowplot::get_legend(ggplot(tibble(x = 1:3, k = factor(c(sprintf("No weather filter (%d sensor-nights)", nofilter$sensor_nights),
                                                                  sprintf("Windy or cloudy nights only (%d)", comp$sensor_nights), "Adopted threshold (10km h⁻¹, 25%)"),
                                                                levels = c(sprintf("No weather filter (%d sensor-nights)", nofilter$sensor_nights), sprintf("Windy or cloudy nights only (%d)", comp$sensor_nights), "Adopted threshold (10km h⁻¹, 25%)"))),
                                    aes(x, x, linetype = k, colour = k)) + geom_line() +
                               scale_linetype_manual(values = c("dashed", "dotted", "solid"), name = NULL) + scale_colour_manual(values = c("#c0392b", "#7f8c8d", "#999999"), name = NULL) +
                               theme_thesis() + theme(legend.position = "bottom", legend.text = element_text(size = 6.3), legend.key.width = unit(0.7, "cm"), legend.spacing.x = unit(0.4, "cm")))
  p <- ((pa | pb | pc) / patchwork::wrap_elements(leg)) + plot_layout(heights = c(1, 0.08))
  save_fig(p, "T14_threshold_sensitivity", width = 7.2, height = 3.1)
}

# ==============================================================================
# T15  BLUPs
figT15_blups <- function() {
  # Purpose: show where the full model under- or over-predicts, i.e. the heat
  # that morphology does not explain.
  b <- read_csv(file.path(TAB, "blups.csv"), show_col_types = FALSE) %>% arrange(blup) %>% mutate(sensor_id = factor(sensor_id, levels = sensor_id), regionf = reg_f(region))
  pa <- ggplot(b, aes(sensor_id, blup, fill = regionf)) + theme_thesis() + geom_col(width = 0.8) +
    { if (all(!is.na(b$blup_sd))) geom_errorbar(aes(ymin = blup - 1.96 * blup_sd, ymax = blup + 1.96 * blup_sd), width = 0, linewidth = 0.25, colour = "#555555") } +
    geom_hline(yintercept = 0, linewidth = 0.3) + scale_fill_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL) +
    labs(x = NULL, y = "Random intercept, BLUP (°C)", title = "(a) Site random intercepts of the full pooled model (± 95% conditional interval)") +
    theme(axis.text.x = element_text(angle = 90, size = 4.8, vjust = 0.5, hjust = 1), legend.position = c(0.08, 0.85))
  bs <- st_transform(st_as_sf(b, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)
  maps <- lapply(seq_along(REGIONS), function(k) {
    region <- REGIONS[k]; bb <- EQB[[region]]
    ggplot() + theme_map() + basemap(bb, island, buildings = bld(region), key = region) +
      geom_sf(data = bs[bs$region == region, ], aes(fill = blup), shape = 21, size = 3, colour = "black", stroke = 0.3) +
      scale_fill_distiller(palette = "PuOr", limits = c(-1, 1), oob = squish, name = "Unexplained site effect, BLUP (°C): warmer (+) or cooler (−) than the morphology predicts") +
      place_labels(region, bb, size = 1.9) + scale_bar(bb, 1000) + north_arrow(bb, loc = c(0.94, 0.86)) +
      panel_tag(bb, sprintf("(%s) %s", letters[k + 1], REG_LABEL[region])) + map_coord(bb, graticule = FALSE) +
      guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
      theme(legend.position = "bottom", legend.title = element_text(size = 6.5), legend.key.width = unit(1.2, "cm"),
            legend.key.height = unit(0.15, "cm"), legend.text = element_text(size = 6), plot.margin = margin(2, 2, 2, 2))
  })
  hmap <- 3.45 * map_aspect(EQB$Honolulu)
  p <- pa / ((maps[[1]] | maps[[2]]) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_layout(heights = c(2.6, hmap + 0.9))
  save_fig(p, "T15_blups", width = 7.2, height = hmap + 3.6)
}

# ==============================================================================
# T16  dT maps
figT16_dT_maps <- function() {
  # Purpose: the observed pattern itself: site-mean dT on the urban-form base map.
  maps <- lapply(seq_along(REGIONS), function(k) {
    region <- REGIONS[k]; b <- EQB[[region]]
    ggplot() + theme_map() + basemap(b, island, buildings = bld(region), key = region) +
      geom_sf(data = sites[sites$region == region, ], aes(fill = dT_site), shape = 21, size = 3.2, colour = "black", stroke = 0.3) +
      scale_fill_distiller(palette = "RdBu", limits = c(-1.6, 1.6), oob = squish, name = "Site-mean nocturnal departure from the network median, ΔT (°C)") +
      place_labels(region, b, size = 1.9) + scale_bar(b, 1000) + north_arrow(b, loc = c(0.94, 0.86)) +
      panel_tag(b, sprintf("(%s) %s", letters[k], REG_LABEL[region])) + graticule_breaks() + map_coord(b, graticule = TRUE) +
      guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
      theme(legend.position = "bottom", legend.title = element_text(size = 6.5), legend.key.width = unit(1.2, "cm"),
            legend.key.height = unit(0.15, "cm"), legend.text = element_text(size = 6), axis.text = element_text(size = 5.5),
            axis.text.y = element_text(angle = 90, hjust = 0.5), plot.margin = margin(2, 2, 2, 2))
  })
  leg <- cowplot::get_legend(ggplot(tibble(k = factor(c("Buildings", "Paved", "Tree canopy", "Water"), levels = c("Buildings", "Paved", "Tree canopy", "Water")), x = 1:4), aes(x, 1, fill = k)) +
                               geom_point(shape = 22, size = 3, colour = NA) + scale_fill_manual(values = c(C_BLDG, C_PAVED, C_TREE, C_WATER), name = NULL) +
                               theme_void() + theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.key.size = unit(0.3, "cm")) + guides(fill = guide_legend(nrow = 1)))
  hmap <- 3.45 * map_aspect(EQB$Honolulu)
  p <- ((maps[[1]] | maps[[2]]) + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / patchwork::wrap_elements(leg) + plot_layout(heights = c(hmap + 0.9, 0.3))
  save_fig(p, "T16_dT_maps", width = 7.2, height = hmap + 1.35)
}

# ==============================================================================
# T17  Shared nights
figT17_shared_nights <- function() {
  # Purpose: compare the two networks in absolute terms on the nights they share.
  rows <- as_tibble(RES$shared_nights$rows); nights <- RES$shared_nights$nights
  hh <- read_csv(file.path(PROC, "halfhourly_calm_clear.csv"), col_types = cols(datetime = col_character(), time_bin = col_character(), night_date = col_character(), .default = col_guess()))
  panels <- lapply(seq_along(nights), function(k) {
    nd <- nights[k]
    g <- hh %>% filter(night_date == nd) %>% mutate(h = (hour_dec - 6) %% 24)
    q <- g %>% group_by(region, h) %>% summarise(lo = quantile(temp_c, 0.1), hi = quantile(temp_c, 0.9), med = first(network_median_temp), .groups = "drop") %>%
      mutate(regionf = reg_f(region))
    labs_ <- setNames(sapply(REGIONS, function(r) sprintf("%s median (night mean %.1f°C)", REG_LABEL[r], rows$T_median_night[rows$night_date == nd & rows$region == r])), REG_LABEL[REGIONS])
    ggplot(q, aes(h)) + theme_thesis() + annotate("rect", xmin = 12, xmax = 24, ymin = -Inf, ymax = Inf, fill = "#f2f2f2") +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = regionf), alpha = 0.18) + geom_line(aes(y = med, colour = regionf), linewidth = 0.7) +
      scale_fill_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), guide = "none") +
      scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), labels = labs_, name = NULL) +
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("06:00", "12:00", "18:00", "00:00", "06:00")) +
      labs(x = "Time (HST)", y = if (k == 1) "Air temperature (°C)" else NULL, title = sprintf("(%s) %s", letters[k], trimws(format(as.Date(nd), "%e %b %Y")))) +
      theme(legend.position = c(0.62, 0.9), legend.text = element_text(size = 6))
  })
  save_fig(wrap_plots(panels, nrow = 1), "T17_shared_nights", width = 7.2, height = 3.0)
}

# ==============================================================================
# T18  LOSO
figT18_loso <- function() {
  # Purpose: report honest out-of-sample skill of the pooled and the district models.
  lp <- read_csv(file.path(TAB, "loso_predictions.csv"), show_col_types = FALSE) %>% mutate(regionf = reg_f(region))
  cv <- as_tibble(RES$loso_cv)
  r <- cv %>% filter(model == "full")
  pa <- ggplot(lp, aes(pred_full, dT_site, colour = regionf)) + theme_thesis() + geom_abline(colour = "#999999", linewidth = 0.3) + geom_point(size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), name = NULL) + coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    labs(x = "Predicted site-mean ΔT (°C), site held out", y = "Observed site-mean ΔT (°C)", title = sprintf("(a) Pooled full model: r = %.2f, RMSE = %.2f°C", r$r, r$rmse)) +
    theme(legend.position = c(0.15, 0.9))
  labs_ <- setNames(sapply(REGIONS, function(reg) {
    bs <- gsub("coast_km", "coast", gsub("imperv", "impervious", RES$regional[[reg]]$best_subset$predictors))
    row <- cv %>% filter(startsWith(model, reg), grepl("best subset", model))
    sprintf("%s: %s (r = %.2f, RMSE = %.2f°C)", REG_LABEL[reg], bs, row$r, row$rmse) }), REG_LABEL[REGIONS])
  pb <- ggplot(lp, aes(pred_regional_best, dT_site, colour = regionf)) + theme_thesis() + geom_abline(colour = "#999999", linewidth = 0.3) + geom_point(size = 1.2, alpha = 0.8) +
    scale_colour_manual(values = setNames(REG_COL[REGIONS], REG_LABEL[REGIONS]), labels = labs_, name = NULL) + coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
    labs(x = "Predicted site-mean ΔT (°C), site held out", y = NULL, title = "(b) District best-subset models") +
    theme(legend.position = c(0.4, 0.92), legend.text = element_text(size = 5.8))
  save_fig(pa | pb, "T18_loso", width = 7.2, height = 3.2)
}

# ==============================================================================
# T19  Night-to-night stability
figT19_night_stability <- function() {
  # Purpose: show how repeatable the spatial pattern is from night to night.
  d <- read_csv(file.path(TAB, "analysis_dataset.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
  panels <- lapply(seq_along(REGIONS), function(j) {
    region <- REGIONS[j]; g <- d %>% filter(region == !!region)
    w <- g %>% select(sensor_id, night_id, dT_night) %>% pivot_wider(names_from = night_id, values_from = dT_night) %>% select(-sensor_id)
    cc <- cor(w, use = "pairwise.complete.obs")
    nights <- g %>% distinct(night_id, night_date, .keep_all = FALSE) %>% arrange(night_date)
    wind <- g %>% group_by(night_id) %>% summarise(w = mean(wind_night))
    lab <- setNames(trimws(format(as.Date(nights$night_date), "%e %b")), nights$night_id)
    dd <- expand.grid(a = colnames(cc), b = colnames(cc), stringsAsFactors = FALSE) %>% mutate(r = mapply(function(x, y) cc[x, y], a, b))
    ggplot(dd, aes(factor(a, levels = nights$night_id), factor(b, levels = rev(nights$night_id)), fill = r)) + theme_thesis() +
      geom_tile(colour = "white", linewidth = 0.3) + geom_text(aes(label = sprintf("%.2f", r), colour = ifelse(r < 0.6, "white", "black")), size = 1.9) + scale_colour_identity() +
      scale_fill_viridis_c(limits = c(-0.2, 1), name = "Correlation of site ΔT\nbetween nights") +
      scale_x_discrete(labels = sprintf("%s\n%.1f", lab[nights$night_id], wind$w[match(nights$night_id, wind$night_id)])) +
      scale_y_discrete(labels = rev(lab[nights$night_id])) + coord_equal() +
      labs(x = "Night (date, mean wind km h⁻¹)", y = NULL, title = sprintf("(%s) %s", letters[j], REG_LABEL[region])) +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 5.5), axis.text.y = element_text(size = 6),
            legend.position = if (j == 2) "right" else "none", legend.key.height = unit(0.7, "cm"), legend.key.width = unit(0.25, "cm"))
  })
  p <- (panels[[1]] | panels[[2]]) + plot_layout(widths = c(7, 4.5))
  save_fig(p, "T19_night_stability", width = 7.2, height = 3.3)
}

FUNCS <- list("01" = figT01_study_area, "02" = figT02_workflow, "03" = figT03_night_selection, "04" = figT04_why_calm_clear,
              "05" = figT05_dT_definition, "06" = figT06_distributions, "07" = figT07_pred_maps_hnl, "08" = figT08_pred_maps_ewa,
              "09" = figT09_scale, "10" = figT10_correlation, "11" = figT11_bivariate, "12" = figT12_coefficients,
              "13" = figT13_interactions, "14" = figT14_thresholds, "15" = figT15_blups, "16" = figT16_dT_maps,
              "17" = figT17_shared_nights, "18" = figT18_loso, "19" = figT19_night_stability)

if (sys.nframe() == 0) {
  which <- commandArgs(trailingOnly = TRUE); if (!length(which)) which <- names(FUNCS)
  for (k in which) FUNCS[[k]]()
}
