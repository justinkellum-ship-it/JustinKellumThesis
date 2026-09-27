#!/usr/bin/env Rscript
# 08_appendix_figures.R ---------------------------------------------------------
# The three maps and diagrams of Appendix B of the thesis.
#
#   B1  land cover of Oʻahu and of the two study domains (impervious surface and
#       tree canopy from the 1 m NOAA C-CAP masks), with the areas and shares
#       quoted in Section 1.2 (results/island_landcover.json, written by
#       T01 in 03_figures.R from the same aggregated masks)
#   B2  terrain: Oʻahu and the Honolulu district from the USGS 10 m DEM
#       (subsamples retrieved from the PacIOOS ERDDAP server and stored as text in
#       data/external/oahu_dem_500m_q5.txt and honolulu_dem_150m_q5.txt)
#   B3  how the sky view factor is calculated, for an example site (the same ray
#       casting as sky_view_factor() in 01_site_predictors.R, with the angle of
#       every ray kept)
#
# Usage: Rscript 08_appendix_figures.R [B1 B2 B3]   (no argument = all three)
#
# Packages: ggplot2, patchwork, sf, terra (as 03_figures.R).
# ------------------------------------------------------------------------------
.here <- if (nzchar(Sys.getenv("THESIS_ROOT"))) file.path(Sys.getenv("THESIS_ROOT"), "analysis") else {
  .f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)); if (length(.f)) dirname(.f[1]) else "." }
source(file.path(.here, "helpers.R"), encoding = "UTF-8"); source(file.path(.here, "maplib.R"), encoding = "UTF-8")

island <- load_island()
sites  <- load_sites(100)
DOM <- list(Honolulu = domain_bounds(sites, "Honolulu", 600), Ewa = domain_bounds(sites, "Ewa", 600))
EQB <- equal_bounds(DOM)                                   # the study domains of Figure 1 and Section 1.2
LC  <- fromJSON(file.path(OUT, "island_landcover.json"))
fmt_km2 <- function(x) ifelse(x >= 100, format(round(x), big.mark = ","), sprintf("%.1f", x))

# ==============================================================================
# B1  Land cover
figB1_landcover <- function() {
  # Purpose: document the island and domain land-cover figures of Section 1.2.
  ib <- st_bbox(island); bT01 <- bounds4(ib$xmin - 1500, ib$ymin - 1500, ib$xmax + 9000, ib$ymax + 1500)   # box of the cached masks (T01)
  bI <- bounds4(ib$xmin - 1500, ib$ymin - 1500, ib$xmax + 1500, ib$ymax + 1500)
  isl <- list(impervious = tint_raster("impervious", 1, bT01, 20, "island"), canopy = tint_raster("tree", 1, bT01, 20, "island"))
  isl <- lapply(isl, function(r) aggregate(r, fact = 3, fun = "mean"))                # 60 m cells for drawing
  boxes <- bind_rows(lapply(names(EQB), function(r) tibble(region = r, xmin = EQB[[r]]["xmin"], xmax = EQB[[r]]["xmax"],
                                                           ymin = EQB[[r]]["ymin"], ymax = EQB[[r]]["ymax"])))
  col <- c(impervious = "#3a3a3a", canopy = "#2e7d32")
  panel <- function(r, kind, bounds, title, show_boxes = FALSE, sbar = 10000) {
    d <- raster_df(r) %>% filter(v > 0.02)
    land <- st_intersection(st_geometry(island), st_as_sfc(st_bbox(bounds, crs = st_crs(CRS_M))))
    ggplot() + theme_map() + theme(panel.background = element_rect(fill = C_OCEAN, colour = NA), legend.position = "none",
                                   plot.title = element_text(size = 7.5, face = "bold")) +
      geom_sf(data = land, fill = "white", colour = NA) +
      geom_raster(data = d, aes(x, y, alpha = v), fill = col[[kind]], interpolate = TRUE) + scale_alpha_identity() +
      geom_sf(data = st_boundary(st_geometry(island)), colour = C_COAST, linewidth = 0.25) +
      { if (show_boxes) geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = region), fill = NA, linewidth = 0.6) } +
      { if (show_boxes) scale_colour_manual(values = REG_COL, guide = "none") } +
      { if (!show_boxes) geom_sf(data = sites[sites$region == names(which(sapply(EQB, identical, bounds))), ], colour = "black", size = 0.35) } +
      scale_bar(bounds, sbar, loc = c(0.05, 0.05)) + labs(title = title) + map_coord(bounds, graticule = FALSE)
  }
  dom <- lapply(names(EQB), function(r) list(
    impervious = aggregate(tint_raster("impervious", 1, EQB[[r]], 4, r), fact = 3, fun = "mean"),
    canopy = aggregate(tint_raster("tree", 1, EQB[[r]], 4, r), fact = 3, fun = "mean")))
  names(dom) <- names(EQB)
  can_km2 <- function(r) LC[[paste0(r, "_canopy_pct_of_land")]] * LC[[paste0(r, "_domain_land_km2")]] / 100
  pA <- panel(isl$impervious, "impervious", bI, sprintf("(a) Oʻahu, impervious surface: %skm², %.1f%% of the land", fmt_km2(LC$impervious_km2), LC$impervious_pct), TRUE)
  pB <- panel(isl$canopy, "canopy", bI, sprintf("(b) Oʻahu, tree canopy: %skm², %.1f%% of the land", fmt_km2(LC$canopy_km2), LC$canopy_pct), TRUE)
  pC <- panel(dom$Honolulu$impervious, "impervious", EQB$Honolulu,
              sprintf("(c) Honolulu domain, impervious: %skm², %.1f%%", fmt_km2(LC$Honolulu_impervious_km2), LC$Honolulu_impervious_pct_of_land), sbar = 1000)
  pD <- panel(dom$Honolulu$canopy, "canopy", EQB$Honolulu,
              sprintf("(d) Honolulu domain, tree canopy: %skm², %.1f%%", fmt_km2(can_km2("Honolulu")), LC$Honolulu_canopy_pct_of_land), sbar = 1000)
  pE <- panel(dom$Ewa$impervious, "impervious", EQB$Ewa,
              sprintf("(e) ʻEwa domain, impervious: %skm², %.1f%%", fmt_km2(LC$Ewa_impervious_km2), LC$Ewa_impervious_pct_of_land), sbar = 1000)
  pF <- panel(dom$Ewa$canopy, "canopy", EQB$Ewa,
              sprintf("(f) ʻEwa domain, tree canopy: %skm², %.1f%%", fmt_km2(can_km2("Ewa")), LC$Ewa_canopy_pct_of_land), sbar = 1000)
  hI <- map_aspect(bI) * 3.55; hD <- map_aspect(EQB$Honolulu) * 3.55
  p <- (pA | pB) / (pC | pD) / (pE | pF) + plot_layout(heights = c(hI, hD, hD))
  save_fig(p, "B1_landcover", width = 7.2, height = hI + 2 * hD + 0.75)
}

# ==============================================================================
# B2  Terrain
read_dem_txt <- function(path) {
  # the compact text grids of data/derived: header with the grid geometry, then one line per
  # row (north to south) with elevation / 5 m, "wN" for N cells without data and a checksum
  L <- readLines(path, encoding = "UTF-8")
  hdr <- paste(L[startsWith(L, "#")], collapse = " ")
  num <- function(pat) as.numeric(sub(pat, "\\1", regmatches(hdr, regexpr(pat, hdr))))
  nrow_ <- num("Grid: ([0-9]+) rows"); ncol_ <- num("x ([0-9]+)\\s*#?\\s*columns")
  lat0 <- num("latitude ([0-9.]+) N"); dlat <- num("descending by ([0-9.]+) deg")
  lon0 <- num("longitude (-[0-9.]+) E"); dlon <- num("ascending by ([0-9.]+) deg")
  rows <- L[!startsWith(L, "#") & nzchar(L)]
  M <- t(sapply(rows, function(line) {
    parts <- strsplit(line, " # ", fixed = TRUE)[[1]]; chk <- as.numeric(strsplit(parts[2], " ")[[1]])
    tok <- strsplit(parts[1], " ")[[1]][-1]
    v <- unlist(lapply(tok, function(t) if (startsWith(t, "w")) rep(NA_real_, as.integer(substring(t, 2))) else as.numeric(t)))
    stopifnot(length(v) == ncol_, sum(!is.na(v)) == chk[1], sum(v, na.rm = TRUE) == chk[2])      # checksum of the row
    v * 5
  }, USE.NAMES = FALSE))
  stopifnot(nrow(M) == nrow_)
  r <- rast(nrows = nrow_, ncols = ncol_, xmin = lon0 - dlon / 2, xmax = lon0 + (ncol_ - 0.5) * dlon,
            ymin = lat0 - (nrow_ - 0.5) * dlat, ymax = lat0 + dlat / 2, crs = "EPSG:4326")
  values(r) <- as.vector(t(M))
  r
}
hypso <- function(z, zmax) {
  # hypsometric tint, multiplied by a hillshade later
  stops <- c(0, 50, 150, 300, 500, 750, 1000, 1250)
  cols <- c("#dfe9cf", "#cfe0b4", "#e5dcaa", "#dcc594", "#c9a77c", "#ad8a67", "#8f7157", "#7a6352")
  ramp <- colorRamp(cols, space = "Lab")
  ramp(approx(stops, seq(0, 1, length.out = length(stops)), xout = pmin(pmax(z, 0), 1250))$y)
}
dem_layers <- function(r_ll, res_m, breaks) {
  r <- project(r_ll, paste0("EPSG:", CRS_M), method = "bilinear", res = res_m)
  sl <- terrain(r, "slope", unit = "radians"); as_ <- terrain(r, "aspect", unit = "radians")
  hs <- shade(sl, as_, angle = 40, direction = 315)
  full <- as.data.frame(c(r, hs), xy = TRUE, na.rm = FALSE); names(full) <- c("x", "y", "z", "hs")
  d <- full %>% filter(!is.na(z))
  hsn <- ifelse(is.na(d$hs), 1, 0.55 + 0.45 * pmin(pmax(d$hs, 0), 1))
  rgbm <- hypso(d$z) * hsn
  d$col <- rgb(rgbm[, 1], rgbm[, 2], rgbm[, 3], maxColorValue = 255)
  list(r = r, d = d, full = full)
}
figB2_terrain <- function() {
  # Purpose: show how steeply the ranges rise behind the leeward coast (Section 1.4) and
  # where the Honolulu sensors sit relative to the valley mouths (Section 5.2).
  I <- dem_layers(read_dem_txt(file.path(EXT, "oahu_dem_500m_q5.txt")), 250)
  H <- dem_layers(read_dem_txt(file.path(EXT, "honolulu_dem_150m_q5.txt")), 75)
  zleg <- scale_colour_gradientn(colours = c("#dfe9cf", "#cfe0b4", "#e5dcaa", "#dcc594", "#c9a77c", "#ad8a67", "#8f7157", "#7a6352"),
                                 values = c(0, 50, 150, 300, 500, 750, 1000, 1250) / 1250, limits = c(0, 1250),
                                 name = "Elevation (m)", breaks = c(0, 250, 500, 750, 1000, 1250))
  ib <- st_bbox(island); bI <- bounds4(ib$xmin - 1500, ib$ymin - 1500, ib$xmax + 1500, ib$ymax + 1500)
  boxes <- bind_rows(lapply(names(EQB), function(r) tibble(region = r, xmin = EQB[[r]]["xmin"], xmax = EQB[[r]]["xmax"],
                                                           ymin = EQB[[r]]["ymin"], ymax = EQB[[r]]["ymax"])))
  lab <- function(tb) bind_cols(tb, as_tibble(t(sapply(seq_len(nrow(tb)), function(i) lonlat_to_xy(tb$lon[i], tb$lat[i]))), .name_repair = ~ c("x", "y")))
  rl <- lab(tibble(name = c("KOʻOLAU  RANGE", "WAIʻANAE  RANGE"), lon = c(-157.905, -158.165), lat = c(21.475, 21.475), rot = c(-47, -70)))
  pk <- lab(tibble(name = c("Kaʻala", "Konahuanui"), lon = c(-158.1453, -157.7939), lat = c(21.5075, 21.3564)))
  pA <- ggplot() + theme_map() + theme(panel.background = element_rect(fill = C_OCEAN, colour = NA), legend.position = "none",
                                        plot.title = element_text(size = 7, face = "bold")) +
    geom_raster(data = I$d, aes(x, y, fill = col)) + scale_fill_identity() +
    geom_contour(data = I$full, aes(x, y, z = z), breaks = seq(200, 1200, 200), colour = "#5b4a3a", linewidth = 0.15, alpha = 0.6) +
    geom_sf(data = st_boundary(st_geometry(island)), colour = C_COAST, linewidth = 0.3) +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = region), fill = NA, linewidth = 0.6, show.legend = FALSE) +
    scale_colour_manual(values = REG_COL, guide = "none")
  pA <- pA + geom_text(data = rl, aes(x, y, label = name, angle = rot), size = 2.4, colour = "#3b2f25", fontface = "bold") +
    geom_point(data = pk, aes(x, y), shape = 17, size = 1.3, colour = "#3b2f25") +
    geom_text(data = pk, aes(x, y, label = name), size = 2.1, colour = "#3b2f25", vjust = -0.8, fontface = "italic") +
    scale_bar(bI, 10000, loc = c(0.05, 0.05)) + north_arrow(bI, loc = c(0.95, 0.80), size = 0.07) +
    labs(title = "(a) Oʻahu (500m grid; contours every 200m)") + map_coord(bI, graticule = FALSE)
  # district panel: sensors coloured by site-mean dT
  bH <- EQB$Honolulu; ext_h <- st_bbox(H$r)
  bH <- bounds4(max(bH["xmin"], ext_h$xmin + 250), max(bH["ymin"], ext_h$ymin + 250), min(bH["xmax"] + 4000, ext_h$xmax - 350),
                min(bH["ymax"] + 5000, ext_h$ymax - 350))     # inside the projected DEM, whose edges are not straight
  s <- sites[sites$region == "Honolulu", ]
  pl <- lab(tibble(name = c("Mānoa\nValley", "Diamond Head", "Waikīkī", "UH Mānoa"), lon = c(-157.806, -157.806, -157.829, -157.8185),
                   lat = c(21.322, 21.2635, 21.2785, 21.2975)))
  cool <- s %>% filter(sensor_id %in% c("HNL01", "HNL04"))
  pB <- ggplot() + theme_map() + theme(panel.background = element_rect(fill = C_OCEAN, colour = NA), legend.position = "none",
                                        plot.title = element_text(size = 7, face = "bold")) +
    geom_raster(data = H$d, aes(x, y, fill = col)) + scale_fill_identity() +
    geom_contour(data = H$full, aes(x, y, z = z), breaks = c(20, 50, 100, 200, 300, 400, 600, 800), colour = "#5b4a3a", linewidth = 0.15, alpha = 0.6) +
    geom_sf(data = st_boundary(st_geometry(island)), colour = C_COAST, linewidth = 0.3) +
    geom_sf(data = s, aes(colour = dT_site), size = 1.9) +
    scale_colour_distiller(palette = "RdBu", limits = c(-1.6, 1.6), oob = scales::squish, name = "Site-mean ΔT (°C)") +
    geom_sf(data = s, shape = 21, size = 1.9, colour = "black", fill = NA, stroke = 0.25) +
    geom_sf_text(data = cool, aes(label = sensor_id), size = 2.0, nudge_x = 420, nudge_y = -120, fontface = "bold") +
    geom_label(data = pl, aes(x, y, label = name), size = 2.1, fontface = "italic", label.size = 0, fill = alpha("white", 0.6),
               label.padding = unit(0.05, "lines"), lineheight = 0.85) +
    scale_bar(bH, 1000, loc = c(0.05, 0.05)) + north_arrow(bH, loc = c(0.94, 0.12), size = 0.07) +
    labs(title = "(b) Honolulu district (150m grid; contours at 20, 50, 100m, then every 100m)") + map_coord(bH, graticule = FALSE)
  # elevation colour bar for both panels (a dummy layer outside the map)
  legz <- cowplot::get_legend(ggplot(tibble(x = 1:2, y = 1:2, z = c(0, 1250)), aes(x, y, colour = z)) + geom_point() + zleg +
                                theme_void() + theme(legend.position = "bottom", legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.2, "cm"),
                                                     legend.title = element_text(size = 6.5), legend.text = element_text(size = 6)) +
                                guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)))
  legt <- cowplot::get_legend(ggplot(tibble(x = 1:2, y = 1:2, z = c(-1.6, 1.6)), aes(x, y, colour = z)) + geom_point() +
                                scale_colour_distiller(palette = "RdBu", limits = c(-1.6, 1.6), name = "Honolulu sites: site-mean ΔT (°C)") +
                                theme_void() + theme(legend.position = "bottom", legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.2, "cm"),
                                                     legend.title = element_text(size = 6.5), legend.text = element_text(size = 6)) +
                                guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)))
  wA <- 3.3; wB <- 3.9
  h <- max(map_aspect(bI) * wA, map_aspect(bH) * wB)
  p <- (pA | pB) / (patchwork::wrap_elements(legz) | patchwork::wrap_elements(legt)) +
    plot_layout(heights = c(h, 0.5)) & theme(plot.margin = margin(2, 3, 2, 3))
  p[[1]] <- p[[1]] + plot_layout(widths = c(wA, wB))
  save_fig(p, "B2_terrain", width = 7.2, height = h + 0.95)
}

# ==============================================================================
# B3  Sky view factor
svf_rays <- function(pt, bldg, z0 = 1.5, R = 200, n_az = 36) {
  # as sky_view_factor() in 01_site_predictors.R, keeping the angle and position of the
  # highest obstruction along each ray
  idx <- st_intersects(st_buffer(pt, R, nQuadSegs = 16), bldg)[[1]]
  cand <- bldg[idx, ]; xy <- st_coordinates(pt)
  bind_rows(lapply(seq_len(n_az) - 1, function(k) {
    th <- 2 * pi * k / n_az; end <- xy + R * c(sin(th), cos(th))
    ray <- st_sfc(st_linestring(rbind(xy, end)), crs = st_crs(pt))
    hit <- cand[st_intersects(ray, cand)[[1]], ]
    beta <- 0; hx <- NA_real_; hy <- NA_real_; hd <- NA_real_; hh <- NA_real_
    if (nrow(hit)) for (j in seq_len(nrow(hit))) {
      seg <- suppressWarnings(st_intersection(ray, st_geometry(hit)[j]))
      d <- max(1, as.numeric(st_distance(pt, seg)))
      b <- atan(max(hit$height_m[j] - z0, 0) / d)
      if (b > beta) { beta <- b; hd <- d; hh <- hit$height_m[j]; p1 <- xy + d * c(sin(th), cos(th)); hx <- p1[1]; hy <- p1[2] }
    }
    tibble(az = 360 * k / n_az, beta = beta, x1 = end[1], y1 = end[2], hx = hx, hy = hy, d = hd, h = hh)
  }))
}
figB3_svf <- function(site_id = "ID26") {
  # Purpose: show how the sky view factor of a site is calculated from the footprints.
  s <- sites[sites$sensor_id == site_id, ]; xy <- st_coordinates(s)
  b <- bounds4(xy[1] - 215, xy[2] - 215, xy[1] + 215, xy[2] + 215)
  bl <- load_buildings(b) %>% mutate(height_m = pmin(pmax(coalesce(height_m, 3), 0.5), 60))
  rays <- svf_rays(st_geometry(s), bl)
  svf <- mean(cos(rays$beta)^2)
  stopifnot(abs(svf - s$svf_point) < 1e-6)                    # same value as the descriptor used in the models
  # (a) the angle beta in a vertical section
  H <- 14; D <- 18; z0 <- 1.5
  pa <- ggplot() + theme_void() + theme(plot.title = element_text(size = 7.5, face = "bold")) +
    annotate("segment", x = -4, xend = 26, y = 0, yend = 0, linewidth = 0.5, colour = "#555555") +
    annotate("rect", xmin = D, xmax = D + 6, ymin = 0, ymax = H, fill = "#9e9e9e", colour = "#555555", linewidth = 0.3) +
    annotate("segment", x = 0, xend = 0, y = 0, yend = z0 + 0.4, linewidth = 0.8, colour = "#333333") +
    annotate("point", x = 0, y = z0, size = 2.2, colour = "#c0392b") +
    annotate("segment", x = 0, xend = D, y = z0, yend = H, linetype = "dashed", linewidth = 0.4) +
    annotate("segment", x = 0, xend = D, y = z0, yend = z0, linewidth = 0.3, colour = "#555555") +
    annotate("curve", x = 5.5, xend = 5.2, y = z0, yend = z0 + 5.2 * (H - z0) / D, curvature = 0.3, linewidth = 0.3) +
    annotate("text", x = 6.6, y = z0 + 1.4, label = "β", size = 3.4, fontface = "italic") +
    annotate("text", x = 0, y = -1.1, label = "sensor\n(1.5m)", size = 2.2, vjust = 1, lineheight = 0.9) +
    annotate("text", x = D / 2, y = z0 - 0.9, label = "d", size = 2.6, fontface = "italic") +
    annotate("text", x = D + 3, y = H + 1.2, label = "building, height H", size = 2.2) +
    annotate("text", x = 11, y = 17.5, label = "tan β = (H − 1.5m) / d", size = 2.5) +
    coord_equal(xlim = c(-5, 27), ylim = c(-4.5, 19.5), expand = FALSE) +
    labs(title = "(a) The obstruction angle β")
  # (b) plan view of the 36 rays
  bb <- st_as_sfc(st_bbox(b, crs = st_crs(CRS_M)))
  circ <- st_buffer(st_geometry(s), 200, nQuadSegs = 32)
  blc <- suppressWarnings(st_intersection(bl, bb))
  pb <- ggplot() + theme_map() + theme(plot.title = element_text(size = 7.5, face = "bold"), legend.position = "bottom",
                                        legend.key.width = unit(0.55, "cm"), legend.key.height = unit(0.15, "cm"),
                                        legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
                                        plot.margin = margin(2, 8, 2, 2)) +
    geom_sf(data = bb, fill = "white", colour = NA) +
    geom_sf(data = blc, aes(fill = pmin(height_m, 40)), colour = "#666666", linewidth = 0.1) +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 40), name = "Building height (m)", oob = scales::squish) +
    geom_sf(data = circ, fill = NA, colour = "#1f3b57", linewidth = 0.3, linetype = "dashed") +
    geom_segment(data = rays, aes(x = xy[1], y = xy[2], xend = x1, yend = y1), colour = "#1f3b57", linewidth = 0.15, alpha = 0.7) +
    geom_point(data = rays %>% filter(!is.na(hx)), aes(hx, hy), colour = "#c0392b", size = 0.7) +
    geom_point(aes(x = xy[1], y = xy[2]), colour = "black", size = 1.4) +
    scale_bar(b, 100, loc = c(0.05, 0.05), label = "100m") +
    labs(title = sprintf("(b) Site %s: 36 directions to 200m", site_id)) + map_coord(b, graticule = FALSE)
  # (c) beta by direction
  pc <- ggplot(rays, aes(az, beta * 180 / pi)) + theme_thesis() +
    geom_col(fill = "#c0392b", width = 7) +
    scale_x_continuous(breaks = seq(0, 330, 90), labels = c("N", "E", "S", "W"), limits = c(-5, 355)) +
    labs(x = "Direction", y = "Obstruction angle β (°)", title = "(c) β in each direction",
         subtitle = sprintf("SVF = mean of cos²β over the 36 directions = %.2f", svf)) +
    theme(plot.subtitle = element_text(size = 6.5))
  p <- (pa | pb | pc) + plot_layout(widths = c(1, 1.05, 0.95))
  save_fig(p, "B3_svf", width = 7.2, height = 3.0)
  invisible(svf)
}

main <- function(which = c("B1", "B2", "B3")) {
  if ("B1" %in% which) figB1_landcover()
  if ("B2" %in% which) figB2_terrain()
  if ("B3" %in% which) figB3_svf()
}
if (sys.nframe() == 0) { a <- commandArgs(trailingOnly = TRUE); main(if (length(a)) a else c("B1", "B2", "B3")) }
