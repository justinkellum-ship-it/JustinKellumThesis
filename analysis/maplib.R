# maplib.R --------------------------------------------------------------------
# Shared cartographic helpers for the thesis maps (ggplot2 + sf + terra).
# All map layers are drawn in EPSG:6634 (NAD83(PA11) / UTM zone 4N, metres).
#
#   load_island()        Oʻahu land polygon (GSHHG full resolution)
#   load_sites()         one sf point per logger site with the 100 m predictors
#   domain_bounds()      bounding box of a district's network plus a margin
#   tint_raster()        C-CAP mask cropped to a box and aggregated to `block` m
#                        (cached as GeoTIFF under results/cache/)
#   basemap()            list of ggplot2 layers: ocean, land, paved / canopy /
#                        water tints, building footprints, coastline
#   scale_bar(), north_arrow(), place_labels(), panel_tag()
#
# Packages: sf (vector data), terra (rasters: vrt, crop, aggregate, subst),
#           ggplot2 (drawing), patchwork (multi-panel layout).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(sf); library(terra); library(ggplot2); library(patchwork) })
Sys.setenv(PROJ_NETWORK = "OFF"); sf_proj_network(FALSE); sf_use_s2(FALSE)

CRS_M <- 6634
C_OCEAN <- "#dbe9f4"; C_LAND <- "#f7f5f0"; C_PAVED <- "#c9c9c9"; C_BLDG <- "#6e6e6e"
C_TREE <- "#8fbf7f"; C_WATER <- "#a9c8e6"; C_COAST <- "#4f6d8a"
REG_COL <- c(Honolulu = "#c0392b", Ewa = "#2471a3")
REG_LABEL <- c(Honolulu = "Honolulu", Ewa = "ʻEwa")
CACHE <- file.path(OUT, "cache"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

theme_thesis <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
          legend.title = element_text(size = base_size - 1), legend.text = element_text(size = base_size - 1.5),
          legend.key.size = unit(0.35, "cm"), legend.background = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = base_size, face = "bold", hjust = 0))
}
theme_map <- function(base_size = 8) {
  theme_void(base_size = base_size) +
    theme(plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
          legend.title = element_text(size = base_size - 1), legend.text = element_text(size = base_size - 1.5),
          panel.border = element_rect(fill = NA, colour = "#444444", linewidth = 0.4))
}

load_island <- function() st_transform(st_read(file.path(EXT, "oahu_gshhs_f.geojson"), quiet = TRUE), CRS_M)

load_sites <- function(radius = 100) {
  sn <- read_csv(file.path(PROC, "sensor_nights.csv"), col_types = cols(night_date = col_character(), .default = col_guess()))
  pr <- read_csv(file.path(DER, "site_predictors_multiscale.csv"), show_col_types = FALSE) %>%
    filter(radius_m == radius) %>% select(-radius_m, -region)
  sm <- read_csv(file.path(TAB, "site_summary.csv"), show_col_types = FALSE) %>% select(sensor_id, dT_site, dT_sd, n_nights)
  s <- sn %>% distinct(sensor_id, .keep_all = TRUE) %>% select(sensor_id, region, latitude, longitude) %>%
    inner_join(pr, by = "sensor_id") %>% left_join(sm, by = "sensor_id")
  st_transform(st_as_sf(s, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_M)
}

bounds4 <- function(xmin, ymin, xmax, ymax) c(xmin = unname(xmin), ymin = unname(ymin), xmax = unname(xmax), ymax = unname(ymax))
domain_bounds <- function(sites, region, pad = 600) {
  b <- st_bbox(sites[sites$region == region, ])
  bounds4(b$xmin - pad, b$ymin - pad, b$xmax + pad, b$ymax + pad)
}
equal_bounds <- function(dom) {
  # both district boxes expanded to the same width and height (same map scale)
  w <- max(sapply(dom, function(b) b["xmax"] - b["xmin"])); h <- max(sapply(dom, function(b) b["ymax"] - b["ymin"]))
  lapply(dom, function(b) { cx <- (b["xmin"] + b["xmax"]) / 2; cy <- (b["ymin"] + b["ymax"]) / 2
    c(xmin = unname(cx - w / 2), ymin = unname(cy - h / 2), xmax = unname(cx + w / 2), ymax = unname(cy + h / 2)) })
}

ccap_vrt <- function(kind) {
  files <- list.files(file.path(EXT, "ccap"), pattern = paste0(kind, ".*\\.tif$"), full.names = TRUE)
  vrt(files, filename = file.path(tempdir(), paste0("ccap_", kind, ".vrt")), overwrite = TRUE)
}

tint_raster <- function(kind, value, bounds, block, key) {
  # fraction of 1 m pixels == value inside `bounds`, on a `block` m grid; cached
  f <- file.path(CACHE, sprintf("%s_%s_%dm.tif", key, kind, block))
  if (file.exists(f)) return(rast(f))
  r <- crop(ccap_vrt(kind), ext(bounds["xmin"], bounds["xmax"], bounds["ymin"], bounds["ymax"]))
  r <- subst(r, NA, 0)                       # the files' no-data value 0 = "not this class"
  r <- aggregate(r == value, fact = block, fun = "mean")
  writeRaster(r, f, overwrite = TRUE, datatype = "FLT4S")
  rast(f)
}
raster_df <- function(r, name = "v") { d <- as.data.frame(r, xy = TRUE); names(d)[3] <- name; d$x <- round(d$x, 3); d$y <- round(d$y, 3); d }

load_buildings <- function(bounds) {
  b <- st_read(file.path(EXT, "buildings_oahu.gpkg"), quiet = TRUE,
               wkt_filter = st_as_text(st_transform(st_as_sfc(st_bbox(c(bounds), crs = st_crs(CRS_M))), 32604)))
  st_transform(b, CRS_M)
}

basemap <- function(bounds, island, buildings = NULL, key = "dom", block = 4,
                    canopy = TRUE, paved = TRUE, water = TRUE, alpha_tree = 0.75) {
  # returns a list of ggplot2 layers; add coord_sf() afterwards
  land <- st_intersection(st_geometry(island), st_as_sfc(st_bbox(bounds, crs = st_crs(CRS_M))))
  layers <- list(
    theme(panel.background = element_rect(fill = C_OCEAN, colour = NA)),
    geom_sf(data = land, fill = C_LAND, colour = NA))
  if (paved) {
    d <- raster_df(tint_raster("impervious", 1, bounds, block, key))
    layers <- c(layers, list(geom_raster(data = d, aes(x, y, alpha = v), fill = C_PAVED, interpolate = TRUE),
                             scale_alpha_identity()))
  }
  if (canopy) {
    d <- raster_df(tint_raster("canopy", 1, bounds, block, key)); d$v <- d$v * alpha_tree
    layers <- c(layers, list(geom_raster(data = d, aes(x, y, alpha = v), fill = C_TREE, interpolate = TRUE)))
  }
  if (water) {
    d <- raster_df(tint_raster("water", 1, bounds, block, key)) %>% filter(v > 0.5)
    layers <- c(layers, list(geom_raster(data = d, aes(x, y), fill = C_WATER)))
  }
  if (!is.null(buildings) && nrow(buildings))
    layers <- c(layers, list(geom_sf(data = st_geometry(buildings), fill = C_BLDG, colour = C_BLDG, linewidth = 0.05)))
  c(layers, list(geom_sf(data = st_boundary(st_geometry(island)), colour = C_COAST, linewidth = 0.3)))
}

map_coord <- function(bounds, graticule = TRUE) {
  if (graticule)
    coord_sf(xlim = c(bounds["xmin"], bounds["xmax"]), ylim = c(bounds["ymin"], bounds["ymax"]), expand = FALSE,
             crs = st_crs(CRS_M), datum = st_crs(4326), label_axes = list(bottom = "E", left = "N"))
  else
    coord_sf(xlim = c(bounds["xmin"], bounds["xmax"]), ylim = c(bounds["ymin"], bounds["ymax"]), expand = FALSE,
             crs = st_crs(CRS_M), datum = NA)
}
graticule_breaks <- function(step = 0.02) list(scale_x_continuous(breaks = seq(-158.30, -157.60, step)), scale_y_continuous(breaks = seq(21.20, 21.50, step)))
map_aspect <- function(bounds) unname((bounds["ymax"] - bounds["ymin"]) / (bounds["xmax"] - bounds["xmin"]))

scale_bar <- function(bounds, length_m = 1000, loc = c(0.05, 0.05), label = NULL) {
  x0 <- bounds["xmin"] + loc[1] * (bounds["xmax"] - bounds["xmin"]); y0 <- bounds["ymin"] + loc[2] * (bounds["ymax"] - bounds["ymin"])
  h <- 0.008 * (bounds["ymax"] - bounds["ymin"])
  list(annotate("segment", x = x0, xend = x0 + length_m, y = y0, yend = y0, linewidth = 0.6),
       annotate("segment", x = c(x0, x0 + length_m), xend = c(x0, x0 + length_m), y = y0, yend = y0 + h, linewidth = 0.6),
       annotate("label", x = x0 + length_m / 2, y = y0 + 1.6 * h, label = label %||% sprintf("%gkm", length_m / 1000),
                size = 2.2, label.size = 0, fill = alpha("white", 0.7), label.padding = unit(0.08, "lines"), vjust = 0))
}
north_arrow <- function(bounds, loc = c(0.94, 0.86), size = 0.07) {
  x <- bounds["xmin"] + loc[1] * (bounds["xmax"] - bounds["xmin"]); y <- bounds["ymin"] + loc[2] * (bounds["ymax"] - bounds["ymin"])
  L <- size * (bounds["ymax"] - bounds["ymin"])
  list(annotate("segment", x = x, xend = x, y = y, yend = y + L, linewidth = 0.5, arrow = arrow(length = unit(0.12, "cm"), type = "closed")),
       annotate("label", x = x, y = y + 1.12 * L, label = "N", size = 2.6, fontface = "bold", label.size = 0,
                fill = alpha("white", 0.7), label.padding = unit(0.05, "lines"), vjust = 0))
}
lonlat_to_xy <- function(lon, lat) st_coordinates(st_transform(st_sfc(st_point(c(lon, lat)), crs = 4326), CRS_M))

PLACES <- list(
  Honolulu = tribble(~name, ~lon, ~lat,
    "Waikīkī", -157.8267, 21.2793, "Ala Moana", -157.8437, 21.2911, "Kakaʻako", -157.8560, 21.2985,
    "McCully–Mōʻiliʻili", -157.8290, 21.2925, "Kapiʻolani\nPark", -157.8225, 21.2690, "UH Mānoa", -157.8170, 21.2975,
    "Diamond Head", -157.8060, 21.2620, "Māmala Bay", -157.8450, 21.2760),
  Ewa = tribble(~name, ~lon, ~lat,
    "Kapolei", -158.0580, 21.3355, "ʻEwa Beach", -158.0080, 21.3155, "Kalaeloa", -158.0700, 21.3120,
    "Honouliuli", -158.0365, 21.3515, "ʻEwa Villages", -158.0330, 21.3350, "Makakilo", -158.0850, 21.3600,
    "West Loch", -158.0175, 21.3585, "Barbers Point\nHarbor", -158.1130, 21.3220, "Campbell\nIndustrial Park", -158.0870, 21.3100))

place_labels <- function(region, bounds, size = 2.1) {
  p <- PLACES[[region]]; xy <- t(sapply(seq_len(nrow(p)), function(i) lonlat_to_xy(p$lon[i], p$lat[i])))
  p$x <- xy[, 1]; p$y <- xy[, 2]
  p <- p %>% filter(x > bounds["xmin"], x < bounds["xmax"], y > bounds["ymin"], y < bounds["ymax"])
  geom_label(data = p, aes(x, y, label = name), size = size, fontface = "italic", colour = "#222222", label.size = 0,
             fill = alpha("white", 0.65), label.padding = unit(0.06, "lines"), lineheight = 0.85)
}
panel_tag <- function(bounds, text, loc = c(0.02, 0.975), size = 3.2) {
  annotate("label", x = bounds["xmin"] + loc[1] * (bounds["xmax"] - bounds["xmin"]), y = bounds["ymin"] + loc[2] * (bounds["ymax"] - bounds["ymin"]),
           label = text, hjust = 0, vjust = 1, size = size, fontface = "bold", label.size = 0, fill = alpha("white", 0.75))
}

save_fig <- function(p, name, width = 7.2, height = 4, dpi = 300) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = width, height = height, dpi = dpi, bg = "white", device = grDevices::png, type = "cairo")
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf, bg = "white")
  message("saved ", name)
}
