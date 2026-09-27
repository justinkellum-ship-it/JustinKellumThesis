#!/usr/bin/env python3
"""
site_properties.py
------------------
Function library for 01_site_predictors.py: the surface-property descriptors
of a source area (Stewart & Oke 2012) from the same public data as the R
pipeline, so that the two implementations can be compared:

  * NOAA C-CAP (2021) 1 m land-cover masks for Hawaiʻi
      - impervious (1 = impervious)
      - canopy     (1 = tree canopy, 2 = shrub/scrub)
      - water      (1 = open water)
  * FEMA / ORNL USA Structures building footprints (Oʻahu subset, with height)

Not run on its own.
"""
import glob
import json
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from shapely.geometry import Point, Polygon, mapping, box
from shapely.ops import unary_union

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
EXT = os.path.join(DATA, "external")
PDATA = os.path.join(ROOT, "replication", "results", "data")

CCAP_DIR = os.path.join(EXT, "ccap")
BUILDINGS = os.path.join(EXT, "buildings_oahu.gpkg")
SENSOR_NIGHTS = os.path.join(PDATA, "sensor_nights.csv")

RASTER_CRS = "EPSG:6634"     # NAD83(PA11) / UTM zone 4N (C-CAP tiles)
RADIUS_M = 200.0             # source-area radius (Stewart & Oke 2012)

TILES = {
    "impervious": sorted(glob.glob(os.path.join(CCAP_DIR, "*impervious*.tif"))),
    "canopy": sorted(glob.glob(os.path.join(CCAP_DIR, "*canopy*.tif"))),
    "water": sorted(glob.glob(os.path.join(CCAP_DIR, "*water*.tif"))),
}


# ----------------------------------------------------------------------
def raster_fractions(geom, tile_paths, values):
    """Fraction (%) of pixels inside `geom` equal to each value in `values`."""
    counts = {v: 0 for v in values}
    total = 0
    for tp in tile_paths:
        with rasterio.open(tp) as src:
            b = src.bounds
            gx0, gy0, gx1, gy1 = geom.bounds
            if gx1 < b.left or gx0 > b.right or gy1 < b.bottom or gy0 > b.top:
                continue
            win = src.window(*geom.bounds).round_offsets().round_lengths()
            # clip window to raster extent
            win = win.intersection(rasterio.windows.Window(0, 0, src.width, src.height))
            if win.width <= 0 or win.height <= 0:
                continue
            data = src.read(1, window=win)
            tr = src.window_transform(win)
            m = geometry_mask([mapping(geom)], out_shape=data.shape, transform=tr, invert=True)
            vals = data[m]
            total += vals.size
            for v in values:
                counts[v] += int((vals == v).sum())
    if total == 0:
        return {v: np.nan for v in values}, 0
    return {v: 100.0 * counts[v] / total for v in values}, total


def building_metrics(geom, bldg, sindex):
    """Building descriptors for the source area `geom` (projected metres)."""
    area = geom.area
    idx = sindex.query(geom, predicate="intersects")
    if len(idx) == 0:
        return dict(bldg_frac=0.0, zH=0.0, zH_arith=0.0, zH_max=0.0,
                    n_bldg=0, mean_footprint=0.0, large_bldg_share=0.0,
                    hw_ratio=0.0)
    sub = bldg.iloc[idx]
    clipped = sub.geometry.intersection(geom)
    a = clipped.area.values
    h = sub["height_m"].values
    full = sub["area_m2"].values
    keep = a > 0
    a, h, full = a[keep], h[keep], full[keep]
    bldg_area = a.sum()
    bldg_frac = 100.0 * bldg_area / area
    zH = float((a * h).sum() / a.sum()) if a.sum() > 0 else 0.0
    zH_arith = float(h.mean()) if len(h) else 0.0
    zH_max = float(h.max()) if len(h) else 0.0
    n = int(keep.sum())
    mean_fp = float(full.mean()) if n else 0.0
    large_share = float(a[full >= 1000].sum() / a.sum()) if a.sum() > 0 else 0.0
    # canyon aspect ratio H/W from mean building spacing (Stewart & Oke 2012
    # style approximation)
    avg_fp_in = bldg_area / n
    spacing = np.sqrt(area / n) - np.sqrt(avg_fp_in)
    W = max(6.0, spacing)
    hw = min(5.0, zH / W)
    return dict(bldg_frac=bldg_frac, zH=zH, zH_arith=zH_arith, zH_max=zH_max,
                n_bldg=n, mean_footprint=mean_fp, large_bldg_share=large_share,
                hw_ratio=hw)


def describe(geom, bldg, sindex):
    fr, npx = raster_fractions(geom, TILES["impervious"], [1])
    imp = fr[1]
    fr, _ = raster_fractions(geom, TILES["canopy"], [1, 2])
    tree, shrub = fr[1], fr[2]
    fr, _ = raster_fractions(geom, TILES["water"], [1])
    water = fr[1]
    bm = building_metrics(geom, bldg, sindex)
    perv = max(0.0, 100.0 - imp - water)
    # approximate sky-view factor at canyon floor for a symmetric canyon
    svf = float(np.cos(np.arctan(2.0 * bm["hw_ratio"]))) if bm["hw_ratio"] > 0 else 1.0
    return dict(imperv_frac=imp, tree_frac=tree, shrub_frac=shrub, water_frac=water,
                perv_frac=perv, svf_est=svf, n_pixels=npx, **bm)


# ----------------------------------------------------------------------
