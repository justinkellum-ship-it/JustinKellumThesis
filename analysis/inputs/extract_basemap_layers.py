#!/usr/bin/env python3
"""
extract_basemap_layers.py
-------------------------
Makes the two Natural Earth base-map layers in data/external from the full 10 m data sets
(https://www.naturalearthdata.com, public domain), which are far too large to commit:

  ne_10m_land_hawaii.geojson   the features of "Land" (10 m physical vectors, ne_10m_land)
                               that intersect a window around the Hawaiian Islands: the
                               land polygons drawn in the study-area map (Figure 1)
  ne_10m_roads_oahu.geojson    the features of "Roads" (10 m cultural vectors, ne_10m_roads)
                               that intersect a window around Oʻahu: the highways drawn in
                               the study-area map

The full layers are read as GeoJSON (convert the downloaded shapefiles once, e.g.
`ogr2ogr -f GeoJSON ne_10m_land.geojson ne_10m_land.shp`); the selected features are written
unchanged, with their attributes.  Nothing is clipped: whole features are kept, so the land
file also carries the other coasts of the same multipolygons.

Usage:  python analysis/inputs/extract_basemap_layers.py <ne_10m_land.geojson> <ne_10m_roads.geojson>
The outputs are committed; run this only to rebuild them.  Requires shapely.
"""
import json
import os
import sys

from shapely.geometry import box, shape

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXT = os.path.join(ROOT, "data", "external")
WINDOWS = {   # (input, output, selection window as lon/lat bounds)
    "land": ("ne_10m_land_hawaii.geojson", box(-160.5, 18.8, -154.6, 22.3)),    # the main Hawaiian Islands
    "roads": ("ne_10m_roads_oahu.geojson", box(-158.35, 21.2, -157.6, 21.75)),  # Oʻahu
}

def extract(src, name, window):
    full = json.load(open(src, encoding="utf-8"))
    keep = [f for f in full["features"] if f["geometry"] and shape(f["geometry"]).intersects(window)]
    out = {"type": "FeatureCollection", "name": name.replace(".geojson", ""), "features": keep}
    if "crs" in full:
        out["crs"] = full["crs"]
    with open(os.path.join(EXT, name), "w", encoding="utf-8") as f:
        json.dump(out, f)
    print(f"{name}: {len(keep)} of {len(full['features'])} features kept")

if len(sys.argv) != 3:
    sys.exit(__doc__)
extract(sys.argv[1], *WINDOWS["land"])
extract(sys.argv[2], *WINDOWS["roads"])
