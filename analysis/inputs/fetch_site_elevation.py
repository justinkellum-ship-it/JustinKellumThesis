#!/usr/bin/env python3
"""
fetch_site_elevation.py
-----------------------
Retrieves data/descriptors/site_elevation.csv: the bare-earth elevation of every logger
position from the USGS 3DEP 1 m lidar DEM, read through the USGS Elevation Point Query
Service (EPQS, https://epqs.nationalmap.gov, API v1), one query per site in
data/raw/sites.csv.  06_elevation_check.R uses the elevations for the terrain check
(Section 4.8, Table 16 and Figure 19 of the thesis).

Each answer carries the elevation in metres (rounded here to 0.01 m), the id of the DEM
raster it was read from and the raster's resolution in metres; the file keeps all three
and the date of retrieval.  When the file already exists, the new elevations are compared
with the stored ones before it is rewritten (the DEM is static, so a rerun reproduces the
values; only the retrieval date in the `source` column changes).

Usage:  python analysis/inputs/fetch_site_elevation.py      (standard library only)
"""
import csv
import datetime as dt
import json
import os
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SITES = os.path.join(ROOT, "data", "raw", "sites.csv")
OUT = os.path.join(ROOT, "data", "descriptors", "site_elevation.csv")
API = "https://epqs.nationalmap.gov/v1/json"

def query(lat, lon):
    q = urllib.parse.urlencode({"x": lon, "y": lat, "wkid": 4326, "units": "Meters", "includeDate": "false"})
    with urllib.request.urlopen(f"{API}?{q}", timeout=60) as r:
        a = json.load(r)
    return float(a["value"]), int(a["rasterId"]), a["resolution"]

with open(SITES, newline="", encoding="utf-8") as f:
    sites = list(csv.DictReader(f))
today = dt.date.today().isoformat()
rows = []
for s in sites:
    lat, lon = float(s["latitude"]), float(s["longitude"])
    elev, raster_id, res = query(lat, lon)
    rows.append({"sensor_id": s["sensor_id"], "region": s["region"], "latitude": s["latitude"], "longitude": s["longitude"],
                 "elev_m": repr(round(elev, 2)), "raster_id": raster_id,
                 "source": f"USGS 3DEP 1 m DEM via EPQS v1 (epqs.nationalmap.gov), resolution={res}, retrieved {today}"})
    print(f"{s['sensor_id']:6s} {s['region']:9s} {elev:8.2f} m  raster {raster_id}  resolution {res} m")

rows.sort(key=lambda r: (r["region"], r["sensor_id"]))     # the order of the committed file

if os.path.exists(OUT):                                   # a rerun must give the stored elevations
    with open(OUT, newline="", encoding="utf-8") as f:
        old = {r["sensor_id"]: r for r in csv.DictReader(f)}
    assert set(old) == {r["sensor_id"] for r in rows}, "the stored file covers different sites"
    diff = max(abs(float(old[r["sensor_id"]]["elev_m"]) - float(r["elev_m"])) for r in rows)
    same_raster = all(int(old[r["sensor_id"]]["raster_id"]) == r["raster_id"] for r in rows)
    print(f"compared with the stored file: largest elevation difference {diff:.2f} m; same rasters: {same_raster}")

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["sensor_id", "region", "latitude", "longitude", "elev_m", "raster_id", "source"], lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT}: {len(rows)} sites")
