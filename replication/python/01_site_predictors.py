#!/usr/bin/env python3
"""
01_site_predictors.py
---------------------
Two- and three-dimensional urban-form descriptors for every sensor site at
three source-area radii (50, 100, 200 m), from the same public inputs used
in site_properties.py (the Python counterpart of analysis/01_site_predictors.R):

  * NOAA C-CAP 2021 1 m masks  -> impervious, tree-canopy and water fractions
  * FEMA USA Structures         -> building footprint fraction, area-weighted
                                   mean building height, canyon aspect ratio
  * 3-D ray casting on the building footprints -> sky view factor at the
    sensor (1.5 m above ground), 36 azimuths, 200 m horizon search
  * GSHHG shoreline             -> distance to the ocean

Output: replication/results/data/site_predictors_multiscale.csv (one row per site x radius)
"""
import importlib.util
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point, LineString
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(ROOT, "data")
EXT = os.path.join(DATA, "external")
PDATA = os.path.join(ROOT, "replication", "results", "data")

spec = importlib.util.spec_from_file_location("sp", os.path.join(HERE, "site_properties.py"))
sp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp)

RADII = [50, 100, 200]
SENSOR_Z = 1.5          # sensor height above ground (m)
N_AZ = 36               # azimuths for the sky-view factor
HORIZON_R = 200.0       # horizon search radius (m)


def sky_view_factor(pt, bldg, sindex, z0=SENSOR_Z, R=HORIZON_R, n_az=N_AZ):
    """Sky view factor of a horizontal surface at height z0 above ground,
    SVF = (1/N) * sum cos^2(beta_i) over N azimuths, where beta_i is the
    elevation angle of the highest building obstruction along azimuth i
    (the cos² β relation of Johnson & Watson 1984 for a horizontal surface)."""
    idx = sindex.query(pt.buffer(R), predicate="intersects")
    if len(idx) == 0:
        return 1.0
    cand = bldg.iloc[idx]
    x, y = pt.x, pt.y
    total = 0.0
    for k in range(n_az):
        th = 2 * np.pi * k / n_az
        ray = LineString([(x, y), (x + R * np.sin(th), y + R * np.cos(th))])
        beta = 0.0
        hit = cand[cand.intersects(ray)]
        for geom, h in zip(hit.geometry, hit.height_m):
            d = max(1.0, pt.distance(ray.intersection(geom)))
            hh = max(h - z0, 0.0)
            beta = max(beta, np.arctan(hh / d))
        total += np.cos(beta) ** 2
    return total / n_az


def main():
    sn = pd.read_csv(os.path.join(PDATA, "sensor_nights.csv"))
    sites = (sn.drop_duplicates("sensor_id")[["sensor_id", "region", "latitude", "longitude"]]
             .sort_values(["region", "sensor_id"]).reset_index(drop=True))
    g = gpd.GeoDataFrame(sites, geometry=[Point(x, y) for x, y in zip(sites.longitude, sites.latitude)],
                         crs="EPSG:4326").to_crs(sp.RASTER_CRS)
    print("Loading buildings ...")
    bldg = gpd.read_file(sp.BUILDINGS).to_crs(sp.RASTER_CRS)
    bldg["height_m"] = bldg["height_m"].fillna(3.0).clip(lower=0.5, upper=60.0)
    sindex = bldg.sindex
    coast = gpd.read_file(os.path.join(EXT, "oahu_gshhs_f.geojson")).to_crs(sp.RASTER_CRS)
    coast_line = unary_union(coast.geometry).boundary

    rows = []
    for _, r in g.iterrows():
        svf = sky_view_factor(r.geometry, bldg, sindex)
        coast_d = float(r.geometry.distance(coast_line))
        for R in RADII:
            d = sp.describe(r.geometry.buffer(R), bldg, sindex)
            rows.append(dict(sensor_id=r.sensor_id, region=r.region, radius_m=R,
                             imperv=d["imperv_frac"], tree=d["tree_frac"], water=d["water_frac"],
                             bldg=d["bldg_frac"], height=d["zH"], height_max=d["zH_max"],
                             n_bldg=d["n_bldg"], aspect=d["hw_ratio"], svf_point=svf, coast_km=coast_d / 1000.0))
        print(f"  {r.sensor_id:8s} svf={svf:.2f} coast={coast_d:6.0f} m")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(PDATA, "site_predictors_multiscale.csv"), index=False)
    print("wrote", len(out), "rows")
    print(out.groupby(["region", "radius_m"])[["imperv", "tree", "bldg", "height", "svf_point", "aspect"]].mean().round(2))


if __name__ == "__main__":
    main()
