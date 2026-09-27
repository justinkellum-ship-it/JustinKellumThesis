#!/usr/bin/env python3
"""
fetch_weather_series.py
-----------------------
Retrieves data/external/weather_hourly_series.csv: hourly 10 m wind speed (km/h) and
total cloud cover (%) for the whole deployment period at the centroid of each district,
from the Open-Meteo historical weather API (https://open-meteo.com/en/docs/historical-weather-api;
Zippenfenig 2023).  The API's default historical dataset combines the ECMWF IFS
high-resolution analysis (9 km), ERA5 and ERA5-Land; for the study period and place the
values it returns are those of the IFS analysis, which was verified hour by hour against each
dataset the API offers (the ERA5 dataset alone differs by up to 21 km/h).

The Open-Meteo Python client (openmeteo-requests) is used because it returns the values at
full (float32) precision; the JSON interface rounds them to one decimal.  Rerunning this
script reproduces the committed file byte for byte (checked 27 September 2026).

The file holds, for every distinct time stamp of a district's loggers (data/raw/logger_readings.csv),
the hourly value of the nearest hour; 00_prepare_inputs.R also accepts a plain hourly series
(one row per hour and district).  Time stamps are local standard time (Pacific/Honolulu has no
daylight saving).

Usage:  pip install openmeteo-requests pandas
        python analysis/inputs/fetch_weather_series.py
"""
import os

import openmeteo_requests
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "data", "external", "weather_hourly_series.csv")
TIMEZONE = "Pacific/Honolulu"

sites = pd.read_csv(os.path.join(RAW, "sites.csv"))
readings = pd.read_csv(os.path.join(RAW, "logger_readings.csv"), usecols=["region", "datetime"])

client = openmeteo_requests.Client()
parts = []
for region, g in sites.groupby("region", sort=False):
    stamps = pd.to_datetime(readings.datetime[readings.region == region].unique())   # every time stamp of the district's loggers
    start, end = stamps.min().strftime("%Y-%m-%d"), stamps.max().strftime("%Y-%m-%d")
    lat, lon = g.latitude.mean(), g.longitude.mean()                    # the district centroid
    r = client.weather_api("https://archive-api.open-meteo.com/v1/archive", params={
        "latitude": lat, "longitude": lon, "start_date": start, "end_date": end,
        "hourly": ["wind_speed_10m", "cloud_cover"], "wind_speed_unit": "kmh", "timezone": TIMEZONE})[0]
    h = r.Hourly()
    t = pd.date_range(start=pd.to_datetime(h.Time(), unit="s", utc=True), end=pd.to_datetime(h.TimeEnd(), unit="s", utc=True),
                      freq=pd.Timedelta(seconds=h.Interval()), inclusive="left").tz_convert(TIMEZONE).tz_localize(None)
    hourly = pd.DataFrame({"hour": t, "wind": h.Variables(0).ValuesAsNumpy(), "cloud": h.Variables(1).ValuesAsNumpy()})
    # the series used by the analysis: the hourly value at the nearest hour, attached to every logger
    # time stamp of the district (pandas rounds a half hour to the even hour, as the original retrieval did)
    at = pd.DataFrame({"datetime": stamps.sort_values(), "hour": stamps.sort_values().round("h")}).merge(hourly, on="hour", how="left")
    parts.append(pd.DataFrame({"datetime": at.datetime.dt.strftime("%Y-%m-%d %H:%M:%S"), "wind": at.wind, "cloud": at.cloud, "region": region}))
    print(f"{region}: centroid {lat:.4f}, {lon:.4f} -> grid point {r.Latitude():.4f}, {r.Longitude():.4f}; "
          f"{len(hourly)} hours {start} to {end}, attached to {len(at)} logger time stamps")
out = pd.concat(parts, ignore_index=True)
assert not out[["wind", "cloud"]].isna().any().any()
out.to_csv(OUT, index=False)                                            # float32 values, shortest exact form
print(f"wrote {OUT}: {len(out)} rows")
