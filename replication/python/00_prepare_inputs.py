#!/usr/bin/env python3
"""
00_prepare_inputs.py
--------------------
Assemble the analysis-ready temperature tables from the complete half-hourly
iButton record (data/raw/logger_readings.csv, every sensor and every
night of both deployments), the hourly wind / cloud series of the two
districts (data/external/weather_hourly_series.csv) and the logger positions
(data/raw/sites.csv).  This is the Python counterpart of analysis/00_prepare_inputs.R.

Outputs (replication/results/data/):
  night_inventory.csv        one row per region x diurnal date (06:00-06:00
                             HST): sensors reporting, night-mean wind and
                             cloud, share of calm & clear night steps, and the
                             inclusion decision with its reason
  halfhourly_calm_clear.csv  one row per sensor x 30-min step for the whole
                             diurnal cycles that contain a qualifying night
  sensor_nights.csv          one row per sensor x qualifying night
                             (18:00-06:00 HST means)

Processing rules (all made explicit here so that they can be varied):
  1. Time binning.  Loggers were started by hand, so their clocks are offset
     from each other by up to +-2 min (readings at :00/:01/:02 and
     :30/:31/:32).  Every reading is assigned to the nearest 30-min bin so
     that all sensors of a network are compared at the same time step (a
     median taken by exact minute stamp would compare each sensor only with
     the sensors that share its cadence).
  2. Reference.  The network median air temperature of all sensors of the
     region reporting in the bin; dT = T_sensor - T_median.
  3. Night.  18:00-06:00 HST; a diurnal date runs 06:00-06:00 and is labelled
     by the date of its first 06:00.
  4. Calm & clear night (adopted from the thesis, after Stewart & Oke 2012 /
     Oke 1982 for the reasoning): ECMWF 10 m wind < 10 km/h AND total cloud
     cover < 25 % in at least 75 % of the 24 night-time bins.  All 24 weather
     bins must be available (partial nights are not classified).
  5. Network completeness.  A night is analysed only if >= 20 sensors of the
     regional network reported (the median of a handful of sensors is not a
     regional reference).
  6. Sensor-night completeness.  A sensor-night is kept only if >= 18 of its
     24 night bins are present.
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "raw")
EXT = os.path.join(DATA, "external")
PDATA = os.path.join(ROOT, "replication", "results", "data")   # this implementation's own processed tables
os.makedirs(PDATA, exist_ok=True)

WIND_MAX_KMH = 10.0       # study-specific threshold (ECMWF 10 m wind)
CLOUD_MAX_PCT = 25.0      # study-specific threshold (ECMWF total cloud cover)
MIN_FRAC_STEPS = 0.75     # share of night bins that must satisfy both
MIN_SENSORS_PER_NIGHT = 20
MIN_STEPS_PER_SENSOR_NIGHT = 18
NIGHT_START, NIGHT_END = 18, 6


def load_binned():
    """Half-hourly record with 30-min bins, full-network median and the weather series joined."""
    hh = pd.read_csv(os.path.join(RAW, "logger_readings.csv"), parse_dates=["datetime"])
    hh["time_bin"] = hh["datetime"].dt.round("30min")
    hh["hour_dec"] = hh["time_bin"].dt.hour + hh["time_bin"].dt.minute / 60.0
    hh["night_date"] = (hh["time_bin"] - pd.Timedelta(hours=NIGHT_END)).dt.strftime("%Y-%m-%d")
    hh["is_night"] = (hh["hour_dec"] >= NIGHT_START) | (hh["hour_dec"] < NIGHT_END)
    hh["hour"] = hh["time_bin"].dt.hour
    # one reading per sensor and bin (duplicates would arise only from clock resets)
    hh = hh.sort_values(["region", "sensor_id", "datetime"]).drop_duplicates(["region", "sensor_id", "time_bin"])
    med = (hh.groupby(["region", "time_bin"])
             .agg(network_median_temp=("temp_c", "median"), n_sensors_active=("sensor_id", "nunique"))
             .reset_index())
    hh = hh.merge(med, on=["region", "time_bin"])
    hh["dT_network"] = hh["temp_c"] - hh["network_median_temp"]

    # The hourly weather values were attached to the logger time stamps by nearest
    # hour, the half-hour stamp being assigned to the following hour (a :31
    # reading carries the value of the next full hour).  The bin value is the
    # one carried by the :01 / :31 stamps of the main logger group.
    wx = pd.read_csv(os.path.join(EXT, "weather_hourly_series.csv"), parse_dates=["datetime"])
    wx["time_bin"] = wx["datetime"].dt.round("30min")
    wx["pref"] = (wx["datetime"].dt.minute % 30 != 1).astype(int)     # :01/:31 rows first
    wx = (wx.sort_values(["region", "time_bin", "pref"]).drop_duplicates(["region", "time_bin"])
              .rename(columns={"wind": "wind_kmh", "cloud": "cloud_pct"})
              [["region", "time_bin", "wind_kmh", "cloud_pct"]])
    hh = hh.merge(wx, on=["region", "time_bin"], how="left")
    return hh, wx


def night_inventory(hh, wx, wind_max=WIND_MAX_KMH, cloud_max=CLOUD_MAX_PCT, min_frac=MIN_FRAC_STEPS,
                    min_sensors=MIN_SENSORS_PER_NIGHT):
    """Classify every diurnal date of each region."""
    e = wx.copy()
    e["hour_dec"] = e.time_bin.dt.hour + e.time_bin.dt.minute / 60.0
    e["night_date"] = (e.time_bin - pd.Timedelta(hours=NIGHT_END)).dt.strftime("%Y-%m-%d")
    e["is_night"] = (e.hour_dec >= NIGHT_START) | (e.hour_dec < NIGHT_END)
    en = e[e.is_night].copy()
    en["ok"] = (en.wind_kmh < wind_max) & (en.cloud_pct < cloud_max)
    inv = (en.groupby(["region", "night_date"])
             .agg(n_weather_bins=("ok", "size"), frac_calm_clear=("ok", "mean"),
                  wind_night=("wind_kmh", "mean"), cloud_night=("cloud_pct", "mean"),
                  frac_calm=("wind_kmh", lambda x: np.mean(x < wind_max)),
                  frac_clear=("cloud_pct", lambda x: np.mean(x < cloud_max)))
             .reset_index())
    sens = (hh[hh.is_night].groupby(["region", "night_date"])
              .agg(n_sensors=("sensor_id", "nunique"), n_bins=("time_bin", "nunique"),
                   T_median_night=("network_median_temp", "mean"))
              .reset_index())
    inv = inv.merge(sens, on=["region", "night_date"], how="outer")
    inv["n_sensors"] = inv.n_sensors.fillna(0).astype(int)
    inv["meets_weather"] = (inv.frac_calm_clear >= min_frac) & (inv.n_weather_bins == 24)
    inv["meets_network"] = inv.n_sensors >= min_sensors
    inv["selected"] = inv.meets_weather & inv.meets_network
    reason = np.where(inv.selected, "selected",
             np.where(~inv.meets_weather & (inv.n_weather_bins < 24), "partial night (deployment/retrieval)",
             np.where(~inv.meets_weather, "not calm and clear",
                      "fewer than %d sensors reporting" % min_sensors)))
    inv["decision"] = reason
    return inv.sort_values(["region", "night_date"]).reset_index(drop=True)


def sensor_nights(hh, inv, min_steps=MIN_STEPS_PER_SENSOR_NIGHT):
    keep = inv.loc[inv.selected, ["region", "night_date"]]
    night = hh[hh.is_night].merge(keep, on=["region", "night_date"])
    sn = (night.groupby(["region", "sensor_id", "night_date"])
          .agg(dT_night=("dT_network", "mean"),
               dT_night_median=("dT_network", "median"),
               temp_night=("temp_c", "mean"),
               ref_night=("network_median_temp", "mean"),
               wind_night=("wind_kmh", "mean"),
               cloud_night=("cloud_pct", "mean"),
               n_steps=("temp_c", "size"))
          .reset_index())
    sn = sn[sn.n_steps >= min_steps].copy()
    sn["night_no"] = sn.groupby("region")["night_date"].rank(method="dense").astype(int)
    sn["night_id"] = sn.region.str[:3] + "_" + sn.night_no.astype(str)
    return sn


def main():
    hh, wx = load_binned()
    inv = night_inventory(hh, wx)
    inv.to_csv(os.path.join(PDATA, "night_inventory.csv"), index=False)
    print(inv[inv.selected | (inv.frac_calm_clear >= MIN_FRAC_STEPS)].to_string(index=False))

    keep = inv.loc[inv.selected, ["region", "night_date"]]
    cyc = hh.merge(keep, on=["region", "night_date"])
    cols = ["region", "sensor_id", "datetime", "time_bin", "night_date", "hour", "hour_dec", "is_night", "temp_c",
            "network_median_temp", "n_sensors_active", "dT_network", "wind_kmh", "cloud_pct"]
    cyc["is_calm"] = cyc.wind_kmh < WIND_MAX_KMH
    cyc["is_clear"] = cyc.cloud_pct < CLOUD_MAX_PCT
    cyc.sort_values(["region", "sensor_id", "datetime"])[cols + ["is_calm", "is_clear"]].to_csv(
        os.path.join(PDATA, "halfhourly_calm_clear.csv"), index=False)

    sn = sensor_nights(hh, inv)
    coords = pd.read_csv(os.path.join(RAW, "sites.csv"), usecols=["sensor_id", "latitude", "longitude"])
    sn = sn.merge(coords, on="sensor_id", how="left").sort_values(["region", "sensor_id", "night_date"])
    sn.to_csv(os.path.join(PDATA, "sensor_nights.csv"), index=False)
    print(sn.groupby("region").agg(sensors=("sensor_id", "nunique"), nights=("night_date", "nunique"),
                                   rows=("sensor_id", "size")))
    print(sn.groupby(["region", "night_date"]).sensor_id.nunique())


if __name__ == "__main__":
    main()
