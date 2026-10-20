"""Build data/mumbai_daily_rain.csv from the raw downloads.

- IMD 0.25 deg grid cell centred on 19.00 N, 72.75 E. It spans 18.875-19.125 N,
  72.625-72.875 E, which contains the Santacruz observatory (19.117 N, 72.85 E).
- NOAA GHCN-Daily station IN012070800 (Bombay/Santacruz) for the cross-check.
  GHCN stores PRCP in tenths of a millimetre.
"""

import imdlib as imd
import numpy as np
import pandas as pd

START, END = 1985, 2024
LAT, LON = 19.0, 72.75


def imd_series():
    ds = imd.open_data("rain", START, END, "yearwise", "data/raw").get_xarray()
    s = ds["rain"].sel(lat=LAT, lon=LON).to_series()
    s = s.where(s != -999.0)
    s.index = pd.to_datetime(s.index).normalize()
    return s.rename("rain_imd")


def gauge_series():
    g = pd.read_csv("data/ghcn_IN012070800.csv", usecols=["DATE", "PRCP"],
                    parse_dates=["DATE"])
    g = g[(g.DATE.dt.year >= START) & (g.DATE.dt.year <= END)]
    return (g.set_index("DATE")["PRCP"] / 10.0).rename("rain_gauge")


def main():
    df = pd.concat([imd_series(), gauge_series()], axis=1)
    df.index.name = "date"
    df = df.loc[f"{START}-01-01":f"{END}-12-31"]
    print("days:", len(df), "| IMD missing:", int(df.rain_imd.isna().sum()),
          "| gauge missing:", int(df.rain_gauge.isna().sum()))
    df.to_csv("data/mumbai_daily_rain.csv", float_format="%.1f")


if __name__ == "__main__":
    main()
