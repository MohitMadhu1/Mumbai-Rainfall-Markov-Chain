"""Download IMD 0.25 deg gridded daily rainfall (1985-2024) from IMD Pune.

Source: India Meteorological Department, https://imdpune.gov.in/
Dataset: Pai et al. (2014), Mausam 65(1), 1-18. Retrieved with imdlib.

Downloads one year per call and retries each year, because the connection can
drop for minutes at a time. Years already on disk are skipped.
"""

import os
import time

import imdlib as imd

START, END = 1985, 2024
RAW_DIR = os.path.join("data", "raw")
MIN_BYTES = 20_000_000   # a full year of the 0.25 deg grid is ~25 MB


def have(year):
    path = os.path.join(RAW_DIR, "rain", f"{year}.grd")
    return os.path.exists(path) and os.path.getsize(path) > MIN_BYTES


def main():
    os.makedirs(RAW_DIR, exist_ok=True)
    for year in range(START, END + 1):
        tries = 0
        while not have(year):
            tries += 1
            try:
                imd.get_data("rain", year, year, fn_format="yearwise", file_dir=RAW_DIR)
            except Exception as exc:  # network drop, server hiccup
                print(f"{year}: attempt {tries} failed ({exc}); retrying in 30 s")
                time.sleep(30)
            if tries >= 20:
                raise SystemExit(f"{year}: giving up after {tries} attempts")
        print(f"{year}: ok")


if __name__ == "__main__":
    main()
