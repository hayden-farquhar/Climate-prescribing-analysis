"""
==============================================================================
Script 02c: ERA5 Download via Google Cloud (ARCO-ERA5)
==============================================================================
Reads ERA5 directly from Google's public ARCO-ERA5 Zarr store — no CDS API
key, no queue, no server-side computation.

Processes day-by-day sequentially for reliability — each .compute() fetches
only 24 hourly chunks per variable, which avoids connection stalls. Progress
is printed every 7 days with ETA.

Skips years where output files already exist.

Requirements:
  pip install gcsfs zarr xarray dask netcdf4 numpy

Outputs (same format as scripts 02/02b):
  data/raw/climate/era5_daily_tmax_YYYY.nc
  data/raw/climate/era5_daily_tmin_YYYY.nc
  data/raw/climate/era5_daily_precip_YYYY.nc
==============================================================================
"""

import sys
import time as time_module
import calendar
from datetime import date, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

try:
    import gcsfs
except ImportError:
    print("Missing dependency. Run:  pip install gcsfs zarr dask")
    sys.exit(1)

try:
    import dask
except ImportError:
    print("Missing dependency. Run:  pip install dask")
    sys.exit(1)

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
CLIMATE_DIR = PROJECT_DIR / "data" / "raw" / "climate"
CLIMATE_DIR.mkdir(parents=True, exist_ok=True)

# Study period
YEARS = list(range(2013, 2024))

# Australia bounding box (ERA5 latitude is descending: 90 to -90)
LAT_SLICE = slice(-10, -44)
LON_SLICE = slice(112, 154)

# ARCO-ERA5 public Zarr store on Google Cloud
ARCO_STORE = "gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"

# Parallel threads for each .compute() call (24 chunks / 8 threads = 3 batches)
DASK_THREADS = 8

# How many days between refreshing the GCS connection
RECONNECT_EVERY = 14


def open_arco_era5():
    """Open the ARCO-ERA5 dataset from Google Cloud (public, no auth)."""
    fs = gcsfs.GCSFileSystem(token="anon")
    mapper = fs.get_mapper(ARCO_STORE)
    ds = xr.open_zarr(mapper, consolidated=True)
    ds = ds[["2m_temperature", "total_precipitation"]]
    return ds


def process_day(ds, d, need_temp=True, need_precip=True):
    """Download and compute daily stats for a single AEST day.

    AEST = UTC+10, so AEST day = UTC 14:00 previous day to UTC 13:00 current day.
    Each variable = 24 hourly chunks — small enough to never stall.
    """
    prev = d - timedelta(days=1)
    utc_start = f"{prev.isoformat()}T14:00"
    utc_end = f"{d.isoformat()}T13:00"

    daily_tmax = daily_tmin = daily_precip = None

    with dask.config.set(scheduler="threads", num_workers=DASK_THREADS):
        if need_temp:
            temp = ds["2m_temperature"].sel(
                time=slice(utc_start, utc_end),
                latitude=LAT_SLICE,
                longitude=LON_SLICE,
            ).compute()
            daily_tmax = temp.max(dim="time")
            daily_tmin = temp.min(dim="time")
            del temp

        if need_precip:
            precip = ds["total_precipitation"].sel(
                time=slice(utc_start, utc_end),
                latitude=LAT_SLICE,
                longitude=LON_SLICE,
            ).compute()
            daily_precip = precip.sum(dim="time")
            del precip

    return daily_tmax, daily_tmin, daily_precip


def process_year(year):
    """Process all days in a year and save yearly files."""
    tmax_path = CLIMATE_DIR / f"era5_daily_tmax_{year}.nc"
    tmin_path = CLIMATE_DIR / f"era5_daily_tmin_{year}.nc"
    precip_path = CLIMATE_DIR / f"era5_daily_precip_{year}.nc"

    need_temp = not (tmax_path.exists() and tmin_path.exists())
    need_precip = not precip_path.exists()

    if not need_temp and not need_precip:
        print(f"  All files exist, skipping.")
        return

    vars_label = []
    if need_temp:
        vars_label.append("temp")
    if need_precip:
        vars_label.append("precip")
    print(f"  Downloading: {' + '.join(vars_label)}")

    # Generate all days in this year
    days = []
    d = date(year, 1, 1)
    end = date(year, 12, 31)
    while d <= end:
        days.append(d)
        d += timedelta(days=1)

    total_days = len(days)
    year_t0 = time_module.time()
    tmax_list = []
    tmin_list = []
    precip_list = []

    ds = open_arco_era5()
    days_since_reconnect = 0

    for i, d in enumerate(days):
        # Refresh connection periodically to avoid stale sockets
        if days_since_reconnect >= RECONNECT_EVERY:
            ds.close()
            ds = open_arco_era5()
            days_since_reconnect = 0

        t0 = time_module.time()
        try:
            daily_tmax, daily_tmin, daily_precip = process_day(
                ds, d, need_temp=need_temp, need_precip=need_precip
            )
        except Exception as e:
            # Reconnect and retry once
            print(f"    {d}: ERROR ({e}), reconnecting...", flush=True)
            ds.close()
            ds = open_arco_era5()
            days_since_reconnect = 0
            daily_tmax, daily_tmin, daily_precip = process_day(
                ds, d, need_temp=need_temp, need_precip=need_precip
            )

        # Assign date coordinate for later concatenation
        time_coord = np.datetime64(d.isoformat())
        if need_temp:
            tmax_list.append(daily_tmax.expand_dims(time=[time_coord]))
            tmin_list.append(daily_tmin.expand_dims(time=[time_coord]))
        if need_precip:
            precip_list.append(daily_precip.expand_dims(time=[time_coord]))

        days_since_reconnect += 1
        elapsed = time_module.time() - t0
        done = i + 1
        so_far = time_module.time() - year_t0
        eta = (so_far / done) * (total_days - done)

        # Print progress every 7 days or on day 1/last
        if done == 1 or done % 7 == 0 or done == total_days:
            print(
                f"    {d}: {elapsed:.0f}s  ({done}/{total_days} days, ETA: {eta / 60:.0f} min)",
                flush=True,
            )

    ds.close()

    # Concatenate and save
    if need_temp:
        yearly_tmax = xr.concat(tmax_list, dim="time")
        yearly_tmin = xr.concat(tmin_list, dim="time")
        yearly_tmax.to_dataset(name="t2m").to_netcdf(tmax_path)
        yearly_tmin.to_dataset(name="t2m").to_netcdf(tmin_path)
        print(f"  -> Saved: {tmax_path.name}, {tmin_path.name}")

    if need_precip:
        yearly_precip = xr.concat(precip_list, dim="time")
        yearly_precip.to_dataset(name="tp").to_netcdf(precip_path)
        print(f"  -> Saved: {precip_path.name}")

    total = time_module.time() - year_t0
    print(f"  Year {year} complete in {total / 60:.1f} min")


def main():
    print("ERA5 Cloud Download (ARCO-ERA5 on Google Cloud)")
    print("=" * 60)
    print("No CDS API — reads directly from public cloud storage.")
    print(f"Day-by-day sequential, {DASK_THREADS} threads per day.\n")

    # Show existing files
    existing = sorted(CLIMATE_DIR.glob("era5_*.nc"))
    print(f"Existing ERA5 files: {len(existing)}")
    for f in existing:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}  ({size_mb:.1f} MB)")

    # Determine what's needed
    needed = []
    for year in YEARS:
        tmax = (CLIMATE_DIR / f"era5_daily_tmax_{year}.nc").exists()
        tmin = (CLIMATE_DIR / f"era5_daily_tmin_{year}.nc").exists()
        precip = (CLIMATE_DIR / f"era5_daily_precip_{year}.nc").exists()
        if not (tmax and tmin and precip):
            needed.append(year)

    if not needed:
        print("\nAll files already exist. Nothing to download.")
        return

    print(f"\nYears to process: {needed}")
    days_total = sum(366 if calendar.isleap(y) else 365 for y in needed)
    print(f"Total days to download: {days_total}\n")

    # Verify we can connect
    print("Testing connection...")
    ds = open_arco_era5()
    print("  Connected OK.\n")
    ds.close()

    # Process each year
    for year in YEARS:
        print(f"\n{'=' * 60}")
        print(f"Year: {year}")
        print(f"{'=' * 60}")
        process_year(year)

    # Summary
    print(f"\n{'=' * 60}")
    print("Download Summary")
    print(f"{'=' * 60}")
    nc_files = sorted(CLIMATE_DIR.glob("era5_*.nc"))
    expected = len(YEARS) * 3
    print(f"Files: {len(nc_files)} / {expected}")
    for f in nc_files:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}  ({size_mb:.1f} MB)")

    if len(nc_files) >= expected:
        print("\nAll ERA5 downloads complete.")
    else:
        print(f"\n{expected - len(nc_files)} files remaining. Re-run to retry.")


if __name__ == "__main__":
    main()
