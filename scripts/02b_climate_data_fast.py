"""
==============================================================================
Script 02b: ERA5 Download — Reliable Sequential Approach
==============================================================================
Downloads hourly ERA5 data from CDS, one full-year request at a time.
CDS does server-side spatial subsetting (Australia only), so each file is
~200-400 MB rather than the ~3 GB per month from cloud sources.

Key design choices for reliability:
  - ONE request at a time (no concurrency — avoids connection storms)
  - Full-year requests (fewer queue waits than monthly)
  - 600s timeout (generous for large downloads)
  - Uses reanalysis-era5-single-levels (pre-computed, no server-side
    aggregation unlike the derived daily-statistics product in script 02)
  - Temperature downloaded once → produces both tmax AND tmin

Expected speed: ~10-20 min per request (queue + download).
  2 requests/year × 8 years = ~16 requests = ~3-5 hours total.

Skips years where output files already exist.

Outputs (same format as script 02):
  data/raw/climate/era5_daily_tmax_YYYY.nc
  data/raw/climate/era5_daily_tmin_YYYY.nc
  data/raw/climate/era5_daily_precip_YYYY.nc
==============================================================================
"""

import sys
import time as time_module
import tempfile
from pathlib import Path

import cdsapi
import xarray as xr
import numpy as np

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
CLIMATE_DIR = PROJECT_DIR / "data" / "raw" / "climate"
CLIMATE_DIR.mkdir(parents=True, exist_ok=True)

# Study period
YEARS = list(range(2013, 2024))

# Australia bounding box [North, West, South, East]
AUSTRALIA_BBOX = [-10, 112, -44, 154]

ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS = [f"{d:02d}" for d in range(1, 32)]
ALL_HOURS = [f"{h:02d}:00" for h in range(24)]

# Generous timeout — large files can take a while to download
CDS_TIMEOUT = 600


def download_year(client, variable, year, output_path):
    """Download a full year of hourly ERA5 data for Australia (single request)."""
    print(f"    Submitting CDS request...", flush=True)
    client.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "variable": variable,
            "year": str(year),
            "month": ALL_MONTHS,
            "day": ALL_DAYS,
            "time": ALL_HOURS,
            "area": AUSTRALIA_BBOX,
            "data_format": "netcdf",
        },
        str(output_path),
    )
    size_mb = output_path.stat().st_size / (1024 * 1024)
    print(f"    Downloaded: {output_path.name} ({size_mb:.0f} MB)", flush=True)


def compute_daily_from_hourly(hourly_path, tmax_path=None, tmin_path=None, precip_path=None):
    """Compute daily statistics from an hourly NetCDF file."""
    ds = xr.open_dataset(hourly_path)

    # Shift to AEST (UTC+10) for Australian daily boundaries
    time_coord = "valid_time" if "valid_time" in ds.coords else "time"
    ds = ds.rename({time_coord: "time"}) if time_coord != "time" else ds
    ds["time"] = ds["time"] + np.timedelta64(10, "h")

    var_name = list(ds.data_vars)[0]
    da = ds[var_name]

    if tmax_path and tmin_path:
        print(f"    Computing daily max/min...", flush=True)
        daily_max = da.resample(time="1D").max()
        daily_min = da.resample(time="1D").min()
        daily_max.to_dataset(name="t2m").to_netcdf(tmax_path)
        daily_min.to_dataset(name="t2m").to_netcdf(tmin_path)
        print(f"    -> {tmax_path.name}, {tmin_path.name}")

    if precip_path:
        print(f"    Computing daily sum...", flush=True)
        daily_sum = da.resample(time="1D").sum()
        daily_sum.to_dataset(name="tp").to_netcdf(precip_path)
        print(f"    -> {precip_path.name}")

    ds.close()


def process_year(client, year):
    """Download and process one year of ERA5 data."""
    tmax_path = CLIMATE_DIR / f"era5_daily_tmax_{year}.nc"
    tmin_path = CLIMATE_DIR / f"era5_daily_tmin_{year}.nc"
    precip_path = CLIMATE_DIR / f"era5_daily_precip_{year}.nc"

    need_temp = not (tmax_path.exists() and tmin_path.exists())
    need_precip = not precip_path.exists()

    if not need_temp and not need_precip:
        print(f"  All files exist, skipping.")
        return

    year_t0 = time_module.time()

    with tempfile.TemporaryDirectory(prefix=f"era5_{year}_") as tmp_dir:
        tmp = Path(tmp_dir)

        if need_temp:
            print(f"  [1/2] 2m_temperature for {year} (1 request → tmax + tmin)")
            hourly_path = tmp / f"hourly_t2m_{year}.nc"
            t0 = time_module.time()
            try:
                download_year(client, "2m_temperature", year, hourly_path)
                compute_daily_from_hourly(hourly_path, tmax_path=tmax_path, tmin_path=tmin_path)
                elapsed = time_module.time() - t0
                print(f"    Temperature done in {elapsed / 60:.1f} min")
            except Exception as e:
                print(f"    ERROR: {e}")
                # Clean up partial outputs
                for p in [tmax_path, tmin_path]:
                    if p.exists():
                        p.unlink()
            finally:
                if hourly_path.exists():
                    hourly_path.unlink()

        if need_precip:
            step = "2/2" if need_temp else "1/1"
            print(f"  [{step}] total_precipitation for {year}")
            hourly_path = tmp / f"hourly_tp_{year}.nc"
            t0 = time_module.time()
            try:
                download_year(client, "total_precipitation", year, hourly_path)
                compute_daily_from_hourly(hourly_path, precip_path=precip_path)
                elapsed = time_module.time() - t0
                print(f"    Precipitation done in {elapsed / 60:.1f} min")
            except Exception as e:
                print(f"    ERROR: {e}")
                if precip_path.exists():
                    precip_path.unlink()
            finally:
                if hourly_path.exists():
                    hourly_path.unlink()

    total = time_module.time() - year_t0
    print(f"  Year {year} total: {total / 60:.1f} min")


def main():
    print("ERA5 Download — Sequential (reanalysis-era5-single-levels)")
    print("=" * 60)
    print("One request at a time. CDS subsets to Australia server-side.\n")

    # Show existing files
    existing = sorted(CLIMATE_DIR.glob("era5_*.nc"))
    print(f"Existing ERA5 files: {len(existing)}")
    for f in existing:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}  ({size_mb:.1f} MB)")

    # Determine what's needed
    work = []
    for year in YEARS:
        tmax = (CLIMATE_DIR / f"era5_daily_tmax_{year}.nc").exists()
        tmin = (CLIMATE_DIR / f"era5_daily_tmin_{year}.nc").exists()
        precip = (CLIMATE_DIR / f"era5_daily_precip_{year}.nc").exists()
        n_requests = (0 if (tmax and tmin) else 1) + (0 if precip else 1)
        if n_requests > 0:
            work.append((year, n_requests))

    if not work:
        print("\nAll files already exist. Nothing to download.")
        return

    total_requests = sum(n for _, n in work)
    print(f"\nWork remaining:")
    for year, n in work:
        print(f"  {year}: {n} request(s)")
    print(f"Total: {total_requests} CDS requests (~{total_requests * 15:.0f} min estimated)\n")

    # Connect to CDS
    print("Connecting to CDS...")
    try:
        client = cdsapi.Client(timeout=CDS_TIMEOUT)
    except Exception as e:
        print(f"ERROR: Could not connect to CDS API.\n{e}")
        sys.exit(1)

    # Process each year sequentially
    for year in YEARS:
        print(f"\n{'=' * 60}")
        print(f"Year: {year}")
        print(f"{'=' * 60}")
        process_year(client, year)

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
