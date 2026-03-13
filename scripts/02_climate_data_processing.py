"""
==============================================================================
Script 02: ERA5 Climate Data Download and Processing
==============================================================================
Downloads ERA5 reanalysis data (daily Tmax, Tmin, total precipitation) for
Australia, 2013-2023, from the Copernicus Climate Data Store (CDS).

SETUP REQUIRED (one-time):
  1. Register for an ECMWF account at https://cds.climate.copernicus.eu/
  2. Log in and go to your profile to find your Personal Access Token
  3. Accept the ERA5 dataset licence on the dataset page
  4. Create ~/.cdsapirc with:
       url: https://cds.climate.copernicus.eu/api
       key: YOUR-PERSONAL-ACCESS-TOKEN
  5. Install dependencies:
       pip install cdsapi xarray netcdf4 geopandas regionmask numpy pandas

Outputs:
  data/raw/climate/era5_daily_tmax_YYYY.nc  (one file per year)
  data/raw/climate/era5_daily_tmin_YYYY.nc
  data/raw/climate/era5_daily_precip_YYYY.nc

Runtime: Several hours (CDS queue times vary). Downloads run year-by-year.
==============================================================================
"""

import os
import sys
from pathlib import Path

import cdsapi

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
CLIMATE_DIR = PROJECT_DIR / "data" / "raw" / "climate"
CLIMATE_DIR.mkdir(parents=True, exist_ok=True)

# Study period
YEARS = list(range(2013, 2024))  # 2013-2023

# Australia bounding box [North, West, South, East]
AUSTRALIA_BBOX = [-10, 112, -44, 154]

# All days and months for requests
ALL_MONTHS = [f"{m:02d}" for m in range(1, 13)]
ALL_DAYS = [f"{d:02d}" for d in range(1, 32)]


def download_era5_daily_stats(client, variable, statistic, year, output_path):
    """Download ERA5 daily statistics for one variable, one year."""
    if output_path.exists():
        print(f"  Already exists: {output_path.name}, skipping.")
        return

    print(f"  Requesting: {variable} ({statistic}) for {year}...")
    print(f"  (CDS may queue this request — check progress at cds.climate.copernicus.eu)")

    try:
        client.retrieve(
            "derived-era5-single-levels-daily-statistics",
            {
                "product_type": "reanalysis",
                "variable": variable,
                "year": str(year),
                "month": ALL_MONTHS,
                "day": ALL_DAYS,
                "daily_statistic": statistic,
                "time_zone": "utc+10:00",  # AEST for Australian daily max/min
                "frequency": "1_hourly",
                "area": AUSTRALIA_BBOX,
            },
            str(output_path),
        )
        print(f"  -> Saved: {output_path.name}")
    except Exception as e:
        print(f"  ERROR downloading {variable} {year}: {e}")
        # Remove partial file if it exists
        if output_path.exists():
            output_path.unlink()
        raise


def main():
    # Test CDS connection
    print("Connecting to Copernicus CDS...")
    try:
        client = cdsapi.Client()
    except Exception as e:
        print(f"\nERROR: Could not connect to CDS API.\n{e}")
        print("\nSetup instructions:")
        print("  1. Register at https://cds.climate.copernicus.eu/")
        print("  2. Create ~/.cdsapirc with your API key")
        print("  See script header for full instructions.")
        sys.exit(1)

    print(f"Connected. Downloading ERA5 data for {YEARS[0]}-{YEARS[-1]}.\n")

    # Download configuration: (variable, statistic, filename_prefix)
    downloads = [
        ("2m_temperature", "daily_maximum", "era5_daily_tmax"),
        ("2m_temperature", "daily_minimum", "era5_daily_tmin"),
        ("total_precipitation", "daily_sum", "era5_daily_precip"),
    ]

    for year in YEARS:
        print(f"\n{'='*60}")
        print(f"Year: {year}")
        print(f"{'='*60}")

        for variable, statistic, prefix in downloads:
            output_path = CLIMATE_DIR / f"{prefix}_{year}.nc"
            download_era5_daily_stats(client, variable, statistic, year, output_path)

    # Verify downloads
    print(f"\n{'='*60}")
    print("Download Summary")
    print(f"{'='*60}")
    nc_files = sorted(CLIMATE_DIR.glob("era5_*.nc"))
    expected = len(YEARS) * 3  # 3 variables per year
    print(f"Files downloaded: {len(nc_files)} / {expected} expected")
    for f in nc_files:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}  ({size_mb:.1f} MB)")

    if len(nc_files) < expected:
        missing = expected - len(nc_files)
        print(f"\nWARNING: {missing} files missing. Re-run to retry failed downloads.")
    else:
        print("\nAll ERA5 downloads complete.")


if __name__ == "__main__":
    main()
