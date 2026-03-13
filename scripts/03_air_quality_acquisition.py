"""
==============================================================================
Script 03: OpenAQ PM2.5 Air Quality Data Download
==============================================================================
Downloads daily PM2.5 data from OpenAQ for all Australian monitoring stations,
covering 2018-2021 (pre/during/post Black Summer bushfire season).

Uses a hybrid approach:
  1. API to discover Australian PM2.5 station locations and sensor IDs
  2. API daily endpoint to download aggregated daily data per sensor

SETUP REQUIRED (one-time):
  1. Register (free) at https://explore.openaq.org/register
  2. Get your API key from your profile
  3. Set environment variable: export OPENAQ_API_KEY="your-key-here"
     Or create a .env file in the project root with OPENAQ_API_KEY=your-key
  4. Install dependencies:
       pip install requests pandas python-dotenv

Outputs:
  data/raw/air_quality/openaq_stations_australia.csv  (station metadata)
  data/raw/air_quality/openaq_pm25_daily.csv          (daily PM2.5 values)

Runtime: ~15-30 minutes (rate-limited to 60 requests/minute)
==============================================================================
"""

import os
import sys
import json
import time
from pathlib import Path

import pandas as pd
import requests

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
AQ_DIR = PROJECT_DIR / "data" / "raw" / "air_quality"
AQ_DIR.mkdir(parents=True, exist_ok=True)

# Study period for bushfire analysis (with buffer)
DATE_FROM = "2018-01-01"
DATE_TO = "2021-12-31"

# OpenAQ API v3 settings
BASE_URL = "https://api.openaq.org/v3"
PM25_PARAMETER_ID = 2
RATE_LIMIT_DELAY = 1.1  # seconds between requests (60/min limit)


def get_api_key():
    """Get OpenAQ API key from environment or .env file."""
    api_key = os.environ.get("OPENAQ_API_KEY")

    if not api_key:
        env_file = PROJECT_DIR / ".env"
        if env_file.exists():
            for line in env_file.read_text().splitlines():
                if line.startswith("OPENAQ_API_KEY="):
                    api_key = line.split("=", 1)[1].strip().strip("'\"")
                    break

    if not api_key:
        print("ERROR: OpenAQ API key not found.")
        print("\nSetup instructions:")
        print("  1. Register at https://explore.openaq.org/register")
        print("  2. Set: export OPENAQ_API_KEY='your-key-here'")
        print("  Or add to .env file in project root: OPENAQ_API_KEY=your-key")
        sys.exit(1)

    return api_key


def api_get(endpoint, params, headers):
    """Make a GET request with rate limiting and error handling."""
    time.sleep(RATE_LIMIT_DELAY)

    url = f"{BASE_URL}/{endpoint}"
    resp = requests.get(url, params=params, headers=headers, timeout=30)

    if resp.status_code == 429:
        # Rate limited — wait and retry
        retry_after = int(resp.headers.get("Retry-After", 60))
        print(f"  Rate limited. Waiting {retry_after}s...")
        time.sleep(retry_after)
        resp = requests.get(url, params=params, headers=headers, timeout=30)

    resp.raise_for_status()
    return resp.json()


def discover_stations(headers):
    """Find all Australian PM2.5 monitoring stations."""
    print("Discovering Australian PM2.5 stations...")

    stations = []
    page = 1

    while True:
        data = api_get(
            "locations",
            params={
                "iso": "AU",
                "parameters_id": PM25_PARAMETER_ID,
                "limit": 1000,
                "page": page,
            },
            headers=headers,
        )

        results = data.get("results", [])
        if not results:
            break

        for loc in results:
            # Extract PM2.5 sensor ID
            pm25_sensor_id = None
            for sensor in loc.get("sensors", []):
                param = sensor.get("parameter", {})
                if param.get("id") == PM25_PARAMETER_ID:
                    pm25_sensor_id = sensor["id"]
                    break

            if pm25_sensor_id is None:
                continue

            stations.append({
                "location_id": loc["id"],
                "location_name": loc.get("name", ""),
                "sensor_id": pm25_sensor_id,
                "latitude": loc.get("coordinates", {}).get("latitude"),
                "longitude": loc.get("coordinates", {}).get("longitude"),
                "datetime_first": loc.get("datetimeFirst"),
                "datetime_last": loc.get("datetimeLast"),
                "provider": loc.get("provider", {}).get("name", ""),
            })

        found = data.get("meta", {}).get("found", 0)
        print(f"  Page {page}: {len(results)} locations (total found: {found})")

        if len(results) < 1000:
            break
        page += 1

    print(f"  Total stations with PM2.5 data: {len(stations)}")
    return pd.DataFrame(stations)


def download_daily_pm25(sensor_id, location_name, headers):
    """Download daily PM2.5 for a single sensor over the study period."""
    all_days = []
    page = 1

    while True:
        try:
            data = api_get(
                f"sensors/{sensor_id}/days",
                params={
                    "date_from": DATE_FROM,
                    "date_to": DATE_TO,
                    "limit": 1000,
                    "page": page,
                },
                headers=headers,
            )
        except requests.exceptions.HTTPError as e:
            print(f"    Error for sensor {sensor_id}: {e}")
            break

        results = data.get("results", [])
        if not results:
            break

        for r in results:
            period = r.get("period", {})
            summary = r.get("summary", {})
            coverage = r.get("coverage", {})

            all_days.append({
                "sensor_id": sensor_id,
                "date": period.get("datetimeFrom", {}).get("local", "")[:10],
                "pm25_mean": r.get("value"),
                "pm25_min": summary.get("min"),
                "pm25_max": summary.get("max"),
                "pm25_sd": summary.get("sd"),
                "expected_hours": coverage.get("expectedCount"),
                "observed_hours": coverage.get("observedCount"),
            })

        if len(results) < 1000:
            break
        page += 1

    return all_days


def main():
    api_key = get_api_key()
    headers = {"X-API-Key": api_key}

    # Step 1: Discover stations
    stations_file = AQ_DIR / "openaq_stations_australia.csv"

    if stations_file.exists():
        print(f"Loading cached station list: {stations_file.name}")
        stations = pd.read_csv(stations_file)
    else:
        stations = discover_stations(headers)
        stations.to_csv(stations_file, index=False)
        print(f"  -> Saved: {stations_file.name}")

    print(f"\n{len(stations)} PM2.5 stations found across Australia.\n")

    # Step 2: Download daily PM2.5 for each station
    print(f"Downloading daily PM2.5 data ({DATE_FROM} to {DATE_TO})...")
    print(f"This will take ~{len(stations) * 2} minutes due to rate limits.\n")

    all_daily = []
    for i, row in stations.iterrows():
        name = row["location_name"]
        sid = row["sensor_id"]
        print(f"  [{i+1}/{len(stations)}] {name} (sensor {sid})...", end=" ")

        days = download_daily_pm25(sid, name, headers)
        all_daily.extend(days)
        print(f"{len(days)} days")

    # Combine and save
    daily_df = pd.DataFrame(all_daily)

    if len(daily_df) > 0:
        # Merge station metadata
        daily_df = daily_df.merge(
            stations[["sensor_id", "location_id", "location_name",
                       "latitude", "longitude", "provider"]],
            on="sensor_id",
            how="left",
        )

        # Sort and save
        daily_df = daily_df.sort_values(["location_name", "date"])
        output_file = AQ_DIR / "openaq_pm25_daily.csv"
        daily_df.to_csv(output_file, index=False)
        print(f"\n-> Saved: {output_file.name}")
        print(f"   {len(daily_df)} daily records across {daily_df['sensor_id'].nunique()} sensors")
        print(f"   Date range: {daily_df['date'].min()} to {daily_df['date'].max()}")
    else:
        print("\nWARNING: No daily PM2.5 data was retrieved.")
        print("Check station coverage and date range.")

    # Step 3: Coverage summary
    print(f"\n{'='*60}")
    print("Coverage Summary")
    print(f"{'='*60}")
    if len(daily_df) > 0:
        coverage = (
            daily_df.groupby("location_name")
            .agg(
                n_days=("date", "nunique"),
                date_min=("date", "min"),
                date_max=("date", "max"),
                mean_pm25=("pm25_mean", "mean"),
            )
            .sort_values("n_days", ascending=False)
        )
        print(coverage.to_string())
        print(f"\nStations with data during Black Summer (Oct 2019 - Feb 2020):")
        bs_data = daily_df[
            (daily_df["date"] >= "2019-10-01") & (daily_df["date"] <= "2020-02-28")
        ]
        bs_stations = bs_data["location_name"].nunique()
        print(f"  {bs_stations} stations with {len(bs_data)} daily records")

    print("\nScript 03 complete.")


if __name__ == "__main__":
    main()
