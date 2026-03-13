"""
==============================================================================
Script 04: Spatial Linkage and Panel Construction
==============================================================================
Aggregates ERA5 gridded climate data to SA4/LGA administrative areas,
computes heatwave metrics, interpolates PM2.5 to areas, merges with
PBS prescribing and socioeconomic data to produce analysis-ready panels.

Inputs:
  data/raw/climate/era5_daily_tmax_YYYY.nc    (ERA5 daily max temperature)
  data/raw/climate/era5_daily_tmin_YYYY.nc    (ERA5 daily min temperature)
  data/raw/climate/era5_daily_precip_YYYY.nc  (ERA5 daily precipitation)
  data/spatial/SA4_2021/SA4_2021_AUST_GDA2020.shp
  data/spatial/LGA_2021/LGA_2021_AUST_GDA2020.shp
  data/raw/air_quality/openaq_pm25_daily.csv
  data/reference/SEIFA_2021_LGA.xlsx
  data/reference/ERP_LGA_2001-2024.xlsx
  data/processed/pbs_weekly_sa4.csv
  data/processed/pbs_monthly_lga.csv

Outputs:
  data/processed/climate_weekly_sa4.csv       (SA4 x weekly climate)
  data/processed/climate_monthly_lga.csv      (LGA x monthly climate)
  data/processed/pm25_weekly_sa4.csv          (SA4 x weekly PM2.5)
  data/processed/panel_weekly_sa4.csv         (full SA4 panel for DLNM)
  data/processed/panel_monthly_lga.csv        (full LGA panel for sensitivity)

Dependencies:
  pip install xarray netcdf4 geopandas regionmask scipy pandas numpy openpyxl

Runtime: ~10-20 minutes depending on how many ERA5 years are available
==============================================================================
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import regionmask
import xarray as xr
from scipy.spatial import cKDTree

# Project paths
PROJECT_DIR = Path(__file__).resolve().parent.parent
CLIMATE_DIR = PROJECT_DIR / "data" / "raw" / "climate"
SPATIAL_DIR = PROJECT_DIR / "data" / "spatial"
AQ_DIR = PROJECT_DIR / "data" / "raw" / "air_quality"
REF_DIR = PROJECT_DIR / "data" / "reference"
PROC_DIR = PROJECT_DIR / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)

YEARS = list(range(2013, 2024))


# ==============================================================================
# PART 1: ERA5 Spatial Aggregation to SA4 and LGA
# ==============================================================================

def load_shapefile(name):
    """Load an ABS shapefile, filter non-geographic regions, fix geometries."""
    if name == "SA4":
        path = SPATIAL_DIR / "SA4_2021" / "SA4_2021_AUST_GDA2020.shp"
        code_col = "SA4_CODE21"
    elif name == "LGA":
        path = SPATIAL_DIR / "LGA_2021" / "LGA_2021_AUST_GDA2020.shp"
        code_col = "LGA_CODE21"
    else:
        raise ValueError(f"Unknown shapefile: {name}")

    gdf = gpd.read_file(path)
    n_raw = len(gdf)

    # Filter out non-geographic regions:
    # - Codes ending in 97 (Migratory/Offshore/Shipping)
    # - Codes ending in 99 (No usual address)
    # - ZZZ (Outside Australia), 901 (Other Territories)
    codes = gdf[code_col].astype(str)
    exclude_mask = (
        codes.str.endswith("97") |
        codes.str.endswith("99") |
        codes.str.contains("ZZZ", case=False) |
        codes.str.startswith("9")  # Other Territories (9xx)
    )
    gdf = gdf[~exclude_mask].copy()

    # Remove empty geometries
    gdf = gdf[~gdf.geometry.is_empty].copy()

    # Fix invalid geometries
    invalid = ~gdf.geometry.is_valid
    if invalid.any():
        gdf.loc[invalid, "geometry"] = gdf.loc[invalid, "geometry"].buffer(0)

    gdf = gdf.reset_index(drop=True)

    print(f"  Loaded {name}: {n_raw} raw -> {len(gdf)} geographic regions")
    print(f"  Columns: {list(gdf.columns)}")
    return gdf


def load_era5_year(variable_prefix, year):
    """Load one ERA5 NetCDF file for a given variable and year."""
    filepath = CLIMATE_DIR / f"{variable_prefix}_{year}.nc"
    if not filepath.exists():
        return None
    ds = xr.open_dataset(filepath)
    # Rename ERA5 coords to standard names expected by regionmask
    rename_map = {}
    if "longitude" in ds.coords:
        rename_map["longitude"] = "lon"
    if "latitude" in ds.coords:
        rename_map["latitude"] = "lat"
    if "valid_time" in ds.coords:
        rename_map["valid_time"] = "time"
    if rename_map:
        ds = ds.rename(rename_map)
    return ds


def aggregate_era5_to_regions(regions_gdf, code_col, name_col):
    """
    Area-weighted aggregation of ERA5 grid data to polygon regions.

    Uses regionmask to create a 2D mask mapping grid cells to regions,
    then computes latitude-weighted means for each region per day.

    Returns a DataFrame with columns:
      region_code, date, tmax_c, tmin_c, precip_mm
    """
    # Create regionmask from geodataframe
    print("  Creating region mask...")
    abbrevs = regions_gdf[code_col].astype(str).tolist()

    # Add an integer index column for regionmask numbering
    regions_gdf = regions_gdf.copy()
    regions_gdf["_rm_idx"] = range(len(regions_gdf))

    region_mask = regionmask.from_geopandas(
        regions_gdf,
        names=name_col,
        abbrevs=code_col,
        numbers="_rm_idx",
    )

    all_results = []

    for year in YEARS:
        # Load all three variables for this year
        ds_tmax = load_era5_year("era5_daily_tmax", year)
        ds_tmin = load_era5_year("era5_daily_tmin", year)
        ds_precip = load_era5_year("era5_daily_precip", year)

        if ds_tmax is None:
            print(f"  Skipping {year} — file not found")
            continue

        print(f"  Processing {year}...", end=" ")

        # Create 3D mask (regions x lat x lon) — True where grid cell is in region
        mask_3d = region_mask.mask_3D(ds_tmax)

        # Latitude weights for area-weighted averaging
        weights = np.cos(np.deg2rad(ds_tmax.lat))

        # Process tmax (Kelvin -> Celsius)
        tmax_k = ds_tmax["t2m"]
        weighted_tmax = tmax_k.weighted(mask_3d * weights)
        tmax_by_region = weighted_tmax.mean(dim=["lat", "lon"]) - 273.15

        # Process tmin if available
        if ds_tmin is not None:
            tmin_k = ds_tmin["t2m"]
            weighted_tmin = tmin_k.weighted(mask_3d * weights)
            tmin_by_region = weighted_tmin.mean(dim=["lat", "lon"]) - 273.15
        else:
            tmin_by_region = None

        # Process precipitation if available (m -> mm)
        if ds_precip is not None:
            precip_var = "tp" if "tp" in ds_precip else list(ds_precip.data_vars)[0]
            precip_m = ds_precip[precip_var]
            weighted_precip = precip_m.weighted(mask_3d * weights)
            precip_by_region = weighted_precip.mean(dim=["lat", "lon"]) * 1000
        else:
            precip_by_region = None

        # Convert to DataFrame
        tmax_df = tmax_by_region.to_dataframe(name="tmax_c").reset_index()
        tmax_df = tmax_df.rename(columns={"time": "date", "region": "region_idx"})

        if tmin_by_region is not None:
            tmin_df = tmin_by_region.to_dataframe(name="tmin_c").reset_index()
            tmin_df = tmin_df.rename(columns={"time": "date", "region": "region_idx"})
            tmax_df = tmax_df.merge(tmin_df[["date", "region_idx", "tmin_c"]],
                                     on=["date", "region_idx"], how="left")

        if precip_by_region is not None:
            precip_df = precip_by_region.to_dataframe(name="precip_mm").reset_index()
            precip_df = precip_df.rename(columns={"time": "date", "region": "region_idx"})
            tmax_df = tmax_df.merge(precip_df[["date", "region_idx", "precip_mm"]],
                                     on=["date", "region_idx"], how="left")

        # Map region index back to region code
        idx_to_code = dict(enumerate(abbrevs))
        tmax_df["region_code"] = tmax_df["region_idx"].map(idx_to_code)
        tmax_df = tmax_df.drop(columns=["region_idx"])

        all_results.append(tmax_df)

        n_days = len(ds_tmax.time)
        print(f"{n_days} days")

        # Close datasets to free memory
        ds_tmax.close()
        if ds_tmin is not None:
            ds_tmin.close()
        if ds_precip is not None:
            ds_precip.close()

    if not all_results:
        print("  WARNING: No ERA5 data found!")
        return pd.DataFrame()

    result = pd.concat(all_results, ignore_index=True)
    result["date"] = pd.to_datetime(result["date"])
    return result


# ==============================================================================
# PART 2: Heatwave Metrics
# ==============================================================================

def compute_heatwave_metrics(daily_climate, region_col="region_code"):
    """
    Compute heatwave indicators from daily temperature data.

    For each region, computes:
    - Local temperature percentiles (90th, 95th, 99th of tmax distribution)
    - Heatwave flag: tmax exceeds local 95th percentile
    - Heatwave day: part of >=3 consecutive days above 95th percentile

    Returns the input DataFrame with additional columns.
    """
    print("  Computing local temperature percentiles...")
    percentiles = (
        daily_climate.groupby(region_col)["tmax_c"]
        .agg(
            tmax_p90=lambda x: np.nanpercentile(x, 90),
            tmax_p95=lambda x: np.nanpercentile(x, 95),
            tmax_p99=lambda x: np.nanpercentile(x, 99),
        )
        .reset_index()
    )

    daily_climate = daily_climate.merge(percentiles, on=region_col, how="left")

    # Flag days exceeding percentile thresholds
    daily_climate["above_p90"] = daily_climate["tmax_c"] > daily_climate["tmax_p90"]
    daily_climate["above_p95"] = daily_climate["tmax_c"] > daily_climate["tmax_p95"]
    daily_climate["above_p99"] = daily_climate["tmax_c"] > daily_climate["tmax_p99"]

    # Heatwave: >=3 consecutive days above 95th percentile
    print("  Identifying heatwave events (>=3 consecutive days above p95)...")
    daily_climate = daily_climate.sort_values([region_col, "date"])

    def mark_heatwaves(group):
        above = group["above_p95"].values
        hw = np.zeros(len(above), dtype=bool)
        # Find runs of consecutive True values
        run_start = None
        for i in range(len(above)):
            if above[i]:
                if run_start is None:
                    run_start = i
            else:
                if run_start is not None and (i - run_start) >= 3:
                    hw[run_start:i] = True
                run_start = None
        # Handle run at end
        if run_start is not None and (len(above) - run_start) >= 3:
            hw[run_start:] = True
        group["heatwave_day"] = hw
        return group

    groups = []
    for _, group in daily_climate.groupby(region_col):
        groups.append(mark_heatwaves(group))
    daily_climate = pd.concat(groups, ignore_index=True)

    n_hw = daily_climate["heatwave_day"].sum()
    print(f"  Total heatwave-days across all regions: {n_hw:,}")

    return daily_climate


# ==============================================================================
# PART 3: Aggregate Daily Climate to Weekly and Monthly
# ==============================================================================

def aggregate_to_weekly(daily_climate, region_col="region_code"):
    """Aggregate daily climate to weekly (Monday-start weeks)."""
    df = daily_climate.copy()
    df["week_start"] = df["date"] - pd.to_timedelta(df["date"].dt.dayofweek, unit="D")

    weekly = (
        df.groupby([region_col, "week_start"])
        .agg(
            tmax_mean=("tmax_c", "mean"),
            tmax_max=("tmax_c", "max"),
            tmin_mean=("tmin_c", "mean"),
            precip_total=("precip_mm", "sum"),
            days_above_p90=("above_p90", "sum"),
            days_above_p95=("above_p95", "sum"),
            days_above_p99=("above_p99", "sum"),
            heatwave_days=("heatwave_day", "sum"),
            n_days=("date", "count"),
        )
        .reset_index()
    )
    return weekly


def aggregate_to_monthly(daily_climate, region_col="region_code"):
    """Aggregate daily climate to monthly."""
    df = daily_climate.copy()
    df["year"] = df["date"].dt.year
    df["month"] = df["date"].dt.month

    monthly = (
        df.groupby([region_col, "year", "month"])
        .agg(
            tmax_mean=("tmax_c", "mean"),
            tmax_max=("tmax_c", "max"),
            tmin_mean=("tmin_c", "mean"),
            precip_total=("precip_mm", "sum"),
            days_above_p90=("above_p90", "sum"),
            days_above_p95=("above_p95", "sum"),
            days_above_p99=("above_p99", "sum"),
            heatwave_days=("heatwave_day", "sum"),
            n_days=("date", "count"),
        )
        .reset_index()
    )
    monthly["date"] = pd.to_datetime(
        monthly[["year", "month"]].assign(day=1)
    )
    return monthly


# ==============================================================================
# PART 4: PM2.5 Interpolation to SA4/LGA (IDW)
# ==============================================================================

def interpolate_pm25_to_regions(regions_gdf, code_col):
    """
    Interpolate station-level PM2.5 to administrative regions using
    inverse distance weighting (IDW).

    Uses region centroids as target points and station locations as
    source points. Power parameter = 2.
    """
    pm25_file = AQ_DIR / "openaq_pm25_daily.csv"
    if not pm25_file.exists():
        print("  WARNING: PM2.5 data not found, skipping interpolation.")
        return pd.DataFrame()

    print("  Loading PM2.5 station data...")
    pm25 = pd.read_csv(pm25_file)
    pm25["date"] = pd.to_datetime(pm25["date"])

    # Get unique station locations
    stations = (
        pm25.groupby("sensor_id")
        .agg(lat=("latitude", "first"), lon=("longitude", "first"),
             name=("location_name", "first"))
        .reset_index()
    )
    stations = stations.dropna(subset=["lat", "lon"])
    print(f"  {len(stations)} stations with coordinates")

    # Region centroids (use representative_point for robustness)
    centroids = regions_gdf.copy()
    centroids["centroid"] = centroids.geometry.representative_point()
    centroids["cent_lon"] = centroids["centroid"].x
    centroids["cent_lat"] = centroids["centroid"].y

    # Build KD-tree of station locations
    station_coords = np.column_stack([stations["lat"].values, stations["lon"].values])
    tree = cKDTree(station_coords)

    # For each region, find nearest K stations (max 10, within 200km)
    K = 10
    max_dist_deg = 2.0  # ~200km
    idw_power = 2

    region_codes = centroids[code_col].astype(str).values
    region_lats = centroids["cent_lat"].values
    region_lons = centroids["cent_lon"].values

    # Pre-compute nearest stations for each region
    region_station_map = {}
    for i, code in enumerate(region_codes):
        dists, idxs = tree.query([region_lats[i], region_lons[i]], k=K)
        # Filter by max distance
        mask = dists < max_dist_deg
        if not np.any(mask):
            # Use nearest single station as fallback
            mask = np.array([True] + [False] * (K - 1))
        region_station_map[code] = {
            "sensor_ids": stations.iloc[idxs[mask]]["sensor_id"].values,
            "distances": dists[mask],
        }

    # IDW interpolation per date per region
    print("  Interpolating PM2.5 to regions (IDW)...")
    dates = pm25["date"].unique()

    # Pivot PM2.5 to wide format for fast lookup
    pm25_pivot = pm25.pivot_table(
        index="date", columns="sensor_id", values="pm25_mean"
    )

    results = []
    for code in region_codes:
        info = region_station_map[code]
        sensor_ids = info["sensor_ids"]
        dists = info["distances"]

        # Get PM2.5 values for these sensors
        available_sensors = [s for s in sensor_ids if s in pm25_pivot.columns]
        if not available_sensors:
            continue

        sensor_vals = pm25_pivot[available_sensors]
        sensor_dists = dists[: len(available_sensors)]

        # IDW weights (inverse distance squared)
        weights = 1.0 / np.maximum(sensor_dists, 0.01) ** idw_power
        weights = weights / weights.sum()

        # Weighted average (handles NaN by re-normalising weights)
        for date_val in sensor_vals.index:
            vals = sensor_vals.loc[date_val].values
            valid = ~np.isnan(vals)
            if not np.any(valid):
                continue
            w = weights[: len(vals)][valid]
            w = w / w.sum()
            pm25_idw = np.sum(vals[valid] * w)
            results.append({
                "region_code": code,
                "date": date_val,
                "pm25_idw": pm25_idw,
                "n_stations": int(valid.sum()),
            })

    result = pd.DataFrame(results)
    print(f"  IDW interpolation complete: {len(result):,} region-day records")
    return result


# ==============================================================================
# PART 5: Load Reference Data (SEIFA, Population, Remoteness)
# ==============================================================================

def load_seifa():
    """Load SEIFA 2021 IRSD at LGA level."""
    filepath = REF_DIR / "SEIFA_2021_LGA.xlsx"
    if not filepath.exists():
        print("  WARNING: SEIFA file not found")
        return pd.DataFrame()

    # SEIFA Excel files typically have a header section to skip
    # Try reading with different skip values
    for skip in [4, 5, 6, 3]:
        try:
            df = pd.read_excel(filepath, sheet_name=1, skiprows=skip,
                               engine="openpyxl")
            # Look for IRSD score column
            cols_lower = [str(c).lower() for c in df.columns]
            if any("irsd" in c or "score" in c or "disadvantage" in c
                   for c in cols_lower):
                break
        except Exception:
            continue

    print(f"  SEIFA columns: {list(df.columns)}")
    print(f"  SEIFA rows: {len(df)}")
    return df


def load_population():
    """Load ERP by LGA."""
    filepath = REF_DIR / "ERP_LGA_2001-2024.xlsx"
    if not filepath.exists():
        print("  WARNING: Population file not found")
        return pd.DataFrame()

    # Try reading — ABS files have complex headers
    for skip in [4, 5, 6, 7, 3]:
        try:
            df = pd.read_excel(filepath, sheet_name=0, skiprows=skip,
                               engine="openpyxl")
            if len(df.columns) > 5:
                break
        except Exception:
            continue

    print(f"  Population columns: {list(df.columns)[:10]}...")
    print(f"  Population rows: {len(df)}")
    return df


# ==============================================================================
# PART 6: SA4-level SEIFA (population-weighted from LGA)
# ==============================================================================

def compute_sa4_seifa(seifa_lga, sa4_gdf, lga_gdf):
    """
    Assign SEIFA IRSD quintiles to SA4 regions by spatial overlay.
    Uses the area-weighted median of LGA IRSD scores within each SA4.
    """
    # This will be refined after inspecting SEIFA column names
    # For now, return placeholder
    print("  SA4-level SEIFA computation — will be refined after data inspection")
    return pd.DataFrame()


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 60)
    print("Script 04: Spatial Linkage and Panel Construction")
    print("=" * 60)

    # Check for ERA5 files
    era5_files = sorted(CLIMATE_DIR.glob("era5_daily_tmax_*.nc"))
    print(f"\nERA5 files available: {len(era5_files)} / {len(YEARS)} expected")
    if len(era5_files) == 0:
        print("ERROR: No ERA5 files found. Run script 02 first.")
        sys.exit(1)
    if len(era5_files) < len(YEARS):
        available_years = [int(f.stem.split("_")[-1]) for f in era5_files]
        print(f"  Available years: {available_years}")
        print("  Proceeding with available data (re-run when download completes)")

    # ------------------------------------------------------------------
    # Load shapefiles
    # ------------------------------------------------------------------
    print("\n--- Loading shapefiles ---")
    sa4_gdf = load_shapefile("SA4")
    lga_gdf = load_shapefile("LGA")

    # Identify code and name columns
    sa4_code_col = [c for c in sa4_gdf.columns if "CODE" in c.upper()][0]
    sa4_name_col = [c for c in sa4_gdf.columns if "NAME" in c.upper()][0]
    lga_code_col = [c for c in lga_gdf.columns if "CODE" in c.upper()][0]
    lga_name_col = [c for c in lga_gdf.columns if "NAME" in c.upper()][0]

    print(f"  SA4: code={sa4_code_col}, name={sa4_name_col}")
    print(f"  LGA: code={lga_code_col}, name={lga_name_col}")

    # ------------------------------------------------------------------
    # ERA5 -> SA4 aggregation (primary)
    # ------------------------------------------------------------------
    print("\n--- ERA5 -> SA4 aggregation ---")
    daily_sa4 = aggregate_era5_to_regions(sa4_gdf, sa4_code_col, sa4_name_col)

    if len(daily_sa4) > 0:
        print(f"\n  Daily SA4 climate: {len(daily_sa4):,} rows")
        print(f"  Date range: {daily_sa4['date'].min()} to {daily_sa4['date'].max()}")
        print(f"  Regions: {daily_sa4['region_code'].nunique()}")

        # Compute heatwave metrics
        print("\n--- Heatwave metrics (SA4) ---")
        daily_sa4 = compute_heatwave_metrics(daily_sa4)

        # Aggregate to weekly
        print("\n--- Aggregating to weekly ---")
        weekly_sa4_climate = aggregate_to_weekly(daily_sa4)
        print(f"  Weekly SA4 climate: {len(weekly_sa4_climate):,} rows")

        weekly_sa4_climate.to_csv(PROC_DIR / "climate_weekly_sa4.csv", index=False)
        print(f"  -> Saved: climate_weekly_sa4.csv")

    # ------------------------------------------------------------------
    # ERA5 -> LGA aggregation (sensitivity)
    # ------------------------------------------------------------------
    print("\n--- ERA5 -> LGA aggregation ---")
    print("  (This takes longer — more regions)")
    daily_lga = aggregate_era5_to_regions(lga_gdf, lga_code_col, lga_name_col)

    if len(daily_lga) > 0:
        print(f"\n  Daily LGA climate: {len(daily_lga):,} rows")

        # Compute heatwave metrics
        print("\n--- Heatwave metrics (LGA) ---")
        daily_lga = compute_heatwave_metrics(daily_lga)

        # Aggregate to monthly
        print("\n--- Aggregating to monthly ---")
        monthly_lga_climate = aggregate_to_monthly(daily_lga)
        print(f"  Monthly LGA climate: {len(monthly_lga_climate):,} rows")

        monthly_lga_climate.to_csv(PROC_DIR / "climate_monthly_lga.csv", index=False)
        print(f"  -> Saved: climate_monthly_lga.csv")

    # ------------------------------------------------------------------
    # PM2.5 interpolation to SA4
    # ------------------------------------------------------------------
    print("\n--- PM2.5 interpolation to SA4 ---")
    pm25_sa4 = interpolate_pm25_to_regions(sa4_gdf, sa4_code_col)

    if len(pm25_sa4) > 0:
        # Aggregate PM2.5 to weekly for SA4
        pm25_sa4["date"] = pd.to_datetime(pm25_sa4["date"])
        pm25_sa4["week_start"] = (
            pm25_sa4["date"]
            - pd.to_timedelta(pm25_sa4["date"].dt.dayofweek, unit="D")
        )
        pm25_weekly_sa4 = (
            pm25_sa4.groupby(["region_code", "week_start"])
            .agg(pm25_mean=("pm25_idw", "mean"), pm25_max=("pm25_idw", "max"),
                 n_stations_mean=("n_stations", "mean"))
            .reset_index()
        )
        pm25_weekly_sa4.to_csv(PROC_DIR / "pm25_weekly_sa4.csv", index=False)
        print(f"  -> Saved: pm25_weekly_sa4.csv ({len(pm25_weekly_sa4):,} rows)")

    # ------------------------------------------------------------------
    # Load reference data
    # ------------------------------------------------------------------
    print("\n--- Loading reference data ---")
    seifa = load_seifa()
    population = load_population()

    # ------------------------------------------------------------------
    # Merge into analysis panels
    # ------------------------------------------------------------------
    print("\n--- Building SA4 weekly panel ---")
    pbs_weekly = PROC_DIR / "pbs_weekly_sa4.csv"
    if pbs_weekly.exists() and len(daily_sa4) > 0:
        pbs = pd.read_csv(pbs_weekly, parse_dates=["week_start"])
        pbs["sa4_code"] = pbs["sa4_code"].astype(str)

        # Ensure matching types for merge
        weekly_sa4_climate["region_code"] = weekly_sa4_climate["region_code"].astype(str)
        weekly_sa4_climate["week_start"] = pd.to_datetime(weekly_sa4_climate["week_start"])

        # Merge PBS with climate
        panel_sa4 = pbs.merge(
            weekly_sa4_climate,
            left_on=["sa4_code", "week_start"],
            right_on=["region_code", "week_start"],
            how="left",
        ).drop(columns=["region_code"], errors="ignore")

        # Merge PM2.5 if available
        if len(pm25_sa4) > 0:
            pm25_weekly_sa4["region_code"] = pm25_weekly_sa4["region_code"].astype(str)
            pm25_weekly_sa4["week_start"] = pd.to_datetime(pm25_weekly_sa4["week_start"])
            panel_sa4 = panel_sa4.merge(
                pm25_weekly_sa4,
                left_on=["sa4_code", "week_start"],
                right_on=["region_code", "week_start"],
                how="left",
            ).drop(columns=["region_code"], errors="ignore")

        panel_sa4.to_csv(PROC_DIR / "panel_weekly_sa4.csv", index=False)
        print(f"  -> Saved: panel_weekly_sa4.csv ({len(panel_sa4):,} rows)")

        # Coverage check
        n_total = len(panel_sa4)
        n_climate = panel_sa4["tmax_mean"].notna().sum()
        print(f"  Climate coverage: {n_climate:,}/{n_total:,} "
              f"({100*n_climate/n_total:.1f}%)")

    print("\n--- Building LGA monthly panel ---")
    pbs_monthly = PROC_DIR / "pbs_monthly_lga.csv"
    if pbs_monthly.exists() and len(daily_lga) > 0:
        pbs_m = pd.read_csv(pbs_monthly, parse_dates=["date"])
        pbs_m["lga_code"] = pbs_m["lga_code"].astype(str)

        monthly_lga_climate["region_code"] = monthly_lga_climate["region_code"].astype(str)
        monthly_lga_climate["date"] = pd.to_datetime(monthly_lga_climate["date"])

        panel_lga = pbs_m.merge(
            monthly_lga_climate,
            left_on=["lga_code", "date"],
            right_on=["region_code", "date"],
            how="left",
        ).drop(columns=["region_code", "year_y", "month_y"], errors="ignore")

        panel_lga.to_csv(PROC_DIR / "panel_monthly_lga.csv", index=False)
        print(f"  -> Saved: panel_monthly_lga.csv ({len(panel_lga):,} rows)")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Script 04 complete")
    print("=" * 60)
    print("\nOutputs in data/processed/:")
    for f in sorted(PROC_DIR.glob("*.csv")):
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  {f.name}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
