# ==============================================================================
# Script 01: PBS Data Acquisition and Import
# ==============================================================================
# Imports AIHW PBS prescribing data from two sources:
#
# PRIMARY: AIHW PHE-338 "Geography and time-specific health data for
#          environmental analysis" — SA4 x weekly, Jul 2002 – Jun 2022
#          Files: aihw-phe-338-pbs-prescriptions-data-table-1-period-{1,2}.csv
#
# SENSITIVITY: AIHW HWE-098 PBS Monthly Data — LGA x monthly, ATC Level 2
#              Files: AIHW-HWE-098-PBS-ATC2-prescriptions-{historical,monthly}-data.xlsx
#
# Outputs:
#   data/processed/pbs_weekly_sa4.csv    (primary: SA4 x weekly, 2013 – Jun 2022)
#   data/processed/pbs_monthly_lga.csv   (sensitivity: LGA x monthly, 2013 – 2023)
# ==============================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(here)

project_dir <- here::here()
raw_dir     <- file.path(project_dir, "data", "raw", "pbs")
proc_dir    <- file.path(project_dir, "data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PART A: Environment-Health Dataset (SA4 x weekly) — PRIMARY
# ==============================================================================

cat("=" |> strrep(70), "\n")
cat("PART A: Environment-Health PBS Data (SA4 x weekly)\n")
cat("=" |> strrep(70), "\n\n")

# Read both period files (SA4-level data)
p1_file <- file.path(raw_dir, "aihw-phe-338-pbs-prescriptions-data-table-1-period-1.csv")
p2_file <- file.path(raw_dir, "aihw-phe-338-pbs-prescriptions-data-table-1-period-2.csv")

if (!file.exists(p1_file) || !file.exists(p2_file)) {
  stop("Missing PHE-338 CSV files in data/raw/pbs/. See script header for download instructions.")
}

cat("Reading period 1 (Jul 2002 – Jun 2012)...\n")
p1 <- read.csv(p1_file, stringsAsFactors = FALSE)
cat("  ", nrow(p1), "rows\n")

cat("Reading period 2 (Jul 2012 – Jun 2022)...\n")
p2 <- read.csv(p2_file, stringsAsFactors = FALSE)
cat("  ", nrow(p2), "rows\n")

# Coerce Count and Crude.rate to character in both before binding
# (period 1 has character due to suppressed cells like "np")
p1$Count <- as.character(p1$Count)
p2$Count <- as.character(p2$Count)
p1$Crude.rate <- as.character(p1$Crude.rate)
p2$Crude.rate <- as.character(p2$Crude.rate)

# Combine and standardise
weekly_raw <- bind_rows(p1, p2)
rm(p1, p2)

cat("\nCombined:", nrow(weekly_raw), "rows\n")
cat("Columns:", paste(names(weekly_raw), collapse = ", "), "\n\n")

# Standardise column names
weekly <- weekly_raw |>
  rename(
    week_start = Week.start.date,
    jurisdiction = Jurisdiction,
    sa4_code = SA4.code,
    sa4_name = SA4.name,
    rx_group = Prescription.group,
    count = Count,
    crude_rate = Crude.rate
  ) |>
  mutate(
    week_start = dmy(week_start),
    count = as.numeric(count),       # NAs from suppressed cells
    crude_rate = as.numeric(crude_rate),
    suppressed = is.na(count)
  )

# Map prescription groups to analysis categories
# The environment-health dataset uses therapeutic groupings, not ATC codes:
#   Cardiovascular prescriptions (aggregate C)
#     - Anti-thrombotic agents (B01 technically, but grouped under CVD here)
#     - Blood pressure lowering medicines (C02/C03/C07/C08/C09)
#     - Lipid-modifying medicines (C10)
#     - Other cardiovascular medicines
#   Mental health prescriptions (aggregate N05+N06)
#     - Anxiolytics (N05)
#     - Antidepressants (N06)
#     - Other mental health
#   Respiratory prescriptions (aggregate R03)
#     - Respiratory - relievers
#     - Respiratory - preventers
#     - Other respiratory (COPD-specific treatments, other asthma)
#     - Oral corticosteroids (H02, not strictly R03)

group_mapping <- c(
  # Aggregates (use these for main analysis)
  "Cardiovascular prescriptions"  = "cvd_total",
  "Mental health prescriptions"   = "mh_total",
  "Respiratory prescriptions"     = "resp_total",
  # Cardiovascular subgroups
  "Anti-thrombotic agents"        = "cvd_antithrombotic",
  "Blood pressure lowering medicines" = "cvd_bp_lowering",
  "Lipid-modifying medicines"     = "cvd_lipid_modifying",
  "Other cardiovascular medicines" = "cvd_other",
  # Mental health subgroups
  "Anxiolytics"                   = "mh_anxiolytics",
  "Antidepressants"               = "mh_antidepressants",
  "Other mental health"           = "mh_other",
  # Respiratory subgroups
  "Respiratory - relievers"       = "resp_relievers",
  "Respiratory - preventers"      = "resp_preventers",
  "Other respiratory (COPD-specific treatments, other asthma)" = "resp_copd_other",
  "Oral corticosteroids"          = "resp_oral_corticosteroids"
)

weekly <- weekly |>
  mutate(analysis_group = group_mapping[rx_group])

# Check for any unmapped groups
unmapped <- weekly |> filter(is.na(analysis_group)) |> pull(rx_group) |> unique()
if (length(unmapped) > 0) {
  cat("WARNING: Unmapped prescription groups:", paste(unmapped, collapse = ", "), "\n")
}

# Filter to study period: Jan 2013 onwards (environment-health ends Jun 2022)
weekly <- weekly |>
  filter(week_start >= as.Date("2013-01-01"))

cat("After filtering to 2013+:", nrow(weekly), "rows\n")
cat("Date range:", as.character(min(weekly$week_start)),
    "to", as.character(max(weekly$week_start)), "\n")
cat("SA4 areas:", n_distinct(weekly$sa4_code), "\n")
cat("Prescription groups:", n_distinct(weekly$rx_group), "\n")

# Suppression summary
n_suppressed <- sum(weekly$suppressed, na.rm = TRUE)
pct_suppressed <- round(100 * n_suppressed / nrow(weekly), 1)
cat("Suppressed cells:", n_suppressed, sprintf("(%.1f%%)\n", pct_suppressed))

# Save
weekly_out <- file.path(proc_dir, "pbs_weekly_sa4.csv")
write.csv(weekly, weekly_out, row.names = FALSE)
cat("\n-> Saved:", weekly_out, "\n")

# Summary table: total prescriptions by group
cat("\nPrescription counts by group (2013 – Jun 2022):\n")
weekly |>
  filter(!suppressed) |>
  group_by(rx_group) |>
  summarise(
    total_scripts = sum(count, na.rm = TRUE),
    n_weeks = n_distinct(week_start),
    n_sa4 = n_distinct(sa4_code),
    .groups = "drop"
  ) |>
  arrange(desc(total_scripts)) |>
  print(n = 20)

# ==============================================================================
# PART B: PBS Monthly ATC Level 2 (LGA x monthly) — SENSITIVITY
# ==============================================================================
# Excel structure (both files identical):
#   Table 6: "PBS prescriptions dispensed and government spending per resident
#             by LGA, for all patients"
#   Row 1 (skip): header row with column labels
#   Columns: State | LGA code | LGA name | Type of script | Value | MMM-YYYY ...
#   "Type of script" = ATC Level 2 group name
#   "Value" = "Benefits" (spending per resident) or "Scripts" (scripts per resident)
#   Month columns are wide format: JAN-2015, FEB-2015, ...
#
# Historical file: May 2002 – Jan 2015 (Table 6 has same structure)
# Current file:    Jan 2015 – present
# ==============================================================================

cat("\n\n", "=" |> strrep(70), "\n")
cat("PART B: PBS Monthly Data (LGA x monthly) — SENSITIVITY\n")
cat("=" |> strrep(70), "\n\n")

hist_file <- file.path(raw_dir, "AIHW-HWE-098-PBS-ATC2-prescriptions-historical-data.xlsx")
curr_file <- file.path(raw_dir, "AIHW-HWE-098-PBS-ATC2-prescriptions-monthly-data.xlsx")

if (!file.exists(hist_file) || !file.exists(curr_file)) {
  cat("WARNING: HWE-098 Excel files not found. Skipping Part B.\n")
  cat("The sensitivity analysis (LGA x monthly) requires these files.\n")
  cat("Download from: https://www.aihw.gov.au/reports/medicines/pbs-monthly-data/data\n")
  cat("\nScript 01 complete.\n")
} else {

  # --------------------------------------------------------------------------
  # Helper: Read and reshape one HWE-098 Excel file (Table 6 = all ages, LGA)
  # --------------------------------------------------------------------------
  read_hwe098_table6 <- function(filepath) {
    fname <- basename(filepath)
    cat("Reading Table 6 from:", fname, "\n")
    cat("  (This file is large — may take a few minutes...)\n")

    # Read Table 6, skip first row (which is the sheet title)
    raw <- read_excel(filepath, sheet = "Table 6", skip = 1,
                      col_names = FALSE, .name_repair = "minimal")

    # First row after skip contains actual headers
    headers <- as.character(raw[1, ])
    headers[1] <- "state"
    headers[2] <- "lga_code"
    headers[3] <- "lga_name"
    headers[4] <- "atc_group"
    headers[5] <- "value_type"

    # Remove header row from data
    raw <- raw[-1, ]

    # Assign column names
    names(raw) <- headers

    cat("  Dimensions:", nrow(raw), "rows x", ncol(raw), "columns\n")

    # Month columns are everything after the first 5
    month_cols <- headers[6:length(headers)]
    cat("  Month range:", month_cols[1], "to", month_cols[length(month_cols)], "\n")

    # Print unique ATC groups (first time only)
    atc_groups <- unique(raw$atc_group)
    cat("  ATC groups found:", length(atc_groups), "\n")
    for (g in sort(atc_groups)) cat("    -", g, "\n")

    # Filter to scripts only (not benefits/spending)
    raw <- raw |>
      filter(tolower(value_type) == "scripts per resident")

    if (nrow(raw) == 0) {
      # Try alternate label
      raw <- read_excel(filepath, sheet = "Table 6", skip = 1,
                        col_names = FALSE, .name_repair = "minimal")
      raw <- raw[-1, ]
      names(raw) <- headers
      # Check what value_type labels exist
      cat("  Value types found:", paste(unique(raw$value_type), collapse = ", "), "\n")
      raw <- raw |> filter(grepl("scri", tolower(value_type)))
    }

    cat("  After filtering to scripts:", nrow(raw), "rows\n")

    # Pivot from wide to long
    cat("  Pivoting to long format...\n")
    long <- raw |>
      tidyr::pivot_longer(
        cols = all_of(month_cols),
        names_to = "month_str",
        values_to = "scripts_per_resident"
      ) |>
      mutate(
        scripts_per_resident = as.numeric(scripts_per_resident),
        # Parse "JAN-2015" format
        date = dmy(paste0("01-", month_str)),
        year = year(date),
        month = month(date)
      ) |>
      select(state, lga_code, lga_name, atc_group, date, year, month,
             scripts_per_resident)

    cat("  Done:", nrow(long), "rows in long format\n")
    return(long)
  }

  # --------------------------------------------------------------------------
  # Read both files
  # --------------------------------------------------------------------------
  cat("\n--- Historical file ---\n")
  hist_long <- read_hwe098_table6(hist_file)

  cat("\n--- Current file ---\n")
  curr_long <- read_hwe098_table6(curr_file)

  # Combine (drop overlapping Jan 2015 from historical if present)
  monthly <- bind_rows(
    hist_long |> filter(date < as.Date("2015-01-01")),
    curr_long
  )

  cat("\n--- Combined ---\n")
  cat("Total rows:", nrow(monthly), "\n")
  cat("Date range:", as.character(min(monthly$date, na.rm = TRUE)),
      "to", as.character(max(monthly$date, na.rm = TRUE)), "\n")
  cat("LGAs:", n_distinct(monthly$lga_code), "\n")

  # Filter to study period 2013-2023
  monthly <- monthly |>
    filter(year >= 2013, year <= 2023)

  cat("After filtering to 2013-2023:", nrow(monthly), "rows\n")

  # Save full dataset (all ATC groups) — we can filter to specific groups later
  monthly_out <- file.path(proc_dir, "pbs_monthly_lga.csv")
  write.csv(monthly, monthly_out, row.names = FALSE)
  cat("\n-> Saved:", monthly_out, "\n")

  # Summary
  cat("\nPrescription groups in LGA monthly data:\n")
  monthly |>
    group_by(atc_group) |>
    summarise(
      n_months = n_distinct(date),
      n_lga = n_distinct(lga_code),
      mean_scripts_per_resident = mean(scripts_per_resident, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(atc_group) |>
    print(n = 30)

  cat("\nScript 01 complete.\n")
}
