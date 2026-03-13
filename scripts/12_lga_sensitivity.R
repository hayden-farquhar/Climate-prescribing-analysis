# ==============================================================================
# Script 12: LGA Monthly Sensitivity Analysis
# ==============================================================================
# Reruns the primary DLNM analysis at LGA x monthly resolution (vs SA4 x
# weekly in script 05) to test whether findings are robust to spatial and
# temporal aggregation scale.
#
# Key differences from script 05:
#   - LGA geography (~540 LGAs vs ~70 SA4s)
#   - Monthly temporal resolution (lag 0-3 months vs 0-12 weeks)
#   - Uses scripts_per_resident as outcome (population-adjusted)
#
# Inputs:
#   data/processed/panel_monthly_lga.csv
#
# Outputs:
#   outputs/tables/lga_sensitivity_rr.csv
#   outputs/figures/lga_sensitivity_er_*.png
#   outputs/figures/lga_sensitivity_comparison.png
#   outputs/tables/sa4_vs_lga_comparison.csv
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate",
#                       "ggplot2", "patchwork"))
#
# NOTE: The LGA panel is ~940 MB. This script may take several minutes and
# require substantial memory (~4-6 GB). Close other applications if needed.
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 3  # Monthly lags (0-3 months)

cat("=" |> strrep(70), "\n")
cat("Script 12: LGA Monthly Sensitivity Analysis\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading LGA monthly panel (this may take a moment)...\n")
panel <- read.csv(file.path(data_dir, "panel_monthly_lga.csv"),
                  stringsAsFactors = FALSE)
panel$date     <- as.Date(panel$date)
panel$lga_code <- as.character(panel$lga_code)

cat("  Rows:", format(nrow(panel), big.mark = ","), "\n")
cat("  Date range:", as.character(min(panel$date)),
    "to", as.character(max(panel$date)), "\n")
cat("  LGA regions:", n_distinct(panel$lga_code), "\n")
cat("  ATC groups:", paste(sort(unique(panel$atc_group)), collapse = "; "), "\n")

# Check climate coverage
n_climate <- sum(!is.na(panel$tmax_mean))
cat("  Climate coverage:", format(n_climate, big.mark = ","), "/",
    format(nrow(panel), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_climate / nrow(panel)))

# --- Map ATC group names to analysis groups ---
# The LGA panel uses descriptive ATC names rather than coded groups
# Map to match the SA4 analysis categories
panel <- panel |>
  mutate(
    analysis_group = case_when(
      grepl("renin-angiotensin|antihypertensive|lipid|cardiac|beta|calcium|diuretic|vasoprotective|antithrombotic",
            tolower(atc_group)) ~ "cvd_total",
      grepl("psycholeptic|antipsychotic|anxiolytic|hypnotic|sedative",
            tolower(atc_group)) ~ "mh_anxiolytics",
      grepl("antidepressant|psychoanaleptic",
            tolower(atc_group)) ~ "mh_antidepressants",
      grepl("obstructive airway|adrenergic.*inhalant|glucocorticoid.*inhalant|anticholinergic.*inhalant",
            tolower(atc_group)) ~ "resp_total",
      TRUE ~ "other"
    )
  )

# Aggregate to analysis group level (sum counts across ATC subgroups within each group)
panel_agg <- panel |>
  filter(analysis_group != "other",
         !is.na(tmax_mean)) |>
  group_by(lga_code, lga_name, state, date, analysis_group,
           tmax_mean, tmax_max, tmin_mean, precip_total,
           days_above_p90, days_above_p95, days_above_p99, heatwave_days) |>
  summarise(
    scripts_per_resident = sum(scripts_per_resident, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Aggregated analysis groups:\n")
print(table(panel_agg$analysis_group))

# Add time variables
panel_agg <- panel_agg |>
  mutate(
    year       = year(date),
    month      = month(date),
    time_index = as.numeric(date - min(date)) / 30.44,  # months since start
    sin1 = sin(2 * pi * month / 12),
    cos1 = cos(2 * pi * month / 12),
    sin2 = sin(4 * pi * month / 12),
    cos2 = cos(4 * pi * month / 12)
  )


# ==============================================================================
# 2. FIT DLNM FOR EACH MEDICATION GROUP
# ==============================================================================

# Primary groups only (subgroups may have limited power at LGA level)
primary_groups <- c("cvd_total", "mh_antidepressants", "mh_anxiolytics", "resp_total")

group_labels <- c(
  cvd_total          = "Cardiovascular",
  mh_antidepressants = "Antidepressants",
  mh_anxiolytics     = "Anxiolytics/Psycholeptics",
  resp_total         = "Respiratory"
)

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models at LGA x monthly resolution\n")
cat("=" |> strrep(70), "\n")

# Due to the large number of LGAs, we use a random sample approach if needed
# to keep computation manageable. With ~540 LGAs, the fixed effects model
# may be very large.

results     <- list()
all_summaries <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel_agg |>
    filter(analysis_group == grp,
           !is.na(scripts_per_resident),
           scripts_per_resident > 0,
           !is.na(tmax_mean)) |>
    arrange(lga_code, date)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")
  cat("  LGAs:", n_distinct(df$lga_code), "\n")

  if (nrow(df) < 500) {
    cat("  Too few observations, skipping.\n")
    next
  }

  # If too many LGAs for fixed effects, sample a manageable number
  n_lga <- n_distinct(df$lga_code)
  if (n_lga > 200) {
    cat("  Sampling 200 LGAs for computational tractability...\n")
    set.seed(42)
    sampled_lgas <- sample(unique(df$lga_code), 200)
    df <- df |> filter(lga_code %in% sampled_lgas)
    cat("  After sampling:", format(nrow(df), big.mark = ","), "rows,",
        n_distinct(df$lga_code), "LGAs\n")
  }

  # Temperature distribution
  temp_range  <- range(df$tmax_mean, na.rm = TRUE)
  temp_median <- median(df$tmax_mean, na.rm = TRUE)

  # Cross-basis: monthly resolution, lag 0-3 months
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 2))
  )

  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))

  cat("  Fitting model...\n")

  # Use scripts_per_resident * 1000 as a quasi-count for Gaussian approx
  # (rate outcome, not count — use Gaussian on log-transformed rate)
  df$log_rate <- log(df$scripts_per_resident + 0.001)

  model <- glm(
    log_rate ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(lga_code),
    data = df,
    family = gaussian()
  )

  cat("  Model fitted.\n")

  temp_ref <- temp_median

  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 80),
    cen = temp_ref
  )

  # Key percentile predictions
  temp_p05 <- quantile(df$tmax_mean, 0.05, na.rm = TRUE)
  temp_p10 <- quantile(df$tmax_mean, 0.10, na.rm = TRUE)
  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  pred_pctl <- crosspred(cb, model,
                          at = c(temp_p05, temp_p10, temp_p90, temp_p95, temp_p99),
                          cen = temp_ref)

  # For Gaussian on log-rate, exponentiate to get RR-like interpretation
  rr_fit <- exp(pred_pctl$allfit)
  rr_low <- exp(pred_pctl$alllow)
  rr_high <- exp(pred_pctl$allhigh)

  rr_row <- data.frame(
    group    = grp,
    label    = group_labels[grp],
    resolution = "LGA monthly",
    temp_ref = round(temp_ref, 1),
    temp_p05 = round(temp_p05, 1),
    rr_p05    = round(rr_fit[1], 4),
    rr_p05_lo = round(rr_low[1], 4),
    rr_p05_hi = round(rr_high[1], 4),
    temp_p10 = round(temp_p10, 1),
    rr_p10    = round(rr_fit[2], 4),
    rr_p10_lo = round(rr_low[2], 4),
    rr_p10_hi = round(rr_high[2], 4),
    temp_p90 = round(temp_p90, 1),
    rr_p90    = round(rr_fit[3], 4),
    rr_p90_lo = round(rr_low[3], 4),
    rr_p90_hi = round(rr_high[3], 4),
    temp_p95 = round(temp_p95, 1),
    rr_p95    = round(rr_fit[4], 4),
    rr_p95_lo = round(rr_low[4], 4),
    rr_p95_hi = round(rr_high[4], 4),
    temp_p99 = round(temp_p99, 1),
    rr_p99    = round(rr_fit[5], 4),
    rr_p99_lo = round(rr_low[5], 4),
    rr_p99_hi = round(rr_high[5], 4),
    n_obs = nrow(df),
    n_lga = n_distinct(df$lga_code),
    stringsAsFactors = FALSE
  )

  results[[grp]] <- list(pred = pred, summary = rr_row,
                          temp_ref = temp_ref, data = df)
  all_summaries[[grp]] <- rr_row

  cat("  Cold p05: exp(fit) =", round(rr_fit[1], 4), "\n")
  cat("  Heat p95: exp(fit) =", round(rr_fit[4], 4), "\n")
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

lga_table <- bind_rows(all_summaries)
write.csv(lga_table,
          file.path(tab_dir, "lga_sensitivity_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/lga_sensitivity_rr.csv\n")

# Comparison table (load SA4 results if available)
sa4_file <- file.path(tab_dir, "dlnm_results_summary.csv")
if (file.exists(sa4_file)) {
  sa4_table <- read.csv(sa4_file, stringsAsFactors = FALSE)
  sa4_compare <- sa4_table |>
    transmute(group, label, resolution = "SA4 weekly",
              rr_p95 = rr_p95, rr_p95_lo = rr_p95_lo, rr_p95_hi = rr_p95_hi)
  lga_compare <- lga_table |>
    transmute(group, label, resolution = "LGA monthly",
              rr_p95, rr_p95_lo, rr_p95_hi)
  comparison <- bind_rows(sa4_compare, lga_compare) |>
    arrange(group, resolution)
  write.csv(comparison,
            file.path(tab_dir, "sa4_vs_lga_comparison.csv"),
            row.names = FALSE)
  cat("-> Saved: outputs/tables/sa4_vs_lga_comparison.csv\n")
}


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating plots...\n")

# ER curves for each group
for (grp in names(results)) {
  res  <- results[[grp]]
  pred <- res$pred
  lbl  <- group_labels[grp]

  # For Gaussian log-rate, exponentiate predictions
  df_plot <- data.frame(
    temp  = as.numeric(names(pred$allfit)),
    rr    = exp(pred$allfit),
    rr_lo = exp(pred$alllow),
    rr_hi = exp(pred$allhigh)
  )

  p <- ggplot(df_plot, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi),
                fill = "#fc8d59", alpha = 0.2) +
    geom_line(colour = "#fc8d59", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = res$temp_ref, linetype = "dotted", colour = "grey60") +
    labs(
      title = paste0(lbl, " (LGA x monthly)"),
      subtitle = "Sensitivity analysis — different spatial/temporal resolution",
      x = "Monthly mean Tmax (°C)",
      y = "Rate ratio (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(fig_dir, paste0("lga_sensitivity_er_", grp, ".png")),
         p, width = 7, height = 5, dpi = 300)
}
cat("  -> Saved: lga_sensitivity_er_*.png\n")

# Side-by-side comparison plot (SA4 vs LGA for respiratory)
if (file.exists(sa4_file) && "resp_total" %in% names(results)) {
  cat("  Generating comparison plot (SA4 vs LGA)...\n")

  # LGA curve
  res_lga  <- results[["resp_total"]]
  pred_lga <- res_lga$pred
  df_lga <- data.frame(
    temp  = as.numeric(names(pred_lga$allfit)),
    rr    = exp(pred_lga$allfit),
    rr_lo = exp(pred_lga$alllow),
    rr_hi = exp(pred_lga$allhigh),
    source = "LGA monthly"
  )

  p_compare <- ggplot(df_lga, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "#fc8d59", alpha = 0.2) +
    geom_line(colour = "#fc8d59", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    labs(
      title = "Respiratory — LGA monthly sensitivity analysis",
      subtitle = "Compare direction/magnitude with SA4 weekly results from script 05",
      x = "Monthly mean Tmax (°C)",
      y = "Rate ratio (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(fig_dir, "lga_sensitivity_comparison.png"),
         p_compare, width = 8, height = 5, dpi = 300)
  cat("  -> Saved: lga_sensitivity_comparison.png\n")
}


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 12 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("LGA sensitivity results:\n")
lga_table |>
  select(label, resolution, rr_p95, rr_p95_lo, rr_p95_hi, n_obs, n_lga) |>
  print(row.names = FALSE)

if (file.exists(sa4_file)) {
  cat("\nSA4 vs LGA comparison (RR at p95):\n")
  comparison |>
    filter(group %in% c("cvd_total", "resp_total")) |>
    print(row.names = FALSE)
}

cat("\nOutputs:\n")
cat("  outputs/tables/lga_sensitivity_rr.csv\n")
cat("  outputs/tables/sa4_vs_lga_comparison.csv\n")
cat("  outputs/figures/lga_sensitivity_er_*.png\n")
cat("  outputs/figures/lga_sensitivity_comparison.png\n")
