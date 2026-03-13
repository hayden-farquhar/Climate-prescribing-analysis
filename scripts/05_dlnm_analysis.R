# ==============================================================================
# Script 05: Distributed Lag Non-Linear Models (DLNM) Analysis
# ==============================================================================
# Fits DLNMs examining the association between weekly heat exposure and
# medication dispensing across Australian SA4 regions.
#
# Models:
#   - Cross-basis: natural cubic spline for temperature x ns for lag (0-12 weeks)
#   - Outcome: weekly prescription count per SA4
#   - Offset: log(population) for rate interpretation
#   - Covariates: long-term trend, seasonality (Fourier), rainfall
#   - Fitted separately for each medication group:
#       cvd_total, mh_total, resp_total (plus subgroups)
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv   (from script 04)
#   data/reference/ERP_LGA_2001-2024.xlsx (for population denominators)
#
# Outputs:
#   outputs/tables/dlnm_results_summary.csv
#   outputs/figures/exposure_response_*.png
#   outputs/figures/lag_response_*.png
#   outputs/figures/overall_cumulative_*.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "mgcv", "ggplot2", "dplyr",
#                       "lubridate", "patchwork", "scales"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(scales)

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("=" |> strrep(70), "\n")
cat("Script 05: DLNM Analysis — Heat and Medication Dispensing\n")
cat("=" |> strrep(70), "\n\n")

# Load panel
cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

cat("  Rows:", format(nrow(panel), big.mark = ","), "\n")
cat("  Date range:", as.character(min(panel$week_start)),
    "to", as.character(max(panel$week_start)), "\n")
cat("  SA4 regions:", n_distinct(panel$sa4_code), "\n")
cat("  Analysis groups:", paste(sort(unique(panel$analysis_group)), collapse = ", "), "\n")

# Check climate coverage
n_total   <- nrow(panel)
n_climate <- sum(!is.na(panel$tmax_mean))
cat("  Climate coverage:", format(n_climate, big.mark = ","), "/",
    format(n_total, big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_climate / n_total))

if (n_climate / n_total < 0.5) {
  cat("\nWARNING: Low climate coverage. Re-run script 04 after ERA5 download.\n")
}

# --- Time variables for trend and seasonality ---
panel <- panel |>
  mutate(
    year       = year(week_start),
    week_of_year = isoweek(week_start),
    time_index = as.numeric(week_start - min(week_start)) / 7,  # weeks since start
    # Fourier terms for seasonality (annual cycle, 52-week period)
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )

# ==============================================================================
# 2. DLNM SPECIFICATION
# ==============================================================================

# Lag range: 0-12 weeks (3 months equivalent)
max_lag <- 12

# Medication groups to model (totals first, then subgroups)
outcome_groups <- c(
  "cvd_total",          # Cardiovascular (primary)
  "mh_total",           # Mental health (primary)
  "resp_total",         # Respiratory (primary)
  "mh_antidepressants", # Antidepressant subgroup
  "mh_anxiolytics",     # Anxiolytic subgroup
  "resp_relievers",     # Respiratory reliever subgroup
  "resp_preventers"     # Respiratory preventer subgroup
)

# Labels for plots
group_labels <- c(
  cvd_total          = "Cardiovascular",
  mh_total           = "Mental Health",
  resp_total         = "Respiratory",
  mh_antidepressants = "Antidepressants",
  mh_anxiolytics     = "Anxiolytics",
  resp_relievers     = "Resp. Relievers",
  resp_preventers    = "Resp. Preventers"
)


# ==============================================================================
# 3. FIT DLNM FOR EACH MEDICATION GROUP
# ==============================================================================

fit_dlnm <- function(data, group_name) {
  #' Fit a DLNM for a single medication group.
  #'
  #' Model: quasi-Poisson GLM with cross-basis for temperature,
  #' natural spline for long-term trend, and Fourier terms for seasonality.
  #'
  #' @param data The full panel dataframe

  #' @param group_name Character: analysis_group value to filter
  #' @return List with model, crossbasis, predictions, and summary stats

  cat("\n--- Fitting DLNM for:", group_labels[group_name], "---\n")

  # Filter to this group and remove missing data
  df <- data |>
    filter(analysis_group == group_name,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")
  cat("  SA4 regions:", n_distinct(df$sa4_code), "\n")
  cat("  Date range:", as.character(min(df$week_start)),
      "to", as.character(max(df$week_start)), "\n")

  if (nrow(df) < 500) {
    cat("  WARNING: Too few observations, skipping.\n")
    return(NULL)
  }

  # Temperature distribution for this dataset
  temp_range <- range(df$tmax_mean, na.rm = TRUE)
  temp_median <- median(df$tmax_mean, na.rm = TRUE)
  cat("  Temperature range:", round(temp_range[1], 1), "to",
      round(temp_range[2], 1), "°C (median:", round(temp_median, 1), ")\n")

  # --- Cross-basis specification ---
  # Temperature dimension: natural cubic spline with knots at 10th, 50th, 90th percentiles
  temp_knots <- quantile(df$tmax_mean,
                         probs = c(0.10, 0.50, 0.90),
                         na.rm = TRUE)

  # Lag dimension: natural cubic spline, 3 equally-spaced knots in log scale
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  cat("  Cross-basis: ns(temp, knots=", paste(round(temp_knots, 1), collapse = ","),
      ") x ns(lag, 3 log-knots)\n")

  # --- Fit quasi-Poisson GLM ---
  # Offset not available without population at SA4 x week level,
  # so use crude_rate as outcome (prescriptions per 100k population)
  # or count with SA4 fixed effects
  cat("  Fitting model...\n")

  model <- glm(
    count ~ cb +
      ns(time_index, df = round(n_distinct(df$year) * 2)) +  # ~2 df per year for trend
      sin1 + cos1 + sin2 + cos2 +                             # Seasonality
      precip_total +                                           # Rainfall control
      factor(sa4_code),                                        # SA4 fixed effects
    data = df,
    family = quasipoisson(link = "log")
  )

  cat("  Dispersion:", round(summary(model)$dispersion, 2), "\n")

  # --- Predictions ---
  # Reference temperature: median (roughly "comfortable" temperature)
  temp_ref <- temp_median

  # Predict over full temperature range
  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 100),
    cen = temp_ref
  )

  # Key percentile predictions
  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  pred_percentiles <- crosspred(
    cb, model,
    at = c(temp_p90, temp_p95, temp_p99),
    cen = temp_ref
  )

  # Overall cumulative RR at key percentiles
  rr_summary <- data.frame(
    group     = group_name,
    label     = group_labels[group_name],
    temp_ref  = round(temp_ref, 1),
    temp_p90  = round(temp_p90, 1),
    temp_p95  = round(temp_p95, 1),
    temp_p99  = round(temp_p99, 1),
    rr_p90    = round(pred_percentiles$allRRfit["temp_p90"], 4),
    rr_p90_lo = round(pred_percentiles$allRRlow["temp_p90"], 4),
    rr_p90_hi = round(pred_percentiles$allRRhigh["temp_p90"], 4),
    rr_p95    = round(pred_percentiles$allRRfit["temp_p95"], 4),
    rr_p95_lo = round(pred_percentiles$allRRlow["temp_p95"], 4),
    rr_p95_hi = round(pred_percentiles$allRRhigh["temp_p95"], 4),
    rr_p99    = round(pred_percentiles$allRRfit["temp_p99"], 4),
    rr_p99_lo = round(pred_percentiles$allRRlow["temp_p99"], 4),
    rr_p99_hi = round(pred_percentiles$allRRhigh["temp_p99"], 4),
    n_obs     = nrow(df),
    n_sa4     = n_distinct(df$sa4_code),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  # Fix row names for percentile lookup
  pctl_names <- as.character(round(c(temp_p90, temp_p95, temp_p99), 1))
  if (length(pred_percentiles$allRRfit) >= 3) {
    rr_summary$rr_p90    <- round(pred_percentiles$allRRfit[1], 4)
    rr_summary$rr_p90_lo <- round(pred_percentiles$allRRlow[1], 4)
    rr_summary$rr_p90_hi <- round(pred_percentiles$allRRhigh[1], 4)
    rr_summary$rr_p95    <- round(pred_percentiles$allRRfit[2], 4)
    rr_summary$rr_p95_lo <- round(pred_percentiles$allRRlow[2], 4)
    rr_summary$rr_p95_hi <- round(pred_percentiles$allRRhigh[2], 4)
    rr_summary$rr_p99    <- round(pred_percentiles$allRRfit[3], 4)
    rr_summary$rr_p99_lo <- round(pred_percentiles$allRRlow[3], 4)
    rr_summary$rr_p99_hi <- round(pred_percentiles$allRRhigh[3], 4)
  }

  cat("  Cumulative RR at 95th pctl:", rr_summary$rr_p95,
      "(", rr_summary$rr_p95_lo, "-", rr_summary$rr_p95_hi, ")\n")

  return(list(
    model      = model,
    cb         = cb,
    pred       = pred,
    pred_pctl  = pred_percentiles,
    summary    = rr_summary,
    data       = df,
    temp_ref   = temp_ref,
    group_name = group_name
  ))
}


# ==============================================================================
# 4. RUN ALL MODELS
# ==============================================================================

cat("\n" , "=" |> strrep(70), "\n")
cat("Fitting DLNM models for", length(outcome_groups), "medication groups\n")
cat("=" |> strrep(70), "\n")

results <- list()
for (grp in outcome_groups) {
  res <- fit_dlnm(panel, grp)
  if (!is.null(res)) {
    results[[grp]] <- res
  }
}

cat("\n\nModels fitted:", length(results), "/", length(outcome_groups), "\n")


# ==============================================================================
# 5. SUMMARY TABLE
# ==============================================================================

if (length(results) > 0) {
  summary_table <- bind_rows(lapply(results, function(r) r$summary))

  write.csv(summary_table,
            file.path(tab_dir, "dlnm_results_summary.csv"),
            row.names = FALSE)
  cat("\n-> Saved: outputs/tables/dlnm_results_summary.csv\n")

  cat("\nCumulative RR Summary (vs median temperature):\n")
  print(
    summary_table |>
      select(label, temp_ref, temp_p95, rr_p95, rr_p95_lo, rr_p95_hi,
             rr_p99, rr_p99_lo, rr_p99_hi, n_obs, dispersion),
    row.names = FALSE
  )
}


# ==============================================================================
# 6. VISUALISATION
# ==============================================================================

# --- 6a: Overall cumulative exposure-response curves ---
plot_exposure_response <- function(res) {
  pred <- res$pred
  grp  <- res$group_name
  lbl  <- group_labels[grp]

  df_plot <- data.frame(
    temp   = as.numeric(names(pred$allRRfit)),
    rr     = pred$allRRfit,
    rr_lo  = pred$allRRlow,
    rr_hi  = pred$allRRhigh
  )

  p <- ggplot(df_plot, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "steelblue", alpha = 0.2) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = res$temp_ref, linetype = "dotted", colour = "grey60") +
    labs(
      title = lbl,
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  return(p)
}

if (length(results) >= 3) {
  cat("\nGenerating exposure-response plots...\n")

  # Panel of primary groups
  primary_groups <- intersect(c("cvd_total", "mh_total", "resp_total"), names(results))

  plots <- lapply(primary_groups, function(g) plot_exposure_response(results[[g]]))
  combined <- wrap_plots(plots, ncol = 3)

  ggsave(file.path(fig_dir, "exposure_response_primary.png"),
         combined, width = 14, height = 5, dpi = 300)
  cat("  -> Saved: outputs/figures/exposure_response_primary.png\n")

  # Individual plots for all groups
  for (grp in names(results)) {
    p <- plot_exposure_response(results[[grp]])
    ggsave(file.path(fig_dir, paste0("exposure_response_", grp, ".png")),
           p, width = 6, height = 5, dpi = 300)
  }
  cat("  -> Saved individual exposure-response plots\n")
}

# --- 6b: Lag-response curves at 95th percentile ---
plot_lag_response <- function(res) {
  pred  <- res$pred
  grp   <- res$group_name
  lbl   <- group_labels[grp]
  temp_p95 <- res$summary$temp_p95

  # Find the closest temperature in pred$predvar to p95
  idx <- which.min(abs(as.numeric(names(pred$allRRfit)) - temp_p95))
  temp_at <- as.numeric(names(pred$allRRfit))[idx]

  # Extract lag-specific RRs at this temperature
  lag_rr    <- pred$matRRfit[idx, ]
  lag_rr_lo <- pred$matRRlow[idx, ]
  lag_rr_hi <- pred$matRRhigh[idx, ]

  df_lag <- data.frame(
    lag    = 0:max_lag,
    rr     = lag_rr,
    rr_lo  = lag_rr_lo,
    rr_hi  = lag_rr_hi
  )

  p <- ggplot(df_lag, aes(x = lag, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "coral", alpha = 0.2) +
    geom_line(colour = "coral", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(breaks = seq(0, max_lag, by = 2)) +
    labs(
      title = paste0(lbl, " — Lag response at ", round(temp_at, 1), "°C (p95)"),
      x = "Lag (weeks)",
      y = "RR"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  return(p)
}

if (length(results) >= 3) {
  cat("Generating lag-response plots...\n")

  primary_groups <- intersect(c("cvd_total", "mh_total", "resp_total"), names(results))
  plots <- lapply(primary_groups, function(g) plot_lag_response(results[[g]]))
  combined <- wrap_plots(plots, ncol = 3)

  ggsave(file.path(fig_dir, "lag_response_primary.png"),
         combined, width = 14, height = 5, dpi = 300)
  cat("  -> Saved: outputs/figures/lag_response_primary.png\n")

  for (grp in names(results)) {
    p <- plot_lag_response(results[[grp]])
    ggsave(file.path(fig_dir, paste0("lag_response_", grp, ".png")),
           p, width = 6, height = 5, dpi = 300)
  }
  cat("  -> Saved individual lag-response plots\n")
}


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 05 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Models fitted:\n")
for (grp in names(results)) {
  r <- results[[grp]]$summary
  cat(sprintf("  %-20s  RR@p95 = %.3f (%.3f-%.3f)  n=%s\n",
              r$label, r$rr_p95, r$rr_p95_lo, r$rr_p95_hi,
              format(r$n_obs, big.mark = ",")))
}

cat("\nOutputs:\n")
cat("  outputs/tables/dlnm_results_summary.csv\n")
cat("  outputs/figures/exposure_response_*.png\n")
cat("  outputs/figures/lag_response_*.png\n")

cat("\nNext steps:\n")
cat("  - Script 06: SEIFA equity stratification (interaction models)\n")
cat("  - Script 07: Black Summer bushfire DiD analysis\n")
