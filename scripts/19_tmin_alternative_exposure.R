# ==============================================================================
# Script 19: Alternative Exposure — Minimum Temperature (Tmin)
# ==============================================================================
# Re-fits the primary DLNM models using weekly mean Tmin instead of Tmax to
# test whether nighttime temperature exposure produces different associations
# with medication dispensing.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/tmin_results_summary.csv
#   outputs/tables/tmax_vs_tmin_comparison.csv
#   outputs/figures/tmin_exposure_response_*.png
#   outputs/figures/tmax_vs_tmin_forest.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate",
#                       "ggplot2", "patchwork"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

outcome_groups <- c("cvd_total", "mh_total", "resp_total",
                    "resp_relievers", "resp_preventers")
group_labels <- c(
  cvd_total       = "Cardiovascular",
  mh_total        = "Mental Health",
  resp_total      = "Respiratory",
  resp_relievers  = "Resp. Relievers",
  resp_preventers = "Resp. Preventers"
)

cat("=" |> strrep(70), "\n")
cat("Script 19: Alternative Exposure — Tmin\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

panel <- panel |>
  mutate(
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    week_of_year = isoweek(week_start),
    year         = year(week_start),
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )

# Check Tmin coverage
n_tmin <- sum(!is.na(panel$tmin_mean))
cat("  Tmin coverage:", format(n_tmin, big.mark = ","), "/",
    format(nrow(panel), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_tmin / nrow(panel)))


# ==============================================================================
# 2. FIT DLNM WITH TMIN
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models with Tmin exposure\n")
cat("=" |> strrep(70), "\n")

results <- list()

for (grp in outcome_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmin_mean)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

  if (nrow(df) < 500) {
    cat("  Too few observations, skipping.\n")
    next
  }

  temp_range  <- range(df$tmin_mean, na.rm = TRUE)
  temp_median <- median(df$tmin_mean, na.rm = TRUE)
  cat("  Tmin range:", round(temp_range[1], 1), "to",
      round(temp_range[2], 1), "°C (median:", round(temp_median, 1), ")\n")

  # Cross-basis on Tmin
  temp_knots <- quantile(df$tmin_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmin_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))

  model <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  temp_ref <- temp_median

  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 80),
    cen = temp_ref
  )

  temp_p05 <- quantile(df$tmin_mean, 0.05, na.rm = TRUE)
  temp_p95 <- quantile(df$tmin_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmin_mean, 0.99, na.rm = TRUE)

  pred_pctl <- crosspred(cb, model,
                          at = c(temp_p05, temp_p95, temp_p99),
                          cen = temp_ref)

  rr_row <- data.frame(
    group     = grp,
    label     = group_labels[grp],
    exposure  = "Tmin",
    temp_ref  = round(temp_ref, 1),
    temp_p05  = round(temp_p05, 1),
    rr_p05    = round(pred_pctl$allRRfit[1], 4),
    rr_p05_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p05_hi = round(pred_pctl$allRRhigh[1], 4),
    temp_p95  = round(temp_p95, 1),
    rr_p95    = round(pred_pctl$allRRfit[2], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[2], 4),
    temp_p99  = round(temp_p99, 1),
    rr_p99    = round(pred_pctl$allRRfit[3], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[3], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[3], 4),
    n_obs     = nrow(df),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  results[[grp]] <- list(pred = pred, summary = rr_row,
                          temp_ref = temp_ref, data = df)

  cat("  RR at p95 (warm nights):", rr_row$rr_p95,
      "(", rr_row$rr_p95_lo, "-", rr_row$rr_p95_hi, ")\n")
  cat("  RR at p05 (cold nights):", rr_row$rr_p05,
      "(", rr_row$rr_p05_lo, "-", rr_row$rr_p05_hi, ")\n")
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

tmin_table <- bind_rows(lapply(results, function(r) r$summary))
write.csv(tmin_table,
          file.path(tab_dir, "tmin_results_summary.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/tmin_results_summary.csv\n")

# Comparison table with Tmax results
tmax_file <- file.path(tab_dir, "dlnm_results_summary.csv")
if (file.exists(tmax_file)) {
  tmax_table <- read.csv(tmax_file, stringsAsFactors = FALSE)
  tmax_compare <- tmax_table |>
    filter(group %in% outcome_groups) |>
    transmute(group, label, exposure = "Tmax",
              rr_p95, rr_p95_lo = rr_p95_lo, rr_p95_hi = rr_p95_hi)
  tmin_compare <- tmin_table |>
    transmute(group, label, exposure = "Tmin",
              rr_p95, rr_p95_lo, rr_p95_hi)
  comparison <- bind_rows(tmax_compare, tmin_compare) |>
    arrange(group, exposure)
  write.csv(comparison,
            file.path(tab_dir, "tmax_vs_tmin_comparison.csv"),
            row.names = FALSE)
  cat("-> Saved: outputs/tables/tmax_vs_tmin_comparison.csv\n")
}


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating Tmin plots...\n")

# ER curves
for (grp in names(results)) {
  res  <- results[[grp]]
  pred <- res$pred
  lbl  <- group_labels[grp]

  df_plot <- data.frame(
    temp  = as.numeric(names(pred$allRRfit)),
    rr    = pred$allRRfit,
    rr_lo = pred$allRRlow,
    rr_hi = pred$allRRhigh
  )

  p <- ggplot(df_plot, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "#fc8d59", alpha = 0.2) +
    geom_line(colour = "#fc8d59", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = res$temp_ref, linetype = "dotted", colour = "grey60") +
    labs(
      title = paste0(lbl, " — Tmin exposure-response"),
      x = "Weekly mean Tmin (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(fig_dir, paste0("tmin_exposure_response_", grp, ".png")),
         p, width = 7, height = 5, dpi = 300)
}
cat("  -> Saved: tmin_exposure_response_*.png\n")

# Comparison forest plot (Tmax vs Tmin)
if (exists("comparison")) {
  p_compare <- ggplot(comparison,
                       aes(x = rr_p95, y = reorder(label, rr_p95),
                           colour = exposure)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.6,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    scale_colour_manual(
      values = c("Tmax" = "#d73027", "Tmin" = "#fc8d59"),
      name = "Exposure"
    ) +
    labs(
      title = "Cumulative RR at 95th percentile: Tmax vs Tmin",
      x = "Cumulative RR (95% CI)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom")

  ggsave(file.path(fig_dir, "tmax_vs_tmin_forest.png"),
         p_compare, width = 9, height = 5, dpi = 300)
  cat("  -> Saved: tmax_vs_tmin_forest.png\n")
}


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 19 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Tmin DLNM results:\n")
print(as.data.frame(
  tmin_table |> select(label, rr_p05, rr_p05_lo, rr_p05_hi,
                        rr_p95, rr_p95_lo, rr_p95_hi)
), row.names = FALSE)

if (exists("comparison")) {
  cat("\nTmax vs Tmin comparison (RR at p95):\n")
  print(as.data.frame(comparison), row.names = FALSE)
}

cat("\nOutputs:\n")
cat("  outputs/tables/tmin_results_summary.csv\n")
cat("  outputs/tables/tmax_vs_tmin_comparison.csv\n")
cat("  outputs/figures/tmin_exposure_response_*.png\n")
cat("  outputs/figures/tmax_vs_tmin_forest.png\n")
