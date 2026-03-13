# ==============================================================================
# Script 24: Minimum Morbidity Temperature (MMT) Analysis
# ==============================================================================
# Identifies the temperature at which medication dispensing is minimised (MMT)
# rather than assuming the median is the optimal reference. This is standard
# practice in climate-health epidemiology (Gasparrini et al. 2015).
#
# For mortality studies the MMT is the "minimum mortality temperature"; here
# we adapt it to find the temperature that minimises prescribing.
#
# Approach:
#   1. Fit DLNM with temperature centred at median (as usual)
#   2. Extract the cumulative exposure-response curve
#   3. Find the temperature that gives the minimum predicted RR
#   4. Re-fit centred at the MMT and report RR at cold/heat percentiles
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/mmt_results.csv
#   outputs/figures/mmt_curves.png
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
cat("Script 24: Minimum Morbidity Temperature (MMT)\n")
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


# ==============================================================================
# 2. FIND MMT FOR EACH MEDICATION GROUP
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Finding minimum morbidity temperature (MMT)\n")
cat("=" |> strrep(70), "\n")

results <- list()

for (grp in outcome_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 500) { cat("  Too few obs, skipping.\n"); next }

  # Fit model
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
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

  temp_median <- median(df$tmax_mean, na.rm = TRUE)
  temp_range  <- range(df$tmax_mean, na.rm = TRUE)

  # Fine-grained prediction to find minimum
  temp_seq <- seq(temp_range[1], temp_range[2], length.out = 200)
  pred_fine <- crosspred(cb, model, at = temp_seq, cen = temp_median)

  # The allRRfit gives RR relative to median.
  # To find MMT, we need the temperature that minimises the overall
  # cumulative association (i.e., minimum of the cumulative fit)
  # This is equivalent to finding the minimum of allRRfit
  min_idx <- which.min(pred_fine$allRRfit)
  mmt <- temp_seq[min_idx]
  mmt_percentile <- ecdf(df$tmax_mean)(mmt) * 100

  cat("  MMT:", round(mmt, 1), "°C (", round(mmt_percentile, 1),
      "th percentile)\n")
  cat("  Median:", round(temp_median, 1), "°C\n")
  cat("  Difference:", round(mmt - temp_median, 1), "°C\n")

  # Re-fit centred at MMT
  pred_mmt <- crosspred(cb, model, at = temp_seq, cen = mmt)

  # RR at key percentiles relative to MMT
  temp_p05 <- quantile(df$tmax_mean, 0.05, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  pred_pctl_mmt <- crosspred(cb, model,
                              at = c(temp_p05, temp_p95, temp_p99),
                              cen = mmt)

  rr_row <- data.frame(
    group          = grp,
    label          = group_labels[grp],
    temp_median    = round(temp_median, 1),
    mmt            = round(mmt, 1),
    mmt_percentile = round(mmt_percentile, 1),
    mmt_minus_median = round(mmt - temp_median, 1),
    # RR at p05 relative to MMT
    rr_p05_mmt     = round(pred_pctl_mmt$allRRfit[1], 4),
    rr_p05_mmt_lo  = round(pred_pctl_mmt$allRRlow[1], 4),
    rr_p05_mmt_hi  = round(pred_pctl_mmt$allRRhigh[1], 4),
    # RR at p95 relative to MMT
    rr_p95_mmt     = round(pred_pctl_mmt$allRRfit[2], 4),
    rr_p95_mmt_lo  = round(pred_pctl_mmt$allRRlow[2], 4),
    rr_p95_mmt_hi  = round(pred_pctl_mmt$allRRhigh[2], 4),
    # RR at p99 relative to MMT
    rr_p99_mmt     = round(pred_pctl_mmt$allRRfit[3], 4),
    rr_p99_mmt_lo  = round(pred_pctl_mmt$allRRlow[3], 4),
    rr_p99_mmt_hi  = round(pred_pctl_mmt$allRRhigh[3], 4),
    n_obs = nrow(df),
    stringsAsFactors = FALSE
  )

  results[[grp]] <- list(
    pred_median = pred_fine,
    pred_mmt    = pred_mmt,
    summary     = rr_row,
    mmt         = mmt,
    temp_median = temp_median,
    temp_seq    = temp_seq
  )

  cat("  RR at p95 (vs MMT):", rr_row$rr_p95_mmt,
      "(", rr_row$rr_p95_mmt_lo, "-", rr_row$rr_p95_mmt_hi, ")\n")
  cat("  RR at p05 (vs MMT):", rr_row$rr_p05_mmt,
      "(", rr_row$rr_p05_mmt_lo, "-", rr_row$rr_p05_mmt_hi, ")\n")
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

mmt_table <- bind_rows(lapply(results, function(r) r$summary))
write.csv(mmt_table,
          file.path(tab_dir, "mmt_results.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/mmt_results.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating MMT plots...\n")

plot_mmt <- function(grp) {
  res <- results[[grp]]
  pred <- res$pred_mmt
  lbl <- group_labels[grp]

  df_plot <- data.frame(
    temp  = as.numeric(names(pred$allRRfit)),
    rr    = pred$allRRfit,
    rr_lo = pred$allRRlow,
    rr_hi = pred$allRRhigh
  )

  p <- ggplot(df_plot, aes(x = temp, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "steelblue", alpha = 0.2) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = res$mmt, linetype = "solid",
               colour = "#d73027", linewidth = 0.7) +
    geom_vline(xintercept = res$temp_median, linetype = "dotted",
               colour = "grey50") +
    annotate("text", x = res$mmt + 0.5, y = Inf,
             label = paste0("MMT = ", round(res$mmt, 1), "°C"),
             vjust = 1.5, hjust = 0, colour = "#d73027", size = 3.5) +
    annotate("text", x = res$temp_median - 0.5, y = Inf,
             label = paste0("Median = ", round(res$temp_median, 1), "°C"),
             vjust = 3, hjust = 1, colour = "grey50", size = 3) +
    labs(
      title = lbl,
      subtitle = paste0("MMT = ", round(res$mmt, 1),
                         "°C (", round(res$summary$mmt_percentile, 0),
                         "th percentile)"),
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs MMT)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  return(p)
}

# Individual plots
for (grp in names(results)) {
  p <- plot_mmt(grp)
  ggsave(file.path(fig_dir, paste0("mmt_curve_", grp, ".png")),
         p, width = 7, height = 5, dpi = 300)
}
cat("  -> Saved: mmt_curve_*.png\n")

# Panel of primary groups
primary <- intersect(c("cvd_total", "mh_total", "resp_total"), names(results))
plots <- lapply(primary, plot_mmt)
combined <- wrap_plots(plots, ncol = 3)
ggsave(file.path(fig_dir, "mmt_curves_panel.png"),
       combined, width = 15, height = 5, dpi = 300)
cat("  -> Saved: mmt_curves_panel.png\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 24 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Minimum Morbidity Temperature (MMT) results:\n")
print(as.data.frame(
  mmt_table |>
    select(label, temp_median, mmt, mmt_percentile, mmt_minus_median,
           rr_p05_mmt, rr_p95_mmt)
), row.names = FALSE)

cat("\nInterpretation:\n")
cat("  MMT = temperature at which prescribing is minimised.\n")
cat("  If MMT >> median: heat is beneficial up to a point, then harmful.\n")
cat("  If MMT ~ median: median is an appropriate reference.\n")
cat("  If MMT << median: cold drives more prescribing than heat.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/mmt_results.csv\n")
cat("  outputs/figures/mmt_curve_*.png\n")
cat("  outputs/figures/mmt_curves_panel.png\n")
