# ==============================================================================
# Script 30: Three Fourier Pairs Sensitivity
# ==============================================================================
# Re-fits primary DLNM models using 3 Fourier pairs for seasonality instead
# of 2, as suggested by the QAIC comparison (script 27).
#
# Compares:
#   - Primary model (2 Fourier pairs: sin1/cos1, sin2/cos2)
#   - Updated model (3 Fourier pairs: + sin3/cos3)
#
# If results are substantively identical, the 2-pair model is more
# parsimonious and standard. If they differ meaningfully, the 3-pair
# model should be considered as the primary specification.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/fourier_sensitivity.csv
#   outputs/figures/fourier_sensitivity_forest.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate", "ggplot2"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)

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
cat("Script 30: Fourier Pairs Sensitivity (2 vs 3)\n")
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
    cos2 = cos(4 * pi * week_of_year / 52),
    sin3 = sin(6 * pi * week_of_year / 52),
    cos3 = cos(6 * pi * week_of_year / 52)
  )


# ==============================================================================
# 2. FIT MODELS AND COMPARE
# ==============================================================================

all_results <- list()

for (grp in outcome_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(precip_total)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))
  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p05 <- quantile(df$tmax_mean, 0.05, na.rm = TRUE)

  # --- Model A: 2 Fourier pairs (primary) ---
  cat("  Fitting 2-pair model...\n")
  model_2f <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  pred_2f <- crosspred(cb, model_2f,
                        at = c(temp_p05, temp_p95),
                        cen = temp_ref)

  disp_2f <- summary(model_2f)$dispersion
  dev_2f  <- deviance(model_2f)
  k_2f    <- length(coef(model_2f))

  # --- Model B: 3 Fourier pairs ---
  cat("  Fitting 3-pair model...\n")
  model_3f <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 + sin3 + cos3 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  pred_3f <- crosspred(cb, model_3f,
                        at = c(temp_p05, temp_p95),
                        cen = temp_ref)

  disp_3f <- summary(model_3f)$dispersion
  dev_3f  <- deviance(model_3f)
  k_3f    <- length(coef(model_3f))

  # QAIC comparison (using 2-pair dispersion as reference)
  qaic_2f <- dev_2f / disp_2f + 2 * k_2f
  qaic_3f <- dev_3f / disp_2f + 2 * k_3f

  cat("  2-pair: RR@p95 =", round(pred_2f$allRRfit[2], 4),
      " disp =", round(disp_2f, 2), " QAIC =", round(qaic_2f, 1), "\n")
  cat("  3-pair: RR@p95 =", round(pred_3f$allRRfit[2], 4),
      " disp =", round(disp_3f, 2), " QAIC =", round(qaic_3f, 1), "\n")
  cat("  Delta QAIC:", round(qaic_2f - qaic_3f, 1),
      "(positive = 3-pair preferred)\n")

  # Percent change in RR
  pct_change_p95 <- 100 * (pred_3f$allRRfit[2] - pred_2f$allRRfit[2]) / pred_2f$allRRfit[2]
  cat("  RR@p95 % change:", round(pct_change_p95, 2), "%\n")

  all_results[[paste0(grp, "_2f")]] <- data.frame(
    group       = grp,
    label       = group_labels[grp],
    fourier     = "2 pairs",
    rr_p05      = round(pred_2f$allRRfit[1], 4),
    rr_p05_lo   = round(pred_2f$allRRlow[1], 4),
    rr_p05_hi   = round(pred_2f$allRRhigh[1], 4),
    rr_p95      = round(pred_2f$allRRfit[2], 4),
    rr_p95_lo   = round(pred_2f$allRRlow[2], 4),
    rr_p95_hi   = round(pred_2f$allRRhigh[2], 4),
    dispersion  = round(disp_2f, 2),
    qaic        = round(qaic_2f, 1),
    n_params    = k_2f,
    stringsAsFactors = FALSE
  )

  all_results[[paste0(grp, "_3f")]] <- data.frame(
    group       = grp,
    label       = group_labels[grp],
    fourier     = "3 pairs",
    rr_p05      = round(pred_3f$allRRfit[1], 4),
    rr_p05_lo   = round(pred_3f$allRRlow[1], 4),
    rr_p05_hi   = round(pred_3f$allRRhigh[1], 4),
    rr_p95      = round(pred_3f$allRRfit[2], 4),
    rr_p95_lo   = round(pred_3f$allRRlow[2], 4),
    rr_p95_hi   = round(pred_3f$allRRhigh[2], 4),
    dispersion  = round(disp_3f, 2),
    qaic        = round(qaic_3f, 1),
    n_params    = k_3f,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

results_table <- bind_rows(all_results)
write.csv(results_table,
          file.path(tab_dir, "fourier_sensitivity.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/fourier_sensitivity.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison forest plot...\n")

p <- ggplot(results_table,
            aes(x = rr_p95, y = label, colour = fourier)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.3, linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("2 pairs" = "steelblue",
                                  "3 pairs" = "#d73027"),
                      name = "Fourier harmonics") +
  labs(
    title = "Cumulative RR at 95th percentile: 2 vs 3 Fourier pairs",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "fourier_sensitivity_forest.png"),
       p, width = 10, height = 5, dpi = 300)
cat("  -> Saved: fourier_sensitivity_forest.png\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 30 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Fourier sensitivity comparison:\n")
print(as.data.frame(
  results_table |>
    select(label, fourier, rr_p95, rr_p95_lo, rr_p95_hi, dispersion, qaic)
), row.names = FALSE)

# Summary assessment
cat("\nAssessment:\n")
for (grp in outcome_groups) {
  r2 <- results_table |> filter(group == grp, fourier == "2 pairs")
  r3 <- results_table |> filter(group == grp, fourier == "3 pairs")
  pct <- abs(100 * (r3$rr_p95 - r2$rr_p95) / r2$rr_p95)
  dq  <- r2$qaic - r3$qaic
  cat("  ", group_labels[grp], ": RR change =", round(pct, 2),
      "%, delta QAIC =", round(dq, 1))
  if (pct < 1 & dq > 0) {
    cat(" -> 3 pairs preferred by QAIC but RR negligibly different\n")
  } else if (pct >= 1) {
    cat(" -> substantive difference, consider 3 pairs\n")
  } else {
    cat(" -> 2 pairs adequate\n")
  }
}

cat("\nOutputs:\n")
cat("  outputs/tables/fourier_sensitivity.csv\n")
cat("  outputs/figures/fourier_sensitivity_forest.png\n")
