# ==============================================================================
# Script 25: Influence Diagnostics
# ==============================================================================
# Identifies influential observations (SA4-weeks) that may disproportionately
# drive DLNM results. Uses Cook's distance and hat values to flag outliers.
#
# Approach:
#   1. Fit primary DLNM models
#   2. Compute Cook's distance and hat (leverage) values
#   3. Flag observations exceeding standard thresholds
#   4. Re-fit models excluding top influential observations
#   5. Compare results with and without influential points
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/influence_diagnostics.csv
#   outputs/tables/influence_sensitivity.csv
#   outputs/figures/influence_cooks_*.png
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

primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

cat("=" |> strrep(70), "\n")
cat("Script 25: Influence Diagnostics\n")
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
# 2. COMPUTE INFLUENCE MEASURES
# ==============================================================================

diag_results <- list()
sensitivity_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

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

  # Cook's distance
  cooks <- cooks.distance(model)
  hat_vals <- hatvalues(model)

  # Standard thresholds
  n <- nrow(df)
  p <- length(coef(model))
  cooks_threshold <- 4 / n
  hat_threshold   <- 2 * p / n

  n_high_cooks <- sum(cooks > cooks_threshold, na.rm = TRUE)
  n_high_hat   <- sum(hat_vals > hat_threshold, na.rm = TRUE)

  cat("  Cook's D threshold (4/n):", format(cooks_threshold, digits = 4), "\n")
  cat("  Obs exceeding Cook's D threshold:", n_high_cooks,
      "(", round(100 * n_high_cooks / n, 2), "%)\n")
  cat("  Hat threshold (2p/n):", round(hat_threshold, 4), "\n")
  cat("  Obs exceeding hat threshold:", n_high_hat,
      "(", round(100 * n_high_hat / n, 2), "%)\n")
  cat("  Max Cook's D:", format(max(cooks, na.rm = TRUE), digits = 4), "\n")

  # Identify top 10 most influential observations
  top_idx <- order(cooks, decreasing = TRUE)[1:min(10, n)]
  top_influential <- df[top_idx, c("sa4_code", "week_start", "count", "tmax_mean")]
  top_influential$cooks_d <- cooks[top_idx]
  top_influential$hat     <- hat_vals[top_idx]
  cat("  Top 10 influential observations:\n")
  print(as.data.frame(top_influential), row.names = FALSE)

  diag_results[[grp]] <- data.frame(
    group              = grp,
    label              = group_labels[grp],
    n_obs              = n,
    n_params           = p,
    cooks_threshold    = round(cooks_threshold, 6),
    n_high_cooks       = n_high_cooks,
    pct_high_cooks     = round(100 * n_high_cooks / n, 2),
    max_cooks          = round(max(cooks, na.rm = TRUE), 6),
    mean_cooks         = round(mean(cooks, na.rm = TRUE), 6),
    hat_threshold      = round(hat_threshold, 4),
    n_high_hat         = n_high_hat,
    pct_high_hat       = round(100 * n_high_hat / n, 2),
    stringsAsFactors   = FALSE
  )

  # --- Cook's D plot ---
  cooks_df <- data.frame(
    idx = seq_along(cooks),
    cooks_d = cooks
  )

  p_cooks <- ggplot(cooks_df, aes(x = idx, y = cooks_d)) +
    geom_point(alpha = 0.15, size = 0.5, colour = "steelblue") +
    geom_hline(yintercept = cooks_threshold, linetype = "dashed",
               colour = "#d73027", linewidth = 0.7) +
    annotate("text", x = n * 0.95, y = cooks_threshold * 1.5,
             label = paste0("4/n = ", format(cooks_threshold, digits = 3)),
             colour = "#d73027", hjust = 1, size = 3.5) +
    labs(
      title = paste0(group_labels[grp], " — Cook's distance"),
      subtitle = paste0(n_high_cooks, " obs (", round(100 * n_high_cooks / n, 1),
                        "%) exceed threshold"),
      x = "Observation index",
      y = "Cook's distance"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, paste0("influence_cooks_", grp, ".png")),
         p_cooks, width = 8, height = 5, dpi = 300)

  # --- Sensitivity: re-fit excluding high-influence observations ---
  cat("  Re-fitting without high-influence observations...\n")
  keep <- cooks <= cooks_threshold | is.na(cooks)
  df_trim <- df[keep, ]
  cat("  Trimmed obs:", format(nrow(df_trim), big.mark = ","),
      "(removed", n - nrow(df_trim), ")\n")

  cb_trim <- crossbasis(
    df_trim$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  model_trim <- glm(
    count ~ cb_trim +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df_trim,
    family = quasipoisson(link = "log")
  )

  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)

  # Original model predictions
  pred_orig <- crosspred(cb, model, at = temp_p95, cen = temp_ref)
  # Trimmed model predictions
  pred_trim <- crosspred(cb_trim, model_trim, at = temp_p95, cen = temp_ref)

  sensitivity_results[[grp]] <- data.frame(
    group      = grp,
    label      = group_labels[grp],
    rr_full    = round(pred_orig$allRRfit[1], 4),
    rr_full_lo = round(pred_orig$allRRlow[1], 4),
    rr_full_hi = round(pred_orig$allRRhigh[1], 4),
    rr_trim    = round(pred_trim$allRRfit[1], 4),
    rr_trim_lo = round(pred_trim$allRRlow[1], 4),
    rr_trim_hi = round(pred_trim$allRRhigh[1], 4),
    n_removed  = n - nrow(df_trim),
    pct_removed = round(100 * (n - nrow(df_trim)) / n, 2),
    stringsAsFactors = FALSE
  )

  cat("  Full model RR at p95:", pred_orig$allRRfit[1],
      "(", pred_orig$allRRlow[1], "-", pred_orig$allRRhigh[1], ")\n")
  cat("  Trimmed RR at p95:", pred_trim$allRRfit[1],
      "(", pred_trim$allRRlow[1], "-", pred_trim$allRRhigh[1], ")\n")
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

diag_table <- bind_rows(diag_results)
write.csv(diag_table,
          file.path(tab_dir, "influence_diagnostics.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/influence_diagnostics.csv\n")

sens_table <- bind_rows(sensitivity_results)
write.csv(sens_table,
          file.path(tab_dir, "influence_sensitivity.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/influence_sensitivity.csv\n")


# ==============================================================================
# 4. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 25 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Influence diagnostics:\n")
print(as.data.frame(
  diag_table |>
    select(label, n_obs, n_high_cooks, pct_high_cooks, max_cooks)
), row.names = FALSE)

cat("\nSensitivity to influential observations (RR at p95):\n")
print(as.data.frame(
  sens_table |>
    select(label, rr_full, rr_trim, n_removed, pct_removed)
), row.names = FALSE)

cat("\nInterpretation:\n")
cat("  If trimmed RR is similar to full RR, results are robust to outliers.\n")
cat("  Large changes suggest specific SA4-weeks are driving the association.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/influence_diagnostics.csv\n")
cat("  outputs/tables/influence_sensitivity.csv\n")
cat("  outputs/figures/influence_cooks_*.png\n")
