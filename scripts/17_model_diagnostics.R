# ==============================================================================
# Script 17: Model Diagnostics
# ==============================================================================
# Performs residual diagnostics on the primary DLNM models:
#   - Overdispersion test and summary
#   - Deviance residuals vs fitted values
#   - Deviance residuals vs time (check for residual trend)
#   - Deviance residuals vs temperature (check for residual non-linearity)
#   - Autocorrelation function (ACF/PACF) of residuals by SA4
#   - Quantile-quantile plots
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/model_diagnostics_summary.csv
#   outputs/figures/diagnostics_residuals_*.png
#   outputs/figures/diagnostics_acf_*.png
#   outputs/figures/diagnostics_qq_*.png
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

outcome_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

cat("=" |> strrep(70), "\n")
cat("Script 17: Model Diagnostics\n")
cat("=" |> strrep(70), "\n\n")

# ==============================================================================
# 1. LOAD DATA AND FIT MODELS
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

diag_results <- list()

for (grp in outcome_groups) {
  cat("\n", "=" |> strrep(50), "\n")
  cat("Diagnostics for:", group_labels[grp], "\n")
  cat("=" |> strrep(50), "\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  # Fit model (same specification as script 05)
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

  # --- Overdispersion ---
  disp <- summary(model)$dispersion
  resid_deviance <- deviance(model)
  resid_df <- model$df.residual
  disp_ratio <- resid_deviance / resid_df
  # Formal test: if Poisson, deviance ~ chi-squared(df)
  disp_pval <- pchisq(resid_deviance, resid_df, lower.tail = FALSE)

  cat("  Overdispersion parameter:", round(disp, 2), "\n")
  cat("  Deviance/df ratio:", round(disp_ratio, 2), "\n")
  cat("  Chi-sq test for overdispersion: p =", format.pval(disp_pval), "\n")

  # --- Deviance residuals ---
  # glm drops rows with NA (from crossbasis lags, precip, etc.)
  # Use na.action to subset df to only the rows the model actually used
  if (!is.null(model$na.action)) {
    df <- df[-model$na.action, ]
  }

  dev_resid <- residuals(model, type = "deviance")
  fitted_vals <- fitted(model)
  cat("  Model used", length(dev_resid), "of", nrow(df) + length(model$na.action), "rows\n")

  df$dev_resid  <- as.numeric(dev_resid)
  df$fitted_val <- as.numeric(fitted_vals)

  # --- 1. Residuals vs Fitted ---
  p1 <- ggplot(df, aes(x = log(fitted_val), y = dev_resid)) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.3, colour = "steelblue",
                se = FALSE, linewidth = 0.8) +
    labs(title = paste0(group_labels[grp], " — Residuals vs Fitted"),
         x = "log(Fitted values)", y = "Deviance residuals") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # --- 2. Residuals vs Time ---
  resid_time <- df |>
    group_by(week_start) |>
    summarise(mean_resid = mean(dev_resid, na.rm = TRUE), .groups = "drop")

  p2 <- ggplot(resid_time, aes(x = week_start, y = mean_resid)) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    geom_smooth(method = "loess", span = 0.2, colour = "steelblue",
                se = FALSE, linewidth = 0.8) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    labs(title = paste0(group_labels[grp], " — Mean residuals over time"),
         x = NULL, y = "Mean deviance residual") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # --- 3. Residuals vs Temperature ---
  p3 <- ggplot(df, aes(x = tmax_mean, y = dev_resid)) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.3, colour = "steelblue",
                se = FALSE, linewidth = 0.8) +
    labs(title = paste0(group_labels[grp], " — Residuals vs Temperature"),
         x = "Weekly mean Tmax (°C)", y = "Deviance residuals") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  # --- 4. QQ Plot ---
  qq_data <- data.frame(
    theoretical = qnorm(ppoints(length(dev_resid))),
    sample      = sort(dev_resid)
  )

  p4 <- ggplot(qq_data, aes(x = theoretical, y = sample)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
    labs(title = paste0(group_labels[grp], " — Q-Q plot"),
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  combined <- (p1 | p2) / (p3 | p4)
  ggsave(file.path(fig_dir, paste0("diagnostics_residuals_", grp, ".png")),
         combined, width = 14, height = 10, dpi = 300)
  cat("  -> Saved: diagnostics_residuals_", grp, ".png\n")

  # --- 5. ACF of residuals (sample of SA4s) ---
  cat("  Computing ACF...\n")
  sa4_sample <- sample(unique(df$sa4_code), min(9, n_distinct(df$sa4_code)))

  acf_plots <- list()
  mean_acf_lag1 <- c()

  for (i in seq_along(sa4_sample)) {
    sa4 <- sa4_sample[i]
    resid_sa4 <- df$dev_resid[df$sa4_code == sa4]
    if (length(resid_sa4) < 20) next

    acf_obj <- acf(resid_sa4, lag.max = 20, plot = FALSE)
    acf_df <- data.frame(
      lag = acf_obj$lag[-1],
      acf = acf_obj$acf[-1]
    )
    mean_acf_lag1 <- c(mean_acf_lag1, acf_obj$acf[2])  # lag-1 ACF

    acf_plots[[i]] <- ggplot(acf_df, aes(x = lag, y = acf)) +
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(length(resid_sa4)),
                 linetype = "dashed", colour = "blue", alpha = 0.5) +
      geom_segment(aes(xend = lag, yend = 0), linewidth = 0.6) +
      labs(title = sa4, x = "Lag", y = "ACF") +
      theme_minimal(base_size = 9) +
      theme(plot.title = element_text(size = 9))
  }

  if (length(acf_plots) >= 4) {
    acf_combined <- wrap_plots(acf_plots[1:min(9, length(acf_plots))],
                                ncol = 3) +
      plot_annotation(
        title = paste0(group_labels[grp],
                       " — ACF of deviance residuals (sample of SA4s)"),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )
    ggsave(file.path(fig_dir, paste0("diagnostics_acf_", grp, ".png")),
           acf_combined, width = 12, height = 10, dpi = 300)
    cat("  -> Saved: diagnostics_acf_", grp, ".png\n")
  }

  # --- Compute ACF lag-1 across ALL SA4s ---
  all_sa4s <- unique(df$sa4_code)
  all_acf1 <- sapply(all_sa4s, function(sa4) {
    r <- df$dev_resid[df$sa4_code == sa4]
    if (length(r) < 20) return(NA)
    acf(r, lag.max = 1, plot = FALSE)$acf[2]
  })

  cat("  Mean lag-1 ACF across all SA4s:", round(mean(all_acf1, na.rm = TRUE), 3), "\n")
  cat("  Median lag-1 ACF:", round(median(all_acf1, na.rm = TRUE), 3), "\n")
  cat("  Range:", round(min(all_acf1, na.rm = TRUE), 3), "to",
      round(max(all_acf1, na.rm = TRUE), 3), "\n")

  # --- Store diagnostics summary ---
  diag_results[[grp]] <- data.frame(
    group             = grp,
    label             = group_labels[grp],
    n_obs             = nrow(df),
    dispersion        = round(disp, 2),
    deviance_df_ratio = round(disp_ratio, 4),
    overdispersion_p  = disp_pval,
    mean_resid        = round(mean(dev_resid), 4),
    sd_resid          = round(sd(dev_resid), 4),
    skewness          = round(moments::skewness(dev_resid), 4) |>
                          tryCatch(error = function(e)
                            round((mean(dev_resid^3)) / (sd(dev_resid)^3), 4)),
    mean_acf_lag1     = round(mean(all_acf1, na.rm = TRUE), 4),
    median_acf_lag1   = round(median(all_acf1, na.rm = TRUE), 4),
    max_acf_lag1      = round(max(all_acf1, na.rm = TRUE), 4),
    stringsAsFactors  = FALSE
  )
}


# ==============================================================================
# 2. SAVE RESULTS
# ==============================================================================

diag_table <- bind_rows(diag_results)
write.csv(diag_table,
          file.path(tab_dir, "model_diagnostics_summary.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/model_diagnostics_summary.csv\n")


# ==============================================================================
# 3. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 17 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Model diagnostics summary:\n")
print(as.data.frame(diag_table), row.names = FALSE)

cat("\nInterpretation guide:\n")
cat("  - Dispersion >> 1: overdispersion present, quasi-Poisson appropriate\n")
cat("  - Mean ACF lag-1 near 0: minimal residual autocorrelation\n")
cat("  - Mean ACF lag-1 > 0.1: consider Newey-West SEs or AR structure\n")
cat("  - Residual plots should show no systematic patterns\n")

cat("\nOutputs:\n")
cat("  outputs/tables/model_diagnostics_summary.csv\n")
cat("  outputs/figures/diagnostics_residuals_*.png\n")
cat("  outputs/figures/diagnostics_acf_*.png\n")
