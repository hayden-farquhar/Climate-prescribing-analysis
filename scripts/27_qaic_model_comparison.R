# ==============================================================================
# Script 27: Quasi-AIC Model Comparison
# ==============================================================================
# Systematic comparison of model specifications using quasi-AIC to identify
# the best-fitting model and verify that the primary specification is
# well-supported relative to alternatives.
#
# Compares:
#   1. Primary model (full specification)
#   2. Without long-term trend
#   3. Without seasonality (Fourier terms)
#   4. Without precipitation
#   5. Without SA4 fixed effects
#   6. Without cross-basis (temperature effect)
#   7. Simplified trend (linear only)
#   8. Additional Fourier harmonics (3 pairs)
#   9. Population offset model
#
# QAIC = -2 * loglik / dispersion + 2 * k
# (Burnham & Anderson 2002; adapted for quasi-likelihood)
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/qaic_model_comparison.csv
#   outputs/figures/qaic_comparison.png
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
cat("Script 27: Quasi-AIC Model Comparison\n")
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
# 2. QAIC FUNCTION
# ==============================================================================

compute_qaic <- function(model) {
  ll <- -0.5 * deviance(model)
  disp <- summary(model)$dispersion
  k <- length(coef(model))
  qaic <- -2 * ll / disp + 2 * k
  return(list(qaic = qaic, k = k, dispersion = disp,
              deviance = deviance(model)))
}


# ==============================================================================
# 3. FIT MODELS AND COMPARE
# ==============================================================================

all_results <- list()

for (grp in primary_groups) {
  cat("\n", "=" |> strrep(50), "\n")
  cat(group_labels[grp], "\n")
  cat("=" |> strrep(50), "\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  n <- nrow(df)
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))
  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)

  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  # --- Model specifications ---
  models <- list()

  # 1. Primary (full)
  cat("  1. Primary model...\n")
  models[["Primary (full)"]] <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 2. Without long-term trend
  cat("  2. Without trend...\n")
  models[["No trend"]] <- glm(
    count ~ cb +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 3. Without seasonality
  cat("  3. Without seasonality...\n")
  models[["No seasonality"]] <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 4. Without precipitation
  cat("  4. Without precipitation...\n")
  models[["No precipitation"]] <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 5. Without SA4 fixed effects
  cat("  5. Without SA4 FE...\n")
  models[["No SA4 FE"]] <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total,
    data = df, family = quasipoisson(link = "log")
  )

  # 6. Without cross-basis (null temperature model)
  cat("  6. Null (no temperature)...\n")
  models[["No temperature"]] <- glm(
    count ~
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 7. Linear trend only
  cat("  7. Linear trend...\n")
  models[["Linear trend"]] <- glm(
    count ~ cb +
      time_index +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 8. Extra Fourier harmonics (3 pairs)
  cat("  8. Extra harmonics...\n")
  models[["3 Fourier pairs"]] <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 + sin3 + cos3 +
      precip_total +
      factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  # 9. With log(count) population offset (use crude_rate denominator approach)
  # Approximate population as count/crude_rate * 100000
  cat("  9. Rate model (offset)...\n")
  if ("crude_rate" %in% names(df) && any(!is.na(df$crude_rate) & df$crude_rate > 0)) {
    df_rate <- df |> filter(!is.na(crude_rate) & crude_rate > 0)
    df_rate$pop_approx <- df_rate$count / df_rate$crude_rate * 100000
    df_rate$log_pop <- log(df_rate$pop_approx)

    cb_rate <- crossbasis(
      df_rate$tmax_mean,
      lag = max_lag,
      argvar = list(fun = "ns", knots = temp_knots),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    )

    models[["With pop offset"]] <- tryCatch(
      glm(
        count ~ cb_rate +
          ns(time_index, df = trend_df) +
          sin1 + cos1 + sin2 + cos2 +
          precip_total +
          factor(sa4_code) +
          offset(log_pop),
        data = df_rate, family = quasipoisson(link = "log")
      ),
      error = function(e) { cat("    Offset model failed:", e$message, "\n"); NULL }
    )
  } else {
    cat("    No crude_rate available, skipping offset model.\n")
  }

  # --- Compute QAIC for all models ---
  cat("\n  Results:\n")

  # Use primary model dispersion for all QAIC calculations (standard practice)
  disp_primary <- summary(models[["Primary (full)"]])$dispersion

  for (mod_name in names(models)) {
    mod <- models[[mod_name]]
    if (is.null(mod)) next

    k <- length(coef(mod))
    dev <- deviance(mod)
    disp_own <- summary(mod)$dispersion
    qaic <- -2 * (-0.5 * dev) / disp_primary + 2 * k

    # Get RR at p95 if model includes temperature
    rr_p95 <- rr_lo <- rr_hi <- NA
    if (mod_name != "No temperature") {
      cb_use <- if (mod_name == "With pop offset") cb_rate else cb
      pred <- tryCatch(
        crosspred(cb_use, mod, at = temp_p95, cen = temp_ref),
        error = function(e) NULL
      )
      if (!is.null(pred)) {
        rr_p95 <- round(pred$allRRfit[1], 4)
        rr_lo  <- round(pred$allRRlow[1], 4)
        rr_hi  <- round(pred$allRRhigh[1], 4)
      }
    }

    all_results[[paste0(grp, "_", mod_name)]] <- data.frame(
      group      = grp,
      label      = group_labels[grp],
      model      = mod_name,
      n_params   = k,
      deviance   = round(dev, 1),
      dispersion = round(disp_own, 2),
      qaic       = round(qaic, 1),
      rr_p95     = rr_p95,
      rr_p95_lo  = rr_lo,
      rr_p95_hi  = rr_hi,
      stringsAsFactors = FALSE
    )

    cat("    ", mod_name, ": QAIC =", round(qaic, 1),
        " k =", k, " RR =", ifelse(is.na(rr_p95), "—", rr_p95), "\n")
  }
}


# ==============================================================================
# 4. COMPUTE DELTA-QAIC AND WEIGHTS
# ==============================================================================

results_table <- bind_rows(all_results)

# Delta QAIC (relative to best model within each group)
results_table <- results_table |>
  group_by(group) |>
  mutate(
    delta_qaic = qaic - min(qaic, na.rm = TRUE),
    # Akaike weights (approximate)
    weight = exp(-0.5 * delta_qaic) / sum(exp(-0.5 * delta_qaic), na.rm = TRUE)
  ) |>
  ungroup()

results_table$delta_qaic <- round(results_table$delta_qaic, 1)
results_table$weight     <- round(results_table$weight, 4)


# ==============================================================================
# 5. SAVE RESULTS
# ==============================================================================

write.csv(as.data.frame(results_table),
          file.path(tab_dir, "qaic_model_comparison.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/qaic_model_comparison.csv\n")


# ==============================================================================
# 6. VISUALISATION
# ==============================================================================

cat("\nGenerating QAIC comparison plots...\n")

p_qaic <- ggplot(results_table,
                  aes(x = delta_qaic,
                      y = reorder(model, -delta_qaic))) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, colour = "grey20") +
  geom_vline(xintercept = 2, linetype = "dashed", colour = "#d73027",
             linewidth = 0.5) +
  annotate("text", x = 2.5, y = 1, label = "Delta = 2",
           colour = "#d73027", hjust = 0, size = 3) +
  facet_wrap(~ label, scales = "free_x", ncol = 1) +
  labs(
    title = "Model comparison: Delta QAIC",
    subtitle = "Lower is better; models within 2 units are considered equivalent",
    x = "Delta QAIC (relative to best model)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "qaic_comparison.png"),
       p_qaic, width = 10, height = 10, dpi = 300)
cat("  -> Saved: qaic_comparison.png\n")


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 27 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("QAIC model comparison:\n")
for (grp in primary_groups) {
  cat("\n", group_labels[grp], ":\n")
  grp_data <- results_table |>
    filter(group == grp) |>
    arrange(delta_qaic) |>
    select(model, n_params, qaic, delta_qaic, weight, rr_p95)
  print(as.data.frame(grp_data), row.names = FALSE)
}

cat("\nInterpretation:\n")
cat("  - Best model has delta QAIC = 0\n")
cat("  - Models within delta < 2 have substantial support\n")
cat("  - Models with delta > 10 have essentially no support\n")
cat("  - If primary model has delta ~ 0, specification is well-supported\n")
cat("  - Weight approximates probability that model is best in the set\n")

cat("\nOutputs:\n")
cat("  outputs/tables/qaic_model_comparison.csv\n")
cat("  outputs/figures/qaic_comparison.png\n")
