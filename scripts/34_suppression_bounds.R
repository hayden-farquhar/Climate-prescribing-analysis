# ==============================================================================
# Script 34: Suppression Bounded Imputation Sensitivity Analysis
# ==============================================================================
# Addresses non-random missingness from AIHW cell suppression (cells < 5).
# Approximately 21% of SA4-week observations are suppressed, concentrated
# in remote/low-population areas.
#
# Approach: probabilistic sensitivity analysis using extreme bounds.
# Since suppressed cells are strictly bounded [0, 4], we create:
#   1. Lower bound dataset: all suppressed cells = 0
#   2. Upper bound dataset: all suppressed cells = 4
#   3. Midpoint dataset: all suppressed cells = 2
# and refit primary DLNMs under each scenario.
#
# If estimates are stable across bounds, suppression bias is negligible.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv
#
# Outputs:
#   outputs/tables/suppression_bounds_rr.csv
#   outputs/tables/suppression_bounds_equity.csv
#   outputs/figures/suppression_bounds_forest.png
#   outputs/figures/suppression_bounds_equity.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate",
#                       "ggplot2", "patchwork", "tidyr"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
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
cat("Script 34: Suppression Bounded Imputation Sensitivity\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA AND CHARACTERISE SUPPRESSION
# ==============================================================================

cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

panel <- panel |>
  mutate(
    year         = year(week_start),
    week_of_year = isoweek(week_start),
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52),
    is_suppressed = (suppressed == "True" | suppressed == TRUE)
  )

# Suppression summary by group
cat("\nSuppression rates by medication group:\n")
panel |>
  group_by(analysis_group) |>
  summarise(
    total = n(),
    n_suppressed = sum(is_suppressed, na.rm = TRUE),
    pct_suppressed = round(100 * mean(is_suppressed, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  filter(analysis_group %in% primary_groups) |>
  print()

# Load SEIFA for equity-specific suppression analysis
seifa <- read.csv(file.path(ref_dir, "seifa_sa4_quintiles.csv"),
                  stringsAsFactors = FALSE)
seifa$sa4_code <- as.character(seifa$sa4_code)

panel <- panel |>
  left_join(seifa |> select(sa4_code, irsd_quintile), by = "sa4_code")

cat("\nSuppression rates by SEIFA quintile (resp_total):\n")
panel |>
  filter(analysis_group == "resp_total", !is.na(irsd_quintile)) |>
  group_by(irsd_quintile) |>
  summarise(
    n = n(),
    pct_suppressed = round(100 * mean(is_suppressed, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  print()


# ==============================================================================
# 2. CREATE BOUNDED DATASETS
# ==============================================================================

cat("\nCreating bounded imputation datasets...\n")

# Lower bound: suppressed cells = 0
panel_lower <- panel |>
  mutate(
    count_imputed = ifelse(is_suppressed, 0, count),
    bound = "lower"
  )

# Upper bound: suppressed cells = 4
panel_upper <- panel |>
  mutate(
    count_imputed = ifelse(is_suppressed, 4, count),
    bound = "upper"
  )

# Midpoint: suppressed cells = 2
panel_mid <- panel |>
  mutate(
    count_imputed = ifelse(is_suppressed, 2, count),
    bound = "midpoint"
  )

# Original (excludes suppressed)
panel_orig <- panel |>
  mutate(
    count_imputed = count,
    bound = "original"
  )

cat("  Lower bound (0):", sum(panel_lower$is_suppressed & panel_lower$analysis_group == "resp_total"),
    "cells imputed\n")
cat("  Upper bound (4):", sum(panel_upper$is_suppressed & panel_upper$analysis_group == "resp_total"),
    "cells imputed\n")


# ==============================================================================
# 3. FIT DLNM UNDER EACH BOUND
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models under each imputation scenario\n")
cat("(3 groups x 4 bounds = 12 models)\n")
cat("=" |> strrep(70), "\n")

fit_bounded_dlnm <- function(data, group_name, bound_name) {
  # For original, exclude suppressed; for bounds, include all with imputed count
  if (bound_name == "original") {
    df <- data |>
      filter(analysis_group == group_name,
             !is_suppressed,
             !is.na(count_imputed),
             !is.na(tmax_mean)) |>
      arrange(sa4_code, week_start)
  } else {
    df <- data |>
      filter(analysis_group == group_name,
             !is.na(count_imputed),
             !is.na(tmax_mean)) |>
      arrange(sa4_code, week_start)
  }

  if (nrow(df) < 500) return(NULL)

  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years  <- n_distinct(year(df$week_start))
  trend_df <- max(2, round(n_years * 2))
  temp_ref <- median(df$tmax_mean, na.rm = TRUE)

  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  model <- glm(
    count_imputed ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  pred <- crosspred(cb, model, at = c(temp_p95, temp_p99), cen = temp_ref)

  data.frame(
    group = group_name,
    label = group_labels[group_name],
    bound = bound_name,
    rr_p95    = round(pred$allRRfit[1], 4),
    rr_p95_lo = round(pred$allRRlow[1], 4),
    rr_p95_hi = round(pred$allRRhigh[1], 4),
    rr_p99    = round(pred$allRRfit[2], 4),
    rr_p99_lo = round(pred$allRRlow[2], 4),
    rr_p99_hi = round(pred$allRRhigh[2], 4),
    n_obs      = nrow(df),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )
}

all_results <- list()
datasets <- list(original = panel_orig, lower = panel_lower,
                 midpoint = panel_mid, upper = panel_upper)

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  for (bnd in c("original", "lower", "midpoint", "upper")) {
    cat("  ", bnd, "... ")
    res <- fit_bounded_dlnm(datasets[[bnd]], grp, bnd)
    if (!is.null(res)) {
      all_results[[paste0(grp, "_", bnd)]] <- res
      cat("RR@p95 =", res$rr_p95,
          "(", res$rr_p95_lo, "-", res$rr_p95_hi, ")  n =",
          format(res$n_obs, big.mark = ","), "\n")
    }
  }
}


# ==============================================================================
# 4. EQUITY-STRATIFIED BOUNDS (Q1 most disadvantaged)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Equity-stratified bounded imputation (Q1 vs Q5)\n")
cat("=" |> strrep(70), "\n")

equity_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  for (quintile in c(1, 5)) {
    qlabel <- ifelse(quintile == 1, "Q1 (most disadvantaged)", "Q5 (least disadvantaged)")
    cat("  ", qlabel, "\n")

    for (bnd in c("original", "lower", "upper")) {
      ds <- datasets[[bnd]]

      if (bnd == "original") {
        df <- ds |>
          filter(analysis_group == grp,
                 irsd_quintile == quintile,
                 !is_suppressed,
                 !is.na(count_imputed),
                 !is.na(tmax_mean)) |>
          arrange(sa4_code, week_start)
      } else {
        df <- ds |>
          filter(analysis_group == grp,
                 irsd_quintile == quintile,
                 !is.na(count_imputed),
                 !is.na(tmax_mean)) |>
          arrange(sa4_code, week_start)
      }

      if (nrow(df) < 300) {
        cat("    ", bnd, ": too few obs\n")
        next
      }

      temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
      cb <- crossbasis(
        df$tmax_mean,
        lag = max_lag,
        argvar = list(fun = "ns", knots = temp_knots),
        arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
      )

      n_years  <- n_distinct(year(df$week_start))
      trend_df <- max(2, round(n_years * 2))
      temp_ref <- median(df$tmax_mean, na.rm = TRUE)
      temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)

      model <- glm(
        count_imputed ~ cb +
          ns(time_index, df = trend_df) +
          sin1 + cos1 + sin2 + cos2 +
          precip_total +
          factor(sa4_code),
        data = df,
        family = quasipoisson(link = "log")
      )

      pred <- crosspred(cb, model, at = temp_p95, cen = temp_ref)

      equity_results[[paste0(grp, "_Q", quintile, "_", bnd)]] <- data.frame(
        group    = grp,
        label    = group_labels[grp],
        quintile = quintile,
        q_label  = qlabel,
        bound    = bnd,
        rr_p95    = round(pred$allRRfit[1], 4),
        rr_p95_lo = round(pred$allRRlow[1], 4),
        rr_p95_hi = round(pred$allRRhigh[1], 4),
        n_obs    = nrow(df),
        stringsAsFactors = FALSE
      )

      cat("    ", bnd, ": RR@p95 =", round(pred$allRRfit[1], 4), "\n")
    }
  }
}


# ==============================================================================
# 5. SAVE RESULTS
# ==============================================================================

bounds_table <- bind_rows(all_results)
write.csv(bounds_table,
          file.path(tab_dir, "suppression_bounds_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/suppression_bounds_rr.csv\n")

equity_table <- bind_rows(equity_results)
write.csv(equity_table,
          file.path(tab_dir, "suppression_bounds_equity.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/suppression_bounds_equity.csv\n")


# ==============================================================================
# 6. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# --- Forest plot: bounds comparison ---
bounds_table$bound <- factor(bounds_table$bound,
                              levels = c("original", "lower", "midpoint", "upper"))
bounds_table$bound_label <- c(
  original = "Original (exclude suppressed)",
  lower    = "Lower bound (impute 0)",
  midpoint = "Midpoint (impute 2)",
  upper    = "Upper bound (impute 4)"
)[as.character(bounds_table$bound)]

bounds_table$y_label <- paste0(bounds_table$label, "\n", bounds_table$bound_label)

p_forest <- ggplot(bounds_table, aes(x = rr_p95, y = y_label, colour = bound)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c(original = "#4575b4", lower = "#91bfdb",
               midpoint = "#fee090", upper = "#d73027"),
    name = "Imputation",
    labels = c("Original", "Lower (0)", "Midpoint (2)", "Upper (4)")
  ) +
  labs(
    title = "Suppression bounded imputation sensitivity",
    subtitle = "Cumulative RR at 95th percentile under extreme imputation scenarios",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "suppression_bounds_forest.png"),
       p_forest, width = 10, height = 10, dpi = 300)
cat("  -> Saved: suppression_bounds_forest.png\n")

# --- Equity bounds plot ---
if (nrow(equity_table) > 0) {
  equity_table$bound <- factor(equity_table$bound,
                                levels = c("original", "lower", "upper"))
  equity_table$y_label <- paste0(equity_table$label, " (", equity_table$q_label,
                                  ")\n", c(original = "Original",
                                           lower = "Lower (0)",
                                           upper = "Upper (4)")[as.character(equity_table$bound)])

  p_equity <- ggplot(equity_table, aes(x = rr_p95, y = y_label,
                                        colour = interaction(bound, quintile))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.6) +
    geom_point(size = 2.5) +
    labs(
      title = "Suppression bounds: equity stratification (Q1 vs Q5)",
      subtitle = "Testing whether suppression disproportionately affects disadvantaged quintile results",
      x = "Cumulative RR (95% CI)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "none"
    )

  ggsave(file.path(fig_dir, "suppression_bounds_equity.png"),
         p_equity, width = 11, height = 12, dpi = 300)
  cat("  -> Saved: suppression_bounds_equity.png\n")
}


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 34 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Bounded imputation results (RR at 95th percentile):\n")
bounds_table |>
  select(label, bound, rr_p95, rr_p95_lo, rr_p95_hi, n_obs) |>
  print(row.names = FALSE)

cat("\nInterpretation:\n")
cat("  - If RR stable across bounds: suppression bias is negligible\n")
cat("  - If RR shifts substantially: suppression bias is material and\n")
cat("    confidence in equity/remoteness stratifications is weakened\n")
cat("  - Upper bound (4) adds more dispensing to low-volume SA4s,\n")
cat("    potentially attenuating effects if suppressed areas differ systematically\n")

cat("\nOutputs:\n")
cat("  outputs/tables/suppression_bounds_rr.csv\n")
cat("  outputs/tables/suppression_bounds_equity.csv\n")
cat("  outputs/figures/suppression_bounds_forest.png\n")
cat("  outputs/figures/suppression_bounds_equity.png\n")
