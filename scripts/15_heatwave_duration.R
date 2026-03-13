# ==============================================================================
# Script 15: Heatwave Duration Analysis
# ==============================================================================
# Tests whether multi-day heatwaves have larger effects on prescribing than
# single hot days, using the heatwave_days variable from the panel.
#
# Approach:
#   1. Categorise weeks by heatwave intensity:
#      - No heatwave days (0 days above p90)
#      - Mild (1-2 heatwave days per week)
#      - Moderate (3-4 heatwave days)
#      - Severe (5-7 heatwave days, i.e. most/all of the week)
#   2. Fit a quasi-Poisson model with heatwave category as exposure
#   3. Compare dose-response across duration categories
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv
#
# Outputs:
#   outputs/tables/heatwave_duration_rr.csv
#   outputs/figures/heatwave_duration_bar.png
#   outputs/figures/heatwave_duration_panel.png
#
# Dependencies:
#   install.packages(c("dplyr", "lubridate", "ggplot2", "patchwork", "splines"))
# ==============================================================================

library(dplyr)
library(lubridate)
library(splines)
library(ggplot2)
library(patchwork)

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

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
cat("Script 15: Heatwave Duration Analysis\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================

cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

# Classify heatwave duration categories
panel <- panel |>
  mutate(
    # days_above_p90 counts days in the week with Tmax > 90th percentile
    hw_category = case_when(
      is.na(days_above_p90) ~ NA_character_,
      days_above_p90 == 0   ~ "None (0 days)",
      days_above_p90 <= 2   ~ "Mild (1-2 days)",
      days_above_p90 <= 4   ~ "Moderate (3-4 days)",
      days_above_p90 >= 5   ~ "Severe (5-7 days)"
    ),
    hw_category = factor(hw_category,
                         levels = c("None (0 days)", "Mild (1-2 days)",
                                    "Moderate (3-4 days)", "Severe (5-7 days)")),
    # Time variables
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    week_of_year = isoweek(week_start),
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )

cat("\n  Heatwave category distribution (all rows):\n")
print(table(panel$hw_category, useNA = "ifany"))


# ==============================================================================
# 2. FIT HEATWAVE DURATION MODELS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting heatwave duration models\n")
cat("=" |> strrep(70), "\n")

all_results <- list()

for (grp in outcome_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(hw_category),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")
  cat("  Category distribution:\n")
  print(table(df$hw_category))

  if (nrow(df) < 500) {
    cat("  Too few observations, skipping.\n")
    next
  }

  n_years  <- n_distinct(year(df$week_start))
  trend_df <- max(2, round(n_years * 2))

  # Model: heatwave category as factor (reference = "None")
  # Also control for mean temperature to isolate duration effect
  model <- glm(
    count ~ hw_category +
      tmax_mean +        # control for temperature level
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  coef_table <- summary(model)$coefficients
  cat("  Dispersion:", round(summary(model)$dispersion, 2), "\n")

  # Extract heatwave category RRs
  hw_terms <- grep("^hw_category", rownames(coef_table), value = TRUE)

  results_rows <- data.frame(
    group    = grp,
    label    = group_labels[grp],
    category = "None (0 days)",
    rr       = 1,
    rr_lo    = 1,
    rr_hi    = 1,
    p_value  = NA,
    stringsAsFactors = FALSE
  )

  for (term in hw_terms) {
    cat_name <- gsub("hw_category", "", term)
    est  <- coef_table[term, "Estimate"]
    se   <- coef_table[term, "Std. Error"]
    pval <- coef_table[term, "Pr(>|t|)"]
    rr   <- exp(est)
    rr_lo <- exp(est - 1.96 * se)
    rr_hi <- exp(est + 1.96 * se)

    results_rows <- bind_rows(results_rows, data.frame(
      group    = grp,
      label    = group_labels[grp],
      category = cat_name,
      rr       = round(rr, 4),
      rr_lo    = round(rr_lo, 4),
      rr_hi    = round(rr_hi, 4),
      p_value  = pval,
      stringsAsFactors = FALSE
    ))

    cat("  ", cat_name, ": RR =", round(rr, 4),
        "(", round(rr_lo, 4), "-", round(rr_hi, 4), ")",
        "p =", format.pval(pval, digits = 3), "\n")
  }

  results_rows$n_obs <- nrow(df)
  all_results[[grp]] <- results_rows
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

results_table <- bind_rows(all_results)
write.csv(results_table,
          file.path(tab_dir, "heatwave_duration_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/heatwave_duration_rr.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating heatwave duration plots...\n")

hw_colours <- c(
  "None (0 days)"     = "#91bfdb",
  "Mild (1-2 days)"   = "#fee090",
  "Moderate (3-4 days)" = "#fc8d59",
  "Severe (5-7 days)" = "#d73027"
)

# --- Bar chart for each medication group ---
plot_hw_duration <- function(grp) {
  df <- results_table |>
    filter(group == grp) |>
    mutate(category = factor(category,
                             levels = c("None (0 days)", "Mild (1-2 days)",
                                        "Moderate (3-4 days)", "Severe (5-7 days)")))

  p <- ggplot(df, aes(x = category, y = rr, fill = category)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = rr_lo, ymax = rr_hi),
                  width = 0.2, linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = hw_colours, guide = "none") +
    labs(
      title = group_labels[grp],
      x = "Heatwave duration (days above p90 per week)",
      y = "Rate Ratio (vs no heatwave)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )

  return(p)
}

# Individual plots
for (grp in outcome_groups) {
  if (grp %in% unique(results_table$group)) {
    p <- plot_hw_duration(grp)
    ggsave(file.path(fig_dir, paste0("heatwave_duration_", grp, ".png")),
           p, width = 7, height = 5, dpi = 300)
  }
}
cat("  -> Saved: heatwave_duration_*.png\n")

# Panel of primary groups
primary <- intersect(c("cvd_total", "mh_total", "resp_total"),
                     unique(results_table$group))
if (length(primary) >= 2) {
  plots <- lapply(primary, plot_hw_duration)
  combined <- wrap_plots(plots, ncol = 3)
  ggsave(file.path(fig_dir, "heatwave_duration_panel.png"),
         combined, width = 16, height = 5, dpi = 300)
  cat("  -> Saved: heatwave_duration_panel.png\n")
}

# --- Dose-response: days_above_p90 as continuous ---
cat("Generating continuous dose-response...\n")

# Use poly() instead of ns() — days_above_p90 is heavily zero-inflated
# so ns() cannot place interior knots. Quadratic polynomial is sufficient
# to capture non-linearity across 0-7 days.
dose_plots <- list()
for (grp in c("cvd_total", "mh_total", "resp_total")) {
  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(days_above_p90),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  n_years  <- n_distinct(year(df$week_start))
  trend_df <- max(2, round(n_years * 2))

  model <- glm(
    count ~ poly(days_above_p90, 2, raw = TRUE) +
      tmax_mean +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  # Predict at each integer value 0-7
  newdata_base <- df[1, , drop = FALSE]
  preds <- list()
  for (d in 0:7) {
    nd <- newdata_base
    nd$days_above_p90 <- d
    p <- predict(model, newdata = nd, type = "link", se.fit = TRUE)
    preds[[d + 1]] <- data.frame(days = d, fit = p$fit, se = p$se.fit)
  }
  pred_df <- bind_rows(preds)

  # Centre on 0 days
  ref_fit <- pred_df$fit[pred_df$days == 0]
  pred_df <- pred_df |>
    mutate(
      rr    = exp(fit - ref_fit),
      rr_lo = exp((fit - ref_fit) - 1.96 * se),
      rr_hi = exp((fit - ref_fit) + 1.96 * se)
    )

  p <- ggplot(pred_df, aes(x = days, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "#fc8d59", alpha = 0.2) +
    geom_line(colour = "#d73027", linewidth = 0.8) +
    geom_point(colour = "#d73027", size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(breaks = 0:7) +
    labs(
      title = group_labels[grp],
      x = "Days above 90th percentile (per week)",
      y = "RR (vs 0 heatwave days)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  dose_plots[[grp]] <- p
}

if (length(dose_plots) >= 2) {
  combined <- wrap_plots(dose_plots, ncol = 3)
  ggsave(file.path(fig_dir, "heatwave_dose_response.png"),
         combined, width = 15, height = 5, dpi = 300)
  cat("  -> Saved: heatwave_dose_response.png\n")
}


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 15 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Heatwave duration RRs (vs no heatwave):\n")
print(
  results_table |>
    filter(category != "None (0 days)") |>
    select(label, category, rr, rr_lo, rr_hi, p_value) |>
    as.data.frame(),
  row.names = FALSE
)

cat("\nOutputs:\n")
cat("  outputs/tables/heatwave_duration_rr.csv\n")
cat("  outputs/figures/heatwave_duration_*.png\n")
cat("  outputs/figures/heatwave_duration_panel.png\n")
cat("  outputs/figures/heatwave_dose_response.png\n")
