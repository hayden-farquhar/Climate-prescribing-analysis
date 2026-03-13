# ==============================================================================
# Script 35: Demand vs Supply Triangulation
# ==============================================================================
# Addresses the concern that dispensing data conflates patient demand (acute
# symptoms) with supply-side artifacts (GP review schedules, stockpiling).
#
# Approach:
#   1. Compare SABA relievers (acute rescue proxy) vs preventer medications
#      (scheduled refills) — relievers should show stronger/different temperature
#      response if the signal reflects genuine symptom exacerbation
#   2. Examine within-week dispensing velocity: acute demand should show
#      faster lag response than supply-driven dispensing
#   3. Quantify the reliever/preventer RR ratio as evidence of demand-driven
#      signal specificity
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/demand_supply_comparison.csv
#   outputs/tables/demand_supply_lag_profiles.csv
#   outputs/figures/demand_supply_forest.png
#   outputs/figures/demand_supply_lag_comparison.png
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
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

# Compare demand-sensitive vs supply-driven medications
comparison_groups <- list(
  respiratory = list(
    demand = "resp_relievers",      # SABAs — acute rescue medication
    supply = "resp_preventers",     # ICS — scheduled refills
    label  = "Respiratory"
  ),
  mental_health = list(
    demand = "mh_anxiolytics",      # PRN anxiolytics — acute distress
    supply = "mh_antidepressants",  # Daily SSRIs — scheduled refills
    label  = "Mental Health"
  )
)

group_labels <- c(
  resp_relievers     = "Resp. Relievers (SABA)",
  resp_preventers    = "Resp. Preventers (ICS)",
  mh_anxiolytics     = "Anxiolytics (PRN)",
  mh_antidepressants = "Antidepressants (daily)",
  cvd_total          = "Cardiovascular (all)"
)

cat("=" |> strrep(70), "\n")
cat("Script 35: Demand vs Supply Triangulation\n")
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
    year         = year(week_start),
    week_of_year = isoweek(week_start),
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )


# ==============================================================================
# 2. FIT DLNM FOR EACH SUBGROUP
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models for demand vs supply medication subgroups\n")
cat("=" |> strrep(70), "\n")

all_groups <- c("resp_relievers", "resp_preventers",
                "mh_anxiolytics", "mh_antidepressants",
                "cvd_total")

all_results <- list()
all_preds   <- list()

for (grp in all_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 500) {
    cat("  Too few observations, skipping.\n")
    next
  }

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

  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  model <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  pred <- crosspred(cb, model,
                     at = seq(min(df$tmax_mean), max(df$tmax_mean), length.out = 80),
                     cen = temp_ref)
  pred_pctl <- crosspred(cb, model,
                          at = c(temp_p90, temp_p95, temp_p99),
                          cen = temp_ref)

  # Lag-specific RRs at 95th percentile
  idx_p95 <- which.min(abs(as.numeric(names(pred$allRRfit)) - temp_p95))
  lag_rr <- data.frame(
    group = grp,
    label = group_labels[grp],
    lag   = 0:max_lag,
    rr    = pred$matRRfit[idx_p95, ],
    rr_lo = pred$matRRlow[idx_p95, ],
    rr_hi = pred$matRRhigh[idx_p95, ]
  )

  # Which lag has the peak effect?
  peak_lag <- which.max(abs(log(lag_rr$rr))) - 1  # 0-indexed

  result_row <- data.frame(
    group = grp,
    label = group_labels[grp],
    type  = ifelse(grp %in% c("resp_relievers", "mh_anxiolytics"), "demand", "supply"),
    rr_p90    = round(pred_pctl$allRRfit[1], 4),
    rr_p90_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p90_hi = round(pred_pctl$allRRhigh[1], 4),
    rr_p95    = round(pred_pctl$allRRfit[2], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[2], 4),
    rr_p99    = round(pred_pctl$allRRfit[3], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[3], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[3], 4),
    peak_lag  = peak_lag,
    peak_lag_rr = round(lag_rr$rr[peak_lag + 1], 4),
    n_obs     = nrow(df),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  cat("  RR@p95:", result_row$rr_p95, "  Peak lag:", peak_lag, "weeks\n")

  all_results[[grp]] <- result_row
  all_preds[[grp]]   <- list(pred = pred, lag_rr = lag_rr,
                              temp_ref = temp_ref)
}


# ==============================================================================
# 3. COMPUTE DEMAND/SUPPLY RR RATIOS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Demand/Supply RR ratios\n")
cat("=" |> strrep(70), "\n")

results_table <- bind_rows(all_results)

for (comp_name in names(comparison_groups)) {
  comp <- comparison_groups[[comp_name]]
  demand_rr <- results_table$rr_p95[results_table$group == comp$demand]
  supply_rr <- results_table$rr_p95[results_table$group == comp$supply]

  if (length(demand_rr) > 0 && length(supply_rr) > 0) {
    ratio <- demand_rr / supply_rr
    cat("  ", comp$label, "— Demand/Supply RR ratio:", round(ratio, 3), "\n")
    cat("    Demand (", group_labels[comp$demand], "): RR =", demand_rr, "\n")
    cat("    Supply (", group_labels[comp$supply], "): RR =", supply_rr, "\n")
  }
}


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

write.csv(results_table,
          file.path(tab_dir, "demand_supply_comparison.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/demand_supply_comparison.csv\n")

lag_profiles <- bind_rows(lapply(all_preds, function(x) x$lag_rr))
write.csv(lag_profiles,
          file.path(tab_dir, "demand_supply_lag_profiles.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/demand_supply_lag_profiles.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# --- A. Forest plot: demand vs supply ---
results_table$type_label <- ifelse(results_table$type == "demand",
                                    "Demand-sensitive (acute/PRN)",
                                    "Supply-driven (scheduled)")
results_table$type_label[results_table$group == "cvd_total"] <- "Reference (CVD)"

p_forest <- ggplot(results_table,
                    aes(x = rr_p95, y = reorder(label, rr_p95),
                        colour = type_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.3, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c("Demand-sensitive (acute/PRN)" = "#d73027",
               "Supply-driven (scheduled)" = "#4575b4",
               "Reference (CVD)" = "grey50"),
    name = "Medication type"
  ) +
  labs(
    title = "Demand-sensitive vs supply-driven medications",
    subtitle = paste0("Cumulative RR at 95th percentile temperature\n",
                      "If signal is demand-driven, acute/PRN medications ",
                      "should show stronger temperature response"),
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "demand_supply_forest.png"),
       p_forest, width = 10, height = 6, dpi = 300)
cat("  -> Saved: demand_supply_forest.png\n")

# --- B. Lag response comparison ---
lag_data <- lag_profiles |>
  filter(group %in% c("resp_relievers", "resp_preventers")) |>
  mutate(label = group_labels[group])

p_lag <- ggplot(lag_data, aes(x = lag, y = rr, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(0, max_lag, 2)) +
  scale_colour_manual(
    values = c("Resp. Relievers (SABA)" = "#d73027",
               "Resp. Preventers (ICS)" = "#4575b4"),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Resp. Relievers (SABA)" = "#d73027",
               "Resp. Preventers (ICS)" = "#4575b4"),
    name = NULL
  ) +
  labs(
    title = "Lag response: relievers vs preventers at 95th percentile",
    subtitle = paste0("Faster response in relievers supports demand-driven ",
                      "signal (acute symptom exacerbation)"),
    x = "Lag (weeks)",
    y = "RR"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

# Mental health comparison
lag_data_mh <- lag_profiles |>
  filter(group %in% c("mh_anxiolytics", "mh_antidepressants")) |>
  mutate(label = group_labels[group])

p_lag_mh <- ggplot(lag_data_mh, aes(x = lag, y = rr, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = seq(0, max_lag, 2)) +
  scale_colour_manual(
    values = c("Anxiolytics (PRN)" = "#d73027",
               "Antidepressants (daily)" = "#4575b4"),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Anxiolytics (PRN)" = "#d73027",
               "Antidepressants (daily)" = "#4575b4"),
    name = NULL
  ) +
  labs(
    title = "Lag response: anxiolytics vs antidepressants at 95th percentile",
    subtitle = "PRN anxiolytics should show faster response than daily antidepressants",
    x = "Lag (weeks)",
    y = "RR"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

combined <- p_lag / p_lag_mh
ggsave(file.path(fig_dir, "demand_supply_lag_comparison.png"),
       combined, width = 10, height = 10, dpi = 300)
cat("  -> Saved: demand_supply_lag_comparison.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 35 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Demand vs supply comparison (RR at 95th percentile):\n")
results_table |>
  select(label, type_label, rr_p95, rr_p95_lo, rr_p95_hi, peak_lag) |>
  print(row.names = FALSE)

cat("\nInterpretation:\n")
cat("  - If demand-sensitive meds show stronger/faster response than supply-driven:\n")
cat("    evidence that the signal reflects genuine symptom exacerbation\n")
cat("  - If both show similar patterns: cannot distinguish demand from supply\n")
cat("  - Relievers (SABA) are the cleanest proxy for acute respiratory distress\n")
cat("  - Antidepressants (daily) should be insensitive to acute temperature\n")
cat("    (negative control for supply vs demand)\n")

cat("\nOutputs:\n")
cat("  outputs/tables/demand_supply_comparison.csv\n")
cat("  outputs/tables/demand_supply_lag_profiles.csv\n")
cat("  outputs/figures/demand_supply_forest.png\n")
cat("  outputs/figures/demand_supply_lag_comparison.png\n")
