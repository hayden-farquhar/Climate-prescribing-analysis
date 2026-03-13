# ==============================================================================
# Script 21: Negative Control Outcome
# ==============================================================================
# Tests the DLNM on a medication class with no plausible heat mechanism to
# rule out systematic bias (e.g., PBS reporting artefacts during hot weeks).
#
# Antidepressants (mh_antidepressants) serve as the negative control:
#   - Chronic medications taken daily regardless of weather
#   - No acute physiological pathway linking heat to antidepressant dispensing
#   - Already showed null effect in script 05 (RR 0.999 at p95)
#
# A null finding for the negative control, combined with a significant
# finding for respiratory, supports causal inference (rules out that heat
# weeks simply have fewer pharmacy visits).
#
# Also tests anxiolytics as a second negative control.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/negative_control_results.csv
#   outputs/figures/negative_control_comparison.png
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

# Test outcome: negative controls + positive control (respiratory)
test_groups <- c("mh_antidepressants", "mh_anxiolytics", "resp_total")
group_labels <- c(
  mh_antidepressants = "Antidepressants (negative control)",
  mh_anxiolytics     = "Anxiolytics (negative control)",
  resp_total         = "Respiratory (positive control)"
)

group_roles <- c(
  mh_antidepressants = "Negative control",
  mh_anxiolytics     = "Negative control",
  resp_total         = "Positive control"
)

cat("=" |> strrep(70), "\n")
cat("Script 21: Negative Control Outcome Analysis\n")
cat("=" |> strrep(70), "\n\n")

cat("Rationale:\n")
cat("  If heat weeks simply had fewer pharmacy visits (bias), ALL medication\n")
cat("  classes would show reduced dispensing. A null finding for chronic\n")
cat("  medications (antidepressants) alongside a significant finding for\n")
cat("  respiratory medications supports a specific biological mechanism.\n\n")


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

results <- list()

for (grp in test_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
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

  model <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_range <- range(df$tmax_mean, na.rm = TRUE)

  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 80),
    cen = temp_ref
  )

  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)
  pred_pctl <- crosspred(cb, model, at = c(temp_p90, temp_p95, temp_p99),
                          cen = temp_ref)

  rr_row <- data.frame(
    group    = grp,
    label    = group_labels[grp],
    role     = group_roles[grp],
    temp_ref = round(temp_ref, 1),
    rr_p90   = round(pred_pctl$allRRfit[1], 4),
    rr_p90_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p90_hi = round(pred_pctl$allRRhigh[1], 4),
    rr_p95   = round(pred_pctl$allRRfit[2], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[2], 4),
    rr_p99   = round(pred_pctl$allRRfit[3], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[3], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[3], 4),
    n_obs    = nrow(df),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  results[[grp]] <- list(pred = pred, summary = rr_row,
                          temp_ref = temp_ref)

  cat("  Role:", group_roles[grp], "\n")
  cat("  RR at p95:", rr_row$rr_p95,
      "(", rr_row$rr_p95_lo, "-", rr_row$rr_p95_hi, ")\n")
  if (group_roles[grp] == "Negative control") {
    if (rr_row$rr_p95_lo <= 1 & rr_row$rr_p95_hi >= 1) {
      cat("  PASS: 95% CI includes 1 (null) — no heat effect detected\n")
    } else {
      cat("  NOTE: 95% CI excludes 1 — unexpected for negative control\n")
    }
  }
}


# ==============================================================================
# 2. SAVE RESULTS
# ==============================================================================

nc_table <- bind_rows(lapply(results, function(r) r$summary))
write.csv(nc_table,
          file.path(tab_dir, "negative_control_results.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/negative_control_results.csv\n")


# ==============================================================================
# 3. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# Side-by-side ER curves
plot_data_all <- list()
for (grp in test_groups) {
  pred <- results[[grp]]$pred
  plot_data_all[[grp]] <- data.frame(
    temp  = as.numeric(names(pred$allRRfit)),
    rr    = pred$allRRfit,
    rr_lo = pred$allRRlow,
    rr_hi = pred$allRRhigh,
    group = group_labels[grp],
    role  = group_roles[grp]
  )
}
df_plot <- bind_rows(plot_data_all)

role_colours <- c("Negative control" = "#4575b4", "Positive control" = "#d73027")

p <- ggplot(df_plot, aes(x = temp, y = rr, colour = role, fill = role)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ group, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = role_colours, name = "Role") +
  scale_fill_manual(values = role_colours, name = "Role") +
  labs(
    title = "Negative control analysis: heat-prescribing exposure-response",
    subtitle = "Negative controls should show flat response; positive control should show effect",
    x = "Weekly mean Tmax (°C)",
    y = "Cumulative RR (vs median)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "negative_control_comparison.png"),
       p, width = 15, height = 5, dpi = 300)
cat("  -> Saved: negative_control_comparison.png\n")

# Forest plot
p_forest <- ggplot(nc_table,
                    aes(x = rr_p95, y = reorder(label, rr_p95),
                        colour = role)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.3, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_colour_manual(values = role_colours, name = "Role") +
  labs(
    title = "Negative control test: RR at 95th percentile",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "negative_control_forest.png"),
       p_forest, width = 9, height = 4, dpi = 300)
cat("  -> Saved: negative_control_forest.png\n")


# ==============================================================================
# 4. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 21 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Negative control results (RR at p95):\n")
print(as.data.frame(
  nc_table |> select(label, role, rr_p95, rr_p95_lo, rr_p95_hi)
), row.names = FALSE)

cat("\nConclusion:\n")
neg <- nc_table |> filter(role == "Negative control")
pos <- nc_table |> filter(role == "Positive control")
all_neg_null <- all(neg$rr_p95_lo <= 1 & neg$rr_p95_hi >= 1)
pos_sig <- any(pos$rr_p95_hi < 1 | pos$rr_p95_lo > 1)

if (all_neg_null & pos_sig) {
  cat("  SUPPORTS CAUSAL INFERENCE: Negative controls show null; positive\n")
  cat("  control shows significant effect. Bias from reduced pharmacy visits\n")
  cat("  during heat is unlikely to explain findings.\n")
} else {
  cat("  Results require careful interpretation — see output table.\n")
}

cat("\nOutputs:\n")
cat("  outputs/tables/negative_control_results.csv\n")
cat("  outputs/figures/negative_control_comparison.png\n")
cat("  outputs/figures/negative_control_forest.png\n")
