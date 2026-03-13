# ==============================================================================
# Script 13: Temporal Adaptation Analysis
# ==============================================================================
# Tests whether the heat-prescribing relationship has changed over the 10-year
# study period by splitting into early (2013-2017) and late (2018-2023) periods
# and comparing the DLNM exposure-response curves.
#
# If the heat effect has attenuated, it suggests population-level adaptation;
# if strengthened, it supports growing climate vulnerability.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/temporal_adaptation_rr.csv
#   outputs/tables/temporal_interaction_tests.csv
#   outputs/figures/temporal_er_overlay_*.png
#   outputs/figures/temporal_forest.png
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

# --- Project paths ---
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

# Period definitions
period_split <- as.Date("2018-01-01")
period_labels <- c("early" = "2013-2017", "late" = "2018-2023")
period_colours <- c("early" = "#4575b4", "late" = "#d73027")

cat("=" |> strrep(70), "\n")
cat("Script 13: Temporal Adaptation Analysis\n")
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
    period = ifelse(week_start < period_split, "early", "late"),
    period = factor(period, levels = c("early", "late")),
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    week_of_year = isoweek(week_start),
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )

cat("  Early period (2013-2017):", sum(panel$period == "early" &
      panel$analysis_group == "resp_total", na.rm = TRUE), "rows\n")
cat("  Late period (2018-2023):", sum(panel$period == "late" &
      panel$analysis_group == "resp_total", na.rm = TRUE), "rows\n")


# ==============================================================================
# 2. STRATIFIED DLNM BY TIME PERIOD
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models by period (3 groups x 2 periods = 6 models)\n")
cat("=" |> strrep(70), "\n")

fit_dlnm_period <- function(data, group_name, prd) {
  df <- data |>
    filter(analysis_group == group_name,
           period == prd,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 300) {
    cat("    ", period_labels[prd], ": too few obs (", nrow(df), "), skipping\n")
    return(NULL)
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

  model <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  temp_ref   <- median(df$tmax_mean, na.rm = TRUE)
  temp_range <- range(df$tmax_mean, na.rm = TRUE)

  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 80),
    cen = temp_ref
  )

  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)
  pred_pctl <- crosspred(cb, model, at = c(temp_p95, temp_p99), cen = temp_ref)

  rr_row <- data.frame(
    group    = group_name,
    period   = prd,
    period_label = period_labels[prd],
    temp_ref = round(temp_ref, 1),
    temp_p95 = round(temp_p95, 1),
    rr_p95    = round(pred_pctl$allRRfit[1], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[1], 4),
    temp_p99 = round(temp_p99, 1),
    rr_p99    = round(pred_pctl$allRRfit[2], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[2], 4),
    n_obs      = nrow(df),
    n_sa4      = n_distinct(df$sa4_code),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  return(list(pred = pred, summary = rr_row, model = model,
              cb = cb, data = df, temp_ref = temp_ref))
}

all_results   <- list()
all_summaries <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  for (prd in c("early", "late")) {
    cat("  ", period_labels[prd], "... ", sep = "")
    res <- fit_dlnm_period(panel, grp, prd)
    if (!is.null(res)) {
      key <- paste0(grp, "_", prd)
      all_results[[key]] <- res
      all_summaries[[key]] <- res$summary
      cat("RR@p95 =", res$summary$rr_p95,
          "(", res$summary$rr_p95_lo, "-", res$summary$rr_p95_hi, ")\n")
    }
  }
}


# ==============================================================================
# 3. INTERACTION TEST (period x temperature)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Testing heat x period interaction\n")
cat("=" |> strrep(70), "\n")

interaction_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  n_years  <- n_distinct(year(df$week_start))
  trend_df <- max(2, round(n_years * 2))

  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  # Null model
  cat("  Fitting null model...\n")
  mod_null <- glm(
    count ~ cb + period +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  # Interaction model
  cat("  Fitting interaction model (tmax_mean x period)...\n")
  mod_int <- glm(
    count ~ cb + period + tmax_mean:period +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  dev_null <- deviance(mod_null)
  dev_int  <- deviance(mod_int)
  df_diff  <- mod_null$df.residual - mod_int$df.residual
  disp     <- summary(mod_int)$dispersion

  f_stat <- ((dev_null - dev_int) / df_diff) / disp
  p_val  <- pf(f_stat, df_diff, mod_int$df.residual, lower.tail = FALSE)

  cat("  F-statistic:", round(f_stat, 3), "  df:", df_diff, "\n")
  cat("  P-value:", format.pval(p_val, digits = 3), "\n")

  interaction_results[[grp]] <- data.frame(
    group       = grp,
    label       = group_labels[grp],
    f_stat      = round(f_stat, 3),
    df_diff     = df_diff,
    p_value     = p_val,
    significant = p_val < 0.05,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

strat_table <- bind_rows(all_summaries)
strat_table$label <- group_labels[strat_table$group]

write.csv(strat_table,
          file.path(tab_dir, "temporal_adaptation_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/temporal_adaptation_rr.csv\n")

int_table <- bind_rows(interaction_results)
write.csv(int_table,
          file.path(tab_dir, "temporal_interaction_tests.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/temporal_interaction_tests.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating temporal comparison plots...\n")

# --- Overlaid ER curves by period ---
plot_temporal_er <- function(grp) {
  lbl <- group_labels[grp]

  plot_data <- list()
  for (prd in c("early", "late")) {
    key <- paste0(grp, "_", prd)
    if (!key %in% names(all_results)) next
    res <- all_results[[key]]
    pred <- res$pred

    df_q <- data.frame(
      temp   = as.numeric(names(pred$allRRfit)),
      rr     = pred$allRRfit,
      rr_lo  = pred$allRRlow,
      rr_hi  = pred$allRRhigh,
      period = prd
    )
    plot_data[[prd]] <- df_q
  }

  if (length(plot_data) == 0) return(NULL)
  df_plot <- bind_rows(plot_data)

  p <- ggplot(df_plot, aes(x = temp, y = rr, colour = period, fill = period)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = period_colours, labels = period_labels,
                        name = "Period") +
    scale_fill_manual(values = period_colours, labels = period_labels,
                      name = "Period") +
    labs(
      title = paste0(lbl, " — Exposure-response by period"),
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom")

  return(p)
}

for (grp in primary_groups) {
  p <- plot_temporal_er(grp)
  if (!is.null(p)) {
    ggsave(file.path(fig_dir, paste0("temporal_er_overlay_", grp, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
}
cat("  -> Saved: temporal_er_overlay_*.png\n")

# Combined panel
plots <- lapply(primary_groups, plot_temporal_er)
plots <- plots[!sapply(plots, is.null)]
if (length(plots) >= 2) {
  combined <- wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "temporal_er_panel.png"),
         combined, width = 9, height = 14, dpi = 300)
  cat("  -> Saved: temporal_er_panel.png\n")
}

# --- Forest plot ---
cat("Generating temporal forest plot...\n")

forest_data <- strat_table |>
  mutate(
    period_f = factor(period_label, levels = rev(c("2013-2017", "2018-2023"))),
    label_f  = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory"))
  )

p_forest <- ggplot(forest_data,
                    aes(x = rr_p95, y = period_f, colour = period)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = period_colours, guide = "none") +
  facet_wrap(~ label_f, scales = "free_x", ncol = 3) +
  labs(
    title = "Cumulative RR at 95th percentile: early vs late period",
    subtitle = "Testing for temporal adaptation in heat-prescribing relationship",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "temporal_forest.png"),
       p_forest, width = 12, height = 4, dpi = 300)
cat("  -> Saved: temporal_forest.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 13 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Period-stratified RRs at 95th percentile:\n")
strat_table |>
  select(label, period_label, rr_p95, rr_p95_lo, rr_p95_hi, n_obs) |>
  print(row.names = FALSE)

cat("\nInteraction tests (heat x period):\n")
int_table |>
  select(label, f_stat, p_value, significant) |>
  print(row.names = FALSE)

cat("\nInterpretation guide:\n")
cat("  - If late RR closer to 1 than early: suggests adaptation\n")
cat("  - If late RR further from 1 than early: suggests increasing vulnerability\n")
cat("  - Non-significant interaction: no evidence of temporal change\n")

cat("\nOutputs:\n")
cat("  outputs/tables/temporal_adaptation_rr.csv\n")
cat("  outputs/tables/temporal_interaction_tests.csv\n")
cat("  outputs/figures/temporal_er_overlay_*.png\n")
cat("  outputs/figures/temporal_forest.png\n")
