# ==============================================================================
# Script 14: Compound Heat-Drought Analysis
# ==============================================================================
# Tests whether combined heat + drought events have amplified effects on
# medication prescribing compared to heat alone.
#
# Approach:
#   1. Classify weeks as "drought" using low rainfall thresholds
#   2. Add heat x drought interaction term to the DLNM model
#   3. Compare exposure-response curves under dry vs wet conditions
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/compound_heat_drought_rr.csv
#   outputs/tables/compound_interaction_tests.csv
#   outputs/figures/compound_er_overlay_*.png
#   outputs/figures/compound_forest.png
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

cat("=" |> strrep(70), "\n")
cat("Script 14: Compound Heat-Drought Analysis\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA AND CLASSIFY DROUGHT
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
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )

# Classify drought: weeks with rainfall below 25th percentile for that SA4
# This is a relative threshold — dry for the local climate
cat("Classifying drought conditions...\n")

# Compute SA4-specific rainfall percentiles
rainfall_thresholds <- panel |>
  filter(!is.na(precip_total)) |>
  group_by(sa4_code) |>
  summarise(
    rain_p25 = quantile(precip_total, 0.25, na.rm = TRUE),
    rain_p50 = quantile(precip_total, 0.50, na.rm = TRUE),
    .groups = "drop"
  )

panel <- panel |>
  left_join(rainfall_thresholds, by = "sa4_code") |>
  mutate(
    drought = ifelse(!is.na(precip_total) & precip_total <= rain_p25, 1L, 0L),
    drought_f = factor(drought, levels = c(0, 1),
                       labels = c("Normal/wet", "Dry (<=p25)"))
  )

n_drought <- sum(panel$drought == 1 & !is.na(panel$drought))
n_total   <- sum(!is.na(panel$drought))
cat("  Drought weeks:", format(n_drought, big.mark = ","), "/",
    format(n_total, big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_drought / n_total))


# ==============================================================================
# 2. STRATIFIED DLNM: DRY vs NORMAL/WET
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting stratified DLNM models (3 groups x 2 conditions = 6 models)\n")
cat("=" |> strrep(70), "\n")

drought_labels  <- c("Normal/wet", "Dry (<=p25)")
drought_colours <- c("Normal/wet" = "#4575b4", "Dry (<=p25)" = "#d73027")

fit_dlnm_drought <- function(data, group_name, drought_val) {
  drought_label <- drought_labels[drought_val + 1]

  df <- data |>
    filter(analysis_group == group_name,
           drought == drought_val,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 300) {
    cat("    ", drought_label, ": too few obs (", nrow(df), "), skipping\n")
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
    group     = group_name,
    condition = drought_label,
    drought   = drought_val,
    temp_ref  = round(temp_ref, 1),
    temp_p95  = round(temp_p95, 1),
    rr_p95    = round(pred_pctl$allRRfit[1], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[1], 4),
    temp_p99  = round(temp_p99, 1),
    rr_p99    = round(pred_pctl$allRRfit[2], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[2], 4),
    n_obs     = nrow(df),
    n_sa4     = n_distinct(df$sa4_code),
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
  for (d in c(0, 1)) {
    cat("  ", drought_labels[d + 1], "... ", sep = "")
    res <- fit_dlnm_drought(panel, grp, d)
    if (!is.null(res)) {
      key <- paste0(grp, "_d", d)
      all_results[[key]] <- res
      all_summaries[[key]] <- res$summary
      cat("RR@p95 =", res$summary$rr_p95,
          "(", res$summary$rr_p95_lo, "-", res$summary$rr_p95_hi, ")\n")
    }
  }
}


# ==============================================================================
# 3. INTERACTION TEST (pooled model with drought x temperature)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Testing heat x drought interaction\n")
cat("=" |> strrep(70), "\n")

interaction_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(drought)) |>
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

  # Null model (drought as additive confounder)
  cat("  Fitting null model...\n")
  mod_null <- glm(
    count ~ cb + drought_f +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  # Interaction model
  cat("  Fitting interaction model (tmax_mean x drought)...\n")
  mod_int <- glm(
    count ~ cb + drought_f + tmax_mean:drought_f +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
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
          file.path(tab_dir, "compound_heat_drought_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/compound_heat_drought_rr.csv\n")

int_table <- bind_rows(interaction_results)
write.csv(int_table,
          file.path(tab_dir, "compound_interaction_tests.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/compound_interaction_tests.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating compound event plots...\n")

# --- Overlaid ER curves: dry vs wet ---
plot_compound_er <- function(grp) {
  lbl <- group_labels[grp]

  plot_data <- list()
  for (d in c(0, 1)) {
    key <- paste0(grp, "_d", d)
    if (!key %in% names(all_results)) next
    res <- all_results[[key]]
    pred <- res$pred

    df_q <- data.frame(
      temp      = as.numeric(names(pred$allRRfit)),
      rr        = pred$allRRfit,
      rr_lo     = pred$allRRlow,
      rr_hi     = pred$allRRhigh,
      condition = drought_labels[d + 1]
    )
    plot_data[[as.character(d)]] <- df_q
  }

  if (length(plot_data) == 0) return(NULL)
  df_plot <- bind_rows(plot_data)

  p <- ggplot(df_plot, aes(x = temp, y = rr, colour = condition, fill = condition)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = drought_colours, name = "Rainfall condition") +
    scale_fill_manual(values = drought_colours, name = "Rainfall condition") +
    labs(
      title = paste0(lbl, " — Heat effect by rainfall condition"),
      subtitle = "Drought defined as weekly rainfall <= 25th percentile (SA4-specific)",
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom")

  return(p)
}

for (grp in primary_groups) {
  p <- plot_compound_er(grp)
  if (!is.null(p)) {
    ggsave(file.path(fig_dir, paste0("compound_er_overlay_", grp, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
}
cat("  -> Saved: compound_er_overlay_*.png\n")

# Combined panel
plots <- lapply(primary_groups, plot_compound_er)
plots <- plots[!sapply(plots, is.null)]
if (length(plots) >= 2) {
  combined <- wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "compound_er_panel.png"),
         combined, width = 9, height = 14, dpi = 300)
  cat("  -> Saved: compound_er_panel.png\n")
}

# --- Forest plot ---
cat("Generating compound forest plot...\n")

forest_data <- strat_table |>
  mutate(
    condition_f = factor(condition,
                         levels = rev(c("Normal/wet", "Dry (<=p25)"))),
    label_f = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory"))
  )

p_forest <- ggplot(forest_data,
                    aes(x = rr_p95, y = condition_f, colour = condition)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = drought_colours, guide = "none") +
  facet_wrap(~ label_f, scales = "free_x", ncol = 3) +
  labs(
    title = "Cumulative RR at 95th percentile: dry vs normal/wet conditions",
    subtitle = "Testing compound heat-drought effects on prescribing",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "compound_forest.png"),
       p_forest, width = 12, height = 4, dpi = 300)
cat("  -> Saved: compound_forest.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 14 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Compound heat-drought RRs at 95th percentile:\n")
strat_table |>
  select(label, condition, rr_p95, rr_p95_lo, rr_p95_hi, n_obs) |>
  print(row.names = FALSE)

cat("\nInteraction tests (heat x drought):\n")
int_table |>
  select(label, f_stat, p_value, significant) |>
  print(row.names = FALSE)

cat("\nOutputs:\n")
cat("  outputs/tables/compound_heat_drought_rr.csv\n")
cat("  outputs/tables/compound_interaction_tests.csv\n")
cat("  outputs/figures/compound_er_overlay_*.png\n")
cat("  outputs/figures/compound_forest.png\n")
