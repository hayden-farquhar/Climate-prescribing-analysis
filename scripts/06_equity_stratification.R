# ==============================================================================
# Script 06: SEIFA Equity Stratification
# ==============================================================================
# Stratifies the DLNM heat-prescribing models by SEIFA IRSD quintile to test
# whether socioeconomically disadvantaged areas show greater sensitivity to
# heat exposure.
#
# Approach:
#   1. Merge SA4-level SEIFA quintiles into the panel
#   2. Fit DLNM separately for each quintile (3 primary medication groups)
#   3. Fit a pooled model with heat × SEIFA interaction to test for
#      differential effects
#   4. Generate comparative exposure-response curves by quintile
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv  (from Python spatial overlay)
#
# Outputs:
#   outputs/tables/equity_stratified_rr.csv
#   outputs/tables/equity_interaction_tests.csv
#   outputs/figures/equity_exposure_response_*.png
#   outputs/figures/equity_lag_response_*.png
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

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

# Medication groups and labels
primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

quintile_labels <- c(
  "1" = "Q1 (Most disadvantaged)",
  "2" = "Q2",
  "3" = "Q3",
  "4" = "Q4",
  "5" = "Q5 (Least disadvantaged)"
)

quintile_colours <- c(
  "1" = "#d73027",  # red
  "2" = "#fc8d59",
  "3" = "#fee090",
  "4" = "#91bfdb",
  "5" = "#4575b4"   # blue
)


# ==============================================================================
# 1. LOAD AND MERGE DATA
# ==============================================================================

cat("=" |> strrep(70), "\n")
cat("Script 06: SEIFA Equity Stratification\n")
cat("=" |> strrep(70), "\n\n")

# Load panel
cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

# Load SEIFA quintiles
cat("Loading SEIFA SA4 quintiles...\n")
seifa <- read.csv(file.path(ref_dir, "seifa_sa4_quintiles.csv"),
                  stringsAsFactors = FALSE)
seifa$sa4_code <- as.character(seifa$sa4_code)
cat("  SA4s with SEIFA:", nrow(seifa), "\n")
cat("  Quintile distribution:\n")
print(table(seifa$irsd_quintile))

# Merge
panel <- panel |>
  left_join(seifa |> select(sa4_code, irsd_score, irsd_quintile, population),
            by = "sa4_code")

n_matched <- sum(!is.na(panel$irsd_quintile))
cat("\n  Panel rows with SEIFA match:", format(n_matched, big.mark = ","),
    "/", format(nrow(panel), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_matched / nrow(panel)))

# Add time variables
panel <- panel |>
  mutate(
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    week_of_year = isoweek(week_start),
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )


# ==============================================================================
# 2. STRATIFIED DLNM BY SEIFA QUINTILE
# ==============================================================================

fit_dlnm_stratified <- function(data, group_name, quintile) {
  #' Fit a DLNM for a single medication group within one SEIFA quintile.

  df <- data |>
    filter(analysis_group == group_name,
           irsd_quintile == quintile,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 200) {
    cat("    Q", quintile, ": too few obs (", nrow(df), "), skipping\n")
    return(NULL)
  }

  # Cross-basis (same spec as script 05)
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years <- n_distinct(year(df$week_start))
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
    quintile = quintile,
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

  return(list(
    pred    = pred,
    summary = rr_row,
    model   = model,
    cb      = cb,
    data    = df,
    temp_ref = temp_ref
  ))
}


cat("\n", "=" |> strrep(70), "\n")
cat("Fitting stratified DLNM models (3 groups x 5 quintiles = 15 models)\n")
cat("=" |> strrep(70), "\n")

all_results <- list()
all_summaries <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  for (q in 1:5) {
    cat("  Quintile ", q, "... ", sep = "")
    res <- fit_dlnm_stratified(panel, grp, q)
    if (!is.null(res)) {
      key <- paste0(grp, "_q", q)
      all_results[[key]] <- res
      all_summaries[[key]] <- res$summary
      cat("RR@p95 =", res$summary$rr_p95,
          "(", res$summary$rr_p95_lo, "-", res$summary$rr_p95_hi, ")",
          "n =", format(res$summary$n_obs, big.mark = ","), "\n")
    }
  }
}


# ==============================================================================
# 3. INTERACTION TEST (pooled model with quintile × temperature)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Testing heat x SEIFA interaction\n")
cat("=" |> strrep(70), "\n")

interaction_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(irsd_quintile)) |>
    mutate(seifa_q = factor(irsd_quintile)) |>
    arrange(sa4_code, week_start)

  n_years <- n_distinct(year(df$week_start))
  trend_df <- max(2, round(n_years * 2))

  # Cross-basis
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  # Model WITHOUT interaction
  cat("  Fitting null model (no interaction)...\n")
  mod_null <- glm(
    count ~ cb + seifa_q +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  # Model WITH interaction (cross-basis × SEIFA quintile)
  # Use a linear temperature term for interaction (full cb interaction has too many params)
  cat("  Fitting interaction model (tmax_mean × SEIFA)...\n")
  mod_int <- glm(
    count ~ cb + seifa_q + tmax_mean:seifa_q +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  # F-test for interaction terms
  # Compare deviance of null vs interaction model
  dev_null <- deviance(mod_null)
  dev_int  <- deviance(mod_int)
  df_diff  <- mod_null$df.residual - mod_int$df.residual
  disp     <- summary(mod_int)$dispersion

  f_stat <- ((dev_null - dev_int) / df_diff) / disp
  p_val  <- pf(f_stat, df_diff, mod_int$df.residual, lower.tail = FALSE)

  cat("  Deviance reduction:", round(dev_null - dev_int, 1), "\n")
  cat("  F-statistic:", round(f_stat, 3), "  df:", df_diff, "\n")
  cat("  P-value:", format.pval(p_val, digits = 3), "\n")

  interaction_results[[grp]] <- data.frame(
    group       = grp,
    label       = group_labels[grp],
    dev_null    = round(dev_null, 1),
    dev_int     = round(dev_int, 1),
    dev_diff    = round(dev_null - dev_int, 1),
    df_diff     = df_diff,
    f_stat      = round(f_stat, 3),
    p_value     = p_val,
    significant = p_val < 0.05,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 4. SAVE RESULTS TABLES
# ==============================================================================

# Stratified RRs
strat_table <- bind_rows(all_summaries)
strat_table$label <- group_labels[strat_table$group]
strat_table$quintile_label <- quintile_labels[as.character(strat_table$quintile)]

write.csv(strat_table,
          file.path(tab_dir, "equity_stratified_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/equity_stratified_rr.csv\n")

# Interaction tests
int_table <- bind_rows(interaction_results)
write.csv(int_table,
          file.path(tab_dir, "equity_interaction_tests.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/equity_interaction_tests.csv\n")


# ==============================================================================
# 5. VISUALISATION: Exposure-Response by SEIFA Quintile
# ==============================================================================

cat("\nGenerating equity plots...\n")

plot_equity_exposure_response <- function(grp) {
  lbl <- group_labels[grp]

  # Collect predictions from each quintile's model
  plot_data <- list()
  for (q in 1:5) {
    key <- paste0(grp, "_q", q)
    if (!key %in% names(all_results)) next
    res <- all_results[[key]]
    pred <- res$pred

    df_q <- data.frame(
      temp     = as.numeric(names(pred$allRRfit)),
      rr       = pred$allRRfit,
      rr_lo    = pred$allRRlow,
      rr_hi    = pred$allRRhigh,
      quintile = as.character(q)
    )
    plot_data[[q]] <- df_q
  }

  if (length(plot_data) == 0) return(NULL)
  df_plot <- bind_rows(plot_data)

  p <- ggplot(df_plot, aes(x = temp, y = rr, colour = quintile, fill = quintile)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(
      values = quintile_colours,
      labels = quintile_labels,
      name = "SEIFA Quintile"
    ) +
    scale_fill_manual(
      values = quintile_colours,
      labels = quintile_labels,
      name = "SEIFA Quintile"
    ) +
    labs(
      title = paste0(lbl, " — Exposure-response by SEIFA quintile"),
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

  return(p)
}

for (grp in primary_groups) {
  p <- plot_equity_exposure_response(grp)
  if (!is.null(p)) {
    ggsave(file.path(fig_dir, paste0("equity_exposure_response_", grp, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
}
cat("  -> Saved: equity_exposure_response_*.png\n")

# Combined panel
plots <- lapply(primary_groups, plot_equity_exposure_response)
plots <- plots[!sapply(plots, is.null)]
if (length(plots) >= 2) {
  combined <- wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "equity_exposure_response_panel.png"),
         combined, width = 9, height = 14, dpi = 300)
  cat("  -> Saved: equity_exposure_response_panel.png\n")
}


# ==============================================================================
# 6. FOREST PLOT: RR at p95 by Quintile
# ==============================================================================

cat("Generating forest plots...\n")

forest_data <- strat_table |>
  mutate(
    quintile_f = factor(quintile, levels = 5:1),
    label_f    = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory"))
  )

p_forest <- ggplot(forest_data, aes(x = rr_p95, y = quintile_f, colour = factor(quintile))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi), height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = quintile_colours, guide = "none") +
  scale_y_discrete(labels = function(x) quintile_labels[x]) +
  facet_wrap(~ label_f, scales = "free_x", ncol = 3) +
  labs(
    title = "Cumulative RR at 95th percentile temperature by SEIFA quintile",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "equity_forest_plot.png"),
       p_forest, width = 12, height = 5, dpi = 300)
cat("  -> Saved: equity_forest_plot.png\n")


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 06 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Stratified RRs at 95th percentile:\n")
strat_table |>
  select(label, quintile_label, rr_p95, rr_p95_lo, rr_p95_hi, n_obs) |>
  print(n = 20, row.names = FALSE)

cat("\nInteraction tests (heat × SEIFA):\n")
int_table |>
  select(label, f_stat, p_value, significant) |>
  print(row.names = FALSE)

cat("\nOutputs:\n")
cat("  outputs/tables/equity_stratified_rr.csv\n")
cat("  outputs/tables/equity_interaction_tests.csv\n")
cat("  outputs/figures/equity_exposure_response_*.png\n")
cat("  outputs/figures/equity_forest_plot.png\n")

cat("\nNext step: Script 07 — Black Summer bushfire DiD analysis\n")
