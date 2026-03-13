# ==============================================================================
# Script 11: Cold Effects Analysis
# ==============================================================================
# Extracts and reports the full temperature-prescribing relationship from
# existing DLNM models, with specific focus on cold extremes (5th, 10th
# percentiles). The original script 05 only reported heat effects; this
# script provides the complete exposure-response characterisation.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/cold_heat_rr_summary.csv
#   outputs/figures/full_er_curve_*.png
#   outputs/figures/full_er_panel.png
#   outputs/figures/cold_heat_forest.png
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

outcome_groups <- c(
  "cvd_total", "mh_total", "resp_total",
  "mh_antidepressants", "mh_anxiolytics",
  "resp_relievers", "resp_preventers"
)

group_labels <- c(
  cvd_total          = "Cardiovascular",
  mh_total           = "Mental Health",
  resp_total         = "Respiratory",
  mh_antidepressants = "Antidepressants",
  mh_anxiolytics     = "Anxiolytics",
  resp_relievers     = "Resp. Relievers",
  resp_preventers    = "Resp. Preventers"
)

cat("=" |> strrep(70), "\n")
cat("Script 11: Cold and Heat Effects — Full Exposure-Response\n")
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
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )


# ==============================================================================
# 2. FIT DLNM AND EXTRACT COLD + HEAT PERCENTILES
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DLNM models with cold + heat percentile extraction\n")
cat("=" |> strrep(70), "\n")

all_results  <- list()
all_summaries <- list()

for (grp in outcome_groups) {
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

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

  # Cross-basis
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

  # Full exposure-response prediction
  pred <- crosspred(
    cb, model,
    at = seq(temp_range[1], temp_range[2], length.out = 100),
    cen = temp_ref
  )

  # Extract at cold AND heat percentiles
  pctls <- c(0.01, 0.05, 0.10, 0.90, 0.95, 0.99)
  temp_pctls <- quantile(df$tmax_mean, probs = pctls, na.rm = TRUE)

  pred_pctl <- crosspred(cb, model, at = temp_pctls, cen = temp_ref)

  rr_row <- data.frame(
    group    = grp,
    label    = group_labels[grp],
    temp_ref = round(temp_ref, 1),
    # Cold percentiles
    temp_p01 = round(temp_pctls[1], 1),
    rr_p01    = round(pred_pctl$allRRfit[1], 4),
    rr_p01_lo = round(pred_pctl$allRRlow[1], 4),
    rr_p01_hi = round(pred_pctl$allRRhigh[1], 4),
    temp_p05 = round(temp_pctls[2], 1),
    rr_p05    = round(pred_pctl$allRRfit[2], 4),
    rr_p05_lo = round(pred_pctl$allRRlow[2], 4),
    rr_p05_hi = round(pred_pctl$allRRhigh[2], 4),
    temp_p10 = round(temp_pctls[3], 1),
    rr_p10    = round(pred_pctl$allRRfit[3], 4),
    rr_p10_lo = round(pred_pctl$allRRlow[3], 4),
    rr_p10_hi = round(pred_pctl$allRRhigh[3], 4),
    # Heat percentiles
    temp_p90 = round(temp_pctls[4], 1),
    rr_p90    = round(pred_pctl$allRRfit[4], 4),
    rr_p90_lo = round(pred_pctl$allRRlow[4], 4),
    rr_p90_hi = round(pred_pctl$allRRhigh[4], 4),
    temp_p95 = round(temp_pctls[5], 1),
    rr_p95    = round(pred_pctl$allRRfit[5], 4),
    rr_p95_lo = round(pred_pctl$allRRlow[5], 4),
    rr_p95_hi = round(pred_pctl$allRRhigh[5], 4),
    temp_p99 = round(temp_pctls[6], 1),
    rr_p99    = round(pred_pctl$allRRfit[6], 4),
    rr_p99_lo = round(pred_pctl$allRRlow[6], 4),
    rr_p99_hi = round(pred_pctl$allRRhigh[6], 4),
    n_obs = nrow(df),
    stringsAsFactors = FALSE
  )

  all_results[[grp]]  <- list(pred = pred, summary = rr_row,
                               temp_ref = temp_ref, data = df)
  all_summaries[[grp]] <- rr_row

  cat("  Cold p05: RR =", rr_row$rr_p05,
      "(", rr_row$rr_p05_lo, "-", rr_row$rr_p05_hi, ")\n")
  cat("  Heat p95: RR =", rr_row$rr_p95,
      "(", rr_row$rr_p95_lo, "-", rr_row$rr_p95_hi, ")\n")
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

summary_table <- bind_rows(all_summaries)
write.csv(summary_table,
          file.path(tab_dir, "cold_heat_rr_summary.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/cold_heat_rr_summary.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating full exposure-response plots...\n")

# --- Full ER curves with cold and heat zones highlighted ---
plot_full_er <- function(grp) {
  res  <- all_results[[grp]]
  pred <- res$pred
  lbl  <- group_labels[grp]

  df_plot <- data.frame(
    temp  = as.numeric(names(pred$allRRfit)),
    rr    = pred$allRRfit,
    rr_lo = pred$allRRlow,
    rr_hi = pred$allRRhigh
  )

  # Temperature percentiles for shading
  df_data <- res$data
  temp_p05 <- quantile(df_data$tmax_mean, 0.05, na.rm = TRUE)
  temp_p10 <- quantile(df_data$tmax_mean, 0.10, na.rm = TRUE)
  temp_p90 <- quantile(df_data$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df_data$tmax_mean, 0.95, na.rm = TRUE)

  p <- ggplot(df_plot, aes(x = temp, y = rr)) +
    # Cold zone shading
    annotate("rect", xmin = -Inf, xmax = temp_p10,
             ymin = -Inf, ymax = Inf,
             fill = "#4575b4", alpha = 0.06) +
    # Heat zone shading
    annotate("rect", xmin = temp_p90, xmax = Inf,
             ymin = -Inf, ymax = Inf,
             fill = "#d73027", alpha = 0.06) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi),
                fill = "steelblue", alpha = 0.2) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = res$temp_ref, linetype = "dotted", colour = "grey60") +
    # Percentile labels
    annotate("text", x = temp_p05, y = Inf, label = "p5",
             vjust = 1.5, hjust = 0.5, size = 2.5, colour = "#4575b4") +
    annotate("text", x = temp_p95, y = Inf, label = "p95",
             vjust = 1.5, hjust = 0.5, size = 2.5, colour = "#d73027") +
    labs(
      title = lbl,
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  return(p)
}

# Individual plots
for (grp in names(all_results)) {
  p <- plot_full_er(grp)
  ggsave(file.path(fig_dir, paste0("full_er_curve_", grp, ".png")),
         p, width = 7, height = 5, dpi = 300)
}
cat("  -> Saved: full_er_curve_*.png\n")

# Panel of primary groups
primary <- intersect(c("cvd_total", "mh_total", "resp_total"), names(all_results))
plots <- lapply(primary, plot_full_er)
combined <- wrap_plots(plots, ncol = 3)
ggsave(file.path(fig_dir, "full_er_panel.png"),
       combined, width = 15, height = 5, dpi = 300)
cat("  -> Saved: full_er_panel.png\n")


# --- Forest plot: cold and heat side-by-side ---
cat("Generating cold/heat comparison forest plot...\n")

forest_data <- summary_table |>
  select(label, rr_p05, rr_p05_lo, rr_p05_hi, rr_p95, rr_p95_lo, rr_p95_hi) |>
  tidyr::pivot_longer(
    cols = starts_with("rr_p"),
    names_to = c(".value", "percentile"),
    names_pattern = "rr_(p\\d+)_(.*)"
  )

# Rebuild properly
forest_long <- bind_rows(
  summary_table |>
    transmute(label, percentile = "5th (cold)",
              rr = rr_p05, rr_lo = rr_p05_lo, rr_hi = rr_p05_hi),
  summary_table |>
    transmute(label, percentile = "95th (heat)",
              rr = rr_p95, rr_lo = rr_p95_lo, rr_hi = rr_p95_hi)
)

forest_long$percentile <- factor(forest_long$percentile,
                                  levels = c("5th (cold)", "95th (heat)"))

p_forest <- ggplot(forest_long,
                    aes(x = rr, y = reorder(label, rr),
                        colour = percentile)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi),
                 height = 0.3, linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    values = c("5th (cold)" = "#4575b4", "95th (heat)" = "#d73027"),
    name = "Temperature Percentile"
  ) +
  labs(
    title = "Cumulative RR at cold (5th) and heat (95th) percentiles",
    subtitle = "vs median temperature (reference)",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "cold_heat_forest.png"),
       p_forest, width = 9, height = 6, dpi = 300)
cat("  -> Saved: cold_heat_forest.png\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 11 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Cold and heat RRs:\n")
summary_table |>
  select(label, rr_p05, rr_p05_lo, rr_p05_hi,
         rr_p95, rr_p95_lo, rr_p95_hi) |>
  print(row.names = FALSE)

cat("\nOutputs:\n")
cat("  outputs/tables/cold_heat_rr_summary.csv\n")
cat("  outputs/figures/full_er_curve_*.png\n")
cat("  outputs/figures/full_er_panel.png\n")
cat("  outputs/figures/cold_heat_forest.png\n")
