# ==============================================================================
# Script 08: Publication-Quality Figures
# ==============================================================================
# Generates polished, journal-ready figures for the primary manuscript.
# Uses consistent styling, colour palettes, and layout for EHP/Lancet format.
#
# Figures:
#   Figure 1: Main exposure-response curves (3 primary groups, panel)
#   Figure 2: Equity forest plot (SEIFA-stratified RRs + interaction)
#   Figure 3: Bushfire DiD forest + event study (combined)
#   Figure 4: Sensitivity analysis summary (forest plot of all specifications)
#   Figure S1: Lag-response curves at p95
#   Figure S2: Remoteness + compound stratification
#   Figure S3: Negative control comparison
#   Figure S4: Temporal adaptation comparison
#
# Inputs:
#   outputs/tables/*.csv (from scripts 05-27)
#   data/processed/panel_weekly_sa4.csv (for re-fitting main models)
#
# Outputs:
#   outputs/figures/fig1_exposure_response.png / .pdf
#   outputs/figures/fig2_equity_forest.png / .pdf
#   outputs/figures/fig3_bushfire.png / .pdf
#   outputs/figures/fig4_sensitivity_forest.png / .pdf
#   outputs/figures/figS1_lag_response.png / .pdf
#   outputs/figures/figS2_stratification.png / .pdf
#   outputs/figures/figS3_negative_control.png / .pdf
#   outputs/figures/figS4_temporal.png / .pdf
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate",
#                       "ggplot2", "patchwork", "scales"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(scales)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

# --- Consistent theme for all figures ---
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2,
                                       hjust = 0, margin = margin(b = 8)),
      plot.subtitle    = element_text(size = base_size, colour = "grey30",
                                       margin = margin(b = 10)),
      strip.text       = element_text(face = "bold", size = base_size),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold", size = base_size),
      legend.text      = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# Colour palette
pal_primary   <- c("Cardiovascular" = "#4575b4",
                    "Mental Health"  = "#fdae61",
                    "Respiratory"    = "#d73027")
pal_seifa     <- c("Q1 (most disadvantaged)" = "#d73027",
                    "Q2" = "#fc8d59", "Q3" = "#fee08b",
                    "Q4" = "#91bfdb", "Q5 (least disadvantaged)" = "#4575b4")

cat("=" |> strrep(70), "\n")
cat("Script 08: Publication-Quality Figures\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 0. LOAD AND FIT PRIMARY MODELS (for ER and lag curves)
# ==============================================================================

cat("Loading data and fitting primary models...\n")
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

primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

# Fit models and extract predictions
er_data <- list()
lag_data <- list()

for (grp in primary_groups) {
  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count), !is.na(tmax_mean), !is.na(precip_total)) |>
    arrange(sa4_code, week_start)

  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean, lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))

  model <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total + factor(sa4_code),
    data = df, family = quasipoisson(link = "log")
  )

  temp_ref   <- median(df$tmax_mean, na.rm = TRUE)
  temp_range <- range(df$tmax_mean, na.rm = TRUE)
  temp_p95   <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_seq   <- seq(temp_range[1], temp_range[2], length.out = 100)

  pred <- crosspred(cb, model, at = temp_seq, cen = temp_ref)

  er_data[[grp]] <- data.frame(
    temp  = as.numeric(names(pred$allRRfit)),
    rr    = pred$allRRfit,
    rr_lo = pred$allRRlow,
    rr_hi = pred$allRRhigh,
    group = group_labels[grp]
  )

  # Lag-response at p95
  pred_lag <- crosspred(cb, model, at = temp_p95, cen = temp_ref)
  lag_data[[grp]] <- data.frame(
    lag   = 0:max_lag,
    rr    = pred_lag$matRRfit[1, ],
    rr_lo = pred_lag$matRRlow[1, ],
    rr_hi = pred_lag$matRRhigh[1, ],
    group = group_labels[grp]
  )

  cat("  ", group_labels[grp], "fitted.\n")
}


# ==============================================================================
# FIGURE 1: Main Exposure-Response Curves
# ==============================================================================

cat("\nFigure 1: Exposure-response curves...\n")

er_df <- bind_rows(er_data)
er_df$group <- factor(er_df$group, levels = c("Cardiovascular", "Mental Health", "Respiratory"))

fig1 <- ggplot(er_df, aes(x = temp, y = rr)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi, fill = group), alpha = 0.15) +
  geom_line(aes(colour = group), linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_wrap(~ group, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = pal_primary) +
  scale_fill_manual(values = pal_primary) +
  labs(
    x = "Weekly mean maximum temperature (\u00b0C)",
    y = "Cumulative relative risk"
  ) +
  theme_pub() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "cm"))

ggsave(file.path(fig_dir, "fig1_exposure_response.png"),
       fig1, width = 13, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "fig1_exposure_response.pdf"),
       fig1, width = 13, height = 4.5)
cat("  -> Saved: fig1_exposure_response.png/pdf\n")


# ==============================================================================
# FIGURE 2: Equity Forest Plot
# ==============================================================================

cat("Figure 2: Equity forest plot...\n")

equity_file <- file.path(tab_dir, "equity_stratified_rr.csv")
if (file.exists(equity_file)) {
  eq <- read.csv(equity_file, stringsAsFactors = FALSE)

  # Filter to primary groups and p95
  eq_plot <- eq |>
    filter(group %in% primary_groups) |>
    mutate(
      label = group_labels[group],
      q_label = paste0("Q", quintile),
      label = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory"))
    )

  fig2 <- ggplot(eq_plot, aes(x = rr_p95, y = q_label, colour = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5) +
    geom_point(size = 2.5) +
    facet_wrap(~ label, ncol = 3, scales = "free_x") +
    scale_colour_manual(values = pal_primary) +
    scale_y_discrete(limits = rev(paste0("Q", 1:5)),
                     labels = rev(c("Q1\n(most\ndisadvantaged)", "Q2", "Q3", "Q4",
                                     "Q5\n(least\ndisadvantaged)"))) +
    labs(
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = "SEIFA IRSD quintile"
    ) +
    theme_pub() +
    theme(legend.position = "none",
          panel.spacing = unit(1.2, "cm"))

  ggsave(file.path(fig_dir, "fig2_equity_forest.png"),
         fig2, width = 13, height = 4.5, dpi = 300)
  ggsave(file.path(fig_dir, "fig2_equity_forest.pdf"),
         fig2, width = 13, height = 4.5)
  cat("  -> Saved: fig2_equity_forest.png/pdf\n")
} else {
  cat("  equity_stratified_rr.csv not found, skipping.\n")
}


# ==============================================================================
# FIGURE 3: Bushfire DiD
# ==============================================================================

cat("Figure 3: Bushfire DiD...\n")

bush_file <- file.path(tab_dir, "bushfire_did_results.csv")
if (file.exists(bush_file)) {
  bush <- read.csv(bush_file, stringsAsFactors = FALSE)

  # Medium and high smoke zones
  bush_plot <- bush |>
    tidyr::pivot_longer(
      cols = c(rr_medium, rr_high),
      names_to = "smoke_level",
      values_to = "rr"
    ) |>
    tidyr::pivot_longer(
      cols = c(rr_medium_lo, rr_high_lo),
      names_to = "lo_name",
      values_to = "rr_lo"
    ) |>
    tidyr::pivot_longer(
      cols = c(rr_medium_hi, rr_high_hi),
      names_to = "hi_name",
      values_to = "rr_hi"
    )

  # Simpler approach: reshape manually
  bush_medium <- bush |>
    transmute(label, smoke = "Medium smoke",
              rr = rr_medium, rr_lo = rr_medium_lo, rr_hi = rr_medium_hi) |>
    filter(!is.na(rr))
  bush_high <- bush |>
    transmute(label, smoke = "High smoke",
              rr = rr_high, rr_lo = rr_high_lo, rr_hi = rr_high_hi) |>
    filter(!is.na(rr))
  bush_long <- bind_rows(bush_medium, bush_high) |>
    mutate(smoke = factor(smoke, levels = c("Medium smoke", "High smoke")))

  fig3 <- ggplot(bush_long,
                  aes(x = rr, y = reorder(label, rr), colour = smoke)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi),
                   height = 0.3, linewidth = 0.5,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c("Medium smoke" = "#fc8d59",
                                    "High smoke" = "#d73027"),
                        name = "Smoke exposure") +
    labs(
      x = "Rate ratio vs low-smoke SA4s (95% CI)",
      y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "fig3_bushfire.png"),
         fig3, width = 8, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "fig3_bushfire.pdf"),
         fig3, width = 8, height = 4)
  cat("  -> Saved: fig3_bushfire.png/pdf\n")
} else {
  cat("  bushfire_did_results.csv not found, skipping.\n")
}


# ==============================================================================
# FIGURE 4: Sensitivity Analysis Summary Forest Plot
# ==============================================================================

cat("Figure 4: Sensitivity summary...\n")

sens_rows <- list()

# Primary results
dlnm_file <- file.path(tab_dir, "dlnm_results_summary.csv")
if (file.exists(dlnm_file)) {
  dlnm <- read.csv(dlnm_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total")
  sens_rows[["Primary (Tmax, SA4 weekly)"]] <- data.frame(
    spec = "Primary (Tmax, SA4 weekly)",
    rr = dlnm$rr_p95, rr_lo = dlnm$rr_p95_lo, rr_hi = dlnm$rr_p95_hi
  )
}

# Tmin alternative
tmin_file <- file.path(tab_dir, "tmin_results_summary.csv")
if (file.exists(tmin_file)) {
  tmin <- read.csv(tmin_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total")
  sens_rows[["Alternative: Tmin"]] <- data.frame(
    spec = "Alternative: Tmin",
    rr = tmin$rr_p95, rr_lo = tmin$rr_p95_lo, rr_hi = tmin$rr_p95_hi
  )
}

# LGA monthly
lga_file <- file.path(tab_dir, "lga_sensitivity_rr.csv")
if (file.exists(lga_file)) {
  lga <- read.csv(lga_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total")
  if (nrow(lga) > 0 && "rr_p95" %in% names(lga)) {
    sens_rows[["LGA monthly resolution"]] <- data.frame(
      spec = "LGA monthly resolution",
      rr = lga$rr_p95, rr_lo = lga$rr_p95_lo, rr_hi = lga$rr_p95_hi
    )
  }
}

# Influence-trimmed
infl_file <- file.path(tab_dir, "influence_sensitivity.csv")
if (file.exists(infl_file)) {
  infl <- read.csv(infl_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total")
  sens_rows[["Excluding influential obs"]] <- data.frame(
    spec = "Excluding influential obs",
    rr = infl$rr_trim, rr_lo = infl$rr_trim_lo, rr_hi = infl$rr_trim_hi
  )
}

# Lag constraint — unconstrained
lag_file <- file.path(tab_dir, "lag_constraint_sensitivity.csv")
if (file.exists(lag_file)) {
  lag_unc <- read.csv(lag_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total", spec_key == "unconstrained")
  if (nrow(lag_unc) > 0) {
    sens_rows[["Unconstrained lag"]] <- data.frame(
      spec = "Unconstrained lag",
      rr = lag_unc$rr_p95, rr_lo = lag_unc$rr_p95_lo, rr_hi = lag_unc$rr_p95_hi
    )
  }
}

# Fourier 3 pairs
four_file <- file.path(tab_dir, "fourier_sensitivity.csv")
if (file.exists(four_file)) {
  four <- read.csv(four_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total", fourier == "3 pairs")
  if (nrow(four) > 0) {
    sens_rows[["3 Fourier pairs"]] <- data.frame(
      spec = "3 Fourier pairs",
      rr = four$rr_p95, rr_lo = four$rr_p95_lo, rr_hi = four$rr_p95_hi
    )
  }
}

# Newey-West
nw_file <- file.path(tab_dir, "newey_west_comparison.csv")
if (file.exists(nw_file)) {
  nw <- read.csv(nw_file, stringsAsFactors = FALSE) |>
    filter(group == "resp_total")
  if (nrow(nw) > 0) {
    sens_rows[["Newey-West HAC SEs"]] <- data.frame(
      spec = "Newey-West HAC SEs",
      rr = nw$rr_hac, rr_lo = nw$rr_hac_lo, rr_hi = nw$rr_hac_hi
    )
  }
}

if (length(sens_rows) > 0) {
  sens_df <- bind_rows(sens_rows)
  sens_df$spec <- factor(sens_df$spec, levels = rev(sens_df$spec))

  fig4 <- ggplot(sens_df, aes(x = rr, y = spec)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi),
                   height = 0.3, linewidth = 0.5, colour = "steelblue") +
    geom_point(size = 2.5, colour = "steelblue") +
    labs(
      title = "Respiratory: sensitivity of RR at 95th percentile",
      x = "Cumulative RR (95% CI)",
      y = NULL
    ) +
    theme_pub() +
    theme(plot.title = element_text(size = 12))

  ggsave(file.path(fig_dir, "fig4_sensitivity_forest.png"),
         fig4, width = 9, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, "fig4_sensitivity_forest.pdf"),
         fig4, width = 9, height = 5)
  cat("  -> Saved: fig4_sensitivity_forest.png/pdf\n")
}


# ==============================================================================
# FIGURE S1: Lag-Response Curves
# ==============================================================================

cat("Figure S1: Lag-response curves...\n")

lag_df <- bind_rows(lag_data)
lag_df$group <- factor(lag_df$group, levels = c("Cardiovascular", "Mental Health", "Respiratory"))

figS1 <- ggplot(lag_df, aes(x = lag, y = rr)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi, fill = group), alpha = 0.15) +
  geom_line(aes(colour = group), linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ group, ncol = 3) +
  scale_colour_manual(values = pal_primary) +
  scale_fill_manual(values = pal_primary) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  labs(
    x = "Lag (weeks)",
    y = "RR per lag week at 95th percentile"
  ) +
  theme_pub() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "cm"))

ggsave(file.path(fig_dir, "figS1_lag_response.png"),
       figS1, width = 13, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "figS1_lag_response.pdf"),
       figS1, width = 13, height = 4.5)
cat("  -> Saved: figS1_lag_response.png/pdf\n")


# ==============================================================================
# FIGURE S2: Remoteness Stratification
# ==============================================================================

cat("Figure S2: Remoteness stratification...\n")

remote_file <- file.path(tab_dir, "remoteness_stratified_rr.csv")
if (file.exists(remote_file)) {
  rem <- read.csv(remote_file, stringsAsFactors = FALSE) |>
    filter(group %in% primary_groups) |>
    mutate(label = group_labels[group],
           label = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory")))

  figS2 <- ggplot(rem, aes(x = rr_p95, y = remoteness, colour = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5) +
    geom_point(size = 2.5) +
    facet_wrap(~ label, ncol = 3, scales = "free_x") +
    scale_colour_manual(values = pal_primary) +
    labs(
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "none",
          panel.spacing = unit(1, "cm"))

  ggsave(file.path(fig_dir, "figS2_remoteness.png"),
         figS2, width = 13, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "figS2_remoteness.pdf"),
         figS2, width = 13, height = 4)
  cat("  -> Saved: figS2_remoteness.png/pdf\n")
}


# ==============================================================================
# FIGURE S3: Negative Control
# ==============================================================================

cat("Figure S3: Negative control...\n")

nc_file <- file.path(tab_dir, "negative_control_results.csv")
if (file.exists(nc_file)) {
  nc <- read.csv(nc_file, stringsAsFactors = FALSE)

  role_colours <- c("Negative control" = "#4575b4", "Positive control" = "#d73027")

  figS3 <- ggplot(nc, aes(x = rr_p95, y = reorder(label, rr_p95), colour = role)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_colour_manual(values = role_colours, name = "Role") +
    labs(
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "figS3_negative_control.png"),
         figS3, width = 8, height = 3.5, dpi = 300)
  ggsave(file.path(fig_dir, "figS3_negative_control.pdf"),
         figS3, width = 8, height = 3.5)
  cat("  -> Saved: figS3_negative_control.png/pdf\n")
}


# ==============================================================================
# FIGURE S4: Temporal Adaptation
# ==============================================================================

cat("Figure S4: Temporal adaptation...\n")

temp_file <- file.path(tab_dir, "temporal_adaptation_rr.csv")
if (file.exists(temp_file)) {
  temp_ad <- read.csv(temp_file, stringsAsFactors = FALSE) |>
    filter(group %in% primary_groups) |>
    mutate(label = group_labels[group],
           label = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory")))

  period_colours <- c("steelblue", "#d73027")
  names(period_colours) <- sort(unique(temp_ad$period))

  figS4 <- ggplot(temp_ad, aes(x = rr_p95, y = period, colour = period)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.6) +
    geom_point(size = 3) +
    facet_wrap(~ label, ncol = 3, scales = "free_x") +
    scale_colour_manual(values = period_colours, name = "Period") +
    labs(
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = NULL
    ) +
    theme_pub() +
    theme(legend.position = "bottom",
          panel.spacing = unit(1, "cm"))

  ggsave(file.path(fig_dir, "figS4_temporal.png"),
         figS4, width = 13, height = 4, dpi = 300)
  ggsave(file.path(fig_dir, "figS4_temporal.pdf"),
         figS4, width = 13, height = 4)
  cat("  -> Saved: figS4_temporal.png/pdf\n")
}


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 08 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Publication figures generated:\n")
cat("  Main text:\n")
cat("    fig1_exposure_response.png/pdf — ER curves (3 groups)\n")
cat("    fig2_equity_forest.png/pdf — SEIFA-stratified forest plot\n")
cat("    fig3_bushfire.png/pdf — Bushfire DiD forest plot\n")
cat("    fig4_sensitivity_forest.png/pdf — Sensitivity analysis summary\n")
cat("  Supplementary:\n")
cat("    figS1_lag_response.png/pdf — Lag-response at p95\n")
cat("    figS2_remoteness.png/pdf — Remoteness stratification\n")
cat("    figS3_negative_control.png/pdf — Negative control comparison\n")
cat("    figS4_temporal.png/pdf — Temporal adaptation\n")
