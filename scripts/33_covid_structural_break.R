# ==============================================================================
# Script 33: COVID-19 Structural Break Analysis
# ==============================================================================
# Addresses the concern that COVID-19 lockdowns and telehealth shifts created
# a structural break in the dispensing time series (March 2020 onwards).
#
# Two approaches:
#   1. Pandemic interaction terms: add COVID period dummy × cross-basis interaction
#      to formally model the interrupted time series
#   2. Pre-pandemic sensitivity: rerun primary models on 2013–2019 data only
#      to confirm findings are not "COVID artifacts"
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/covid_prepandemic_rr.csv
#   outputs/tables/covid_interaction_results.csv
#   outputs/figures/covid_prepandemic_forest.png
#   outputs/figures/covid_er_comparison.png
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

primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

# COVID period definition
covid_start <- as.Date("2020-03-01")
covid_end   <- as.Date("2021-12-31")

# Pre-pandemic cutoff
prepandemic_end <- as.Date("2019-12-31")

cat("=" |> strrep(70), "\n")
cat("Script 33: COVID-19 Structural Break Analysis\n")
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
    cos2 = cos(4 * pi * week_of_year / 52),
    # COVID dummy
    covid_period = as.integer(week_start >= covid_start & week_start <= covid_end),
    # Pre-pandemic flag
    pre_pandemic = week_start <= prepandemic_end
  )

n_covid <- sum(panel$covid_period == 1 & panel$analysis_group == "resp_total", na.rm = TRUE)
n_pre   <- sum(panel$pre_pandemic & panel$analysis_group == "resp_total", na.rm = TRUE)
cat("  COVID period rows (resp_total):", format(n_covid, big.mark = ","), "\n")
cat("  Pre-pandemic rows (resp_total):", format(n_pre, big.mark = ","), "\n")


# ==============================================================================
# 2. APPROACH 1: COVID INTERACTION TERMS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Approach 1: COVID interaction terms\n")
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
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  # --- Model A: Original (no COVID term) ---
  cat("  Model A: No COVID term...\n")
  mod_a <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )
  pred_a <- crosspred(cb, mod_a, at = c(temp_p95, temp_p99), cen = temp_ref)

  # --- Model B: COVID dummy (level shift only) ---
  cat("  Model B: COVID dummy...\n")
  mod_b <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      covid_period +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )
  pred_b <- crosspred(cb, mod_b, at = c(temp_p95, temp_p99), cen = temp_ref)

  # COVID level shift coefficient
  covid_coef <- coef(mod_b)["covid_period"]
  covid_se   <- summary(mod_b)$coefficients["covid_period", "Std. Error"]
  covid_pval <- summary(mod_b)$coefficients["covid_period", "Pr(>|t|)"]

  # --- Model C: COVID × temperature interaction ---
  cat("  Model C: COVID × temperature interaction...\n")
  mod_c <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      covid_period + tmax_mean:covid_period +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )
  pred_c <- crosspred(cb, mod_c, at = c(temp_p95, temp_p99), cen = temp_ref)

  # Interaction test: model C vs model B
  dev_b <- deviance(mod_b)
  dev_c <- deviance(mod_c)
  df_diff <- mod_b$df.residual - mod_c$df.residual
  disp_c  <- summary(mod_c)$dispersion
  f_stat  <- ((dev_b - dev_c) / df_diff) / disp_c
  int_pval <- pf(f_stat, df_diff, mod_c$df.residual, lower.tail = FALSE)

  cat("  COVID level shift: RR =", round(exp(covid_coef), 4),
      "  p =", format.pval(covid_pval, digits = 3), "\n")
  cat("  COVID × temp interaction: F =", round(f_stat, 3),
      "  p =", format.pval(int_pval, digits = 3), "\n")

  interaction_results[[grp]] <- data.frame(
    group = grp,
    label = group_labels[grp],
    # COVID level shift
    covid_rr      = round(exp(covid_coef), 4),
    covid_rr_lo   = round(exp(covid_coef - 1.96 * covid_se), 4),
    covid_rr_hi   = round(exp(covid_coef + 1.96 * covid_se), 4),
    covid_p       = covid_pval,
    # Interaction test
    int_f_stat    = round(f_stat, 3),
    int_p         = int_pval,
    # RR@p95 from model A (no COVID)
    rr_p95_nocovid    = round(pred_a$allRRfit[1], 4),
    rr_p95_nocovid_lo = round(pred_a$allRRlow[1], 4),
    rr_p95_nocovid_hi = round(pred_a$allRRhigh[1], 4),
    # RR@p95 from model B (COVID dummy)
    rr_p95_coviddum    = round(pred_b$allRRfit[1], 4),
    rr_p95_coviddum_lo = round(pred_b$allRRlow[1], 4),
    rr_p95_coviddum_hi = round(pred_b$allRRhigh[1], 4),
    # RR@p95 from model C (COVID interaction)
    rr_p95_covidint    = round(pred_c$allRRfit[1], 4),
    rr_p95_covidint_lo = round(pred_c$allRRlow[1], 4),
    rr_p95_covidint_hi = round(pred_c$allRRhigh[1], 4),
    n_obs = nrow(df),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. APPROACH 2: PRE-PANDEMIC SENSITIVITY ANALYSIS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Approach 2: Pre-pandemic (2013-2019) sensitivity analysis\n")
cat("=" |> strrep(70), "\n")

prepandemic_results <- list()
prepandemic_preds   <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  # Full period model
  df_full <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  # Pre-pandemic only
  df_pre <- df_full |>
    filter(pre_pandemic)

  cat("  Full period:", format(nrow(df_full), big.mark = ","), "obs\n")
  cat("  Pre-pandemic:", format(nrow(df_pre), big.mark = ","), "obs\n")

  # --- Full period model ---
  temp_knots <- quantile(df_full$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb_full <- crossbasis(
    df_full$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years_full <- n_distinct(df_full$year)
  trend_df_full <- max(2, round(n_years_full * 2))
  temp_ref <- median(df_full$tmax_mean, na.rm = TRUE)
  temp_range <- range(df_full$tmax_mean, na.rm = TRUE)

  mod_full <- glm(
    count ~ cb_full +
      ns(time_index, df = trend_df_full) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df_full,
    family = quasipoisson(link = "log")
  )

  pred_full <- crosspred(cb_full, mod_full,
                          at = seq(temp_range[1], temp_range[2], length.out = 80),
                          cen = temp_ref)
  pred_full_pctl <- crosspred(cb_full, mod_full,
                               at = c(quantile(df_full$tmax_mean, c(0.95, 0.99))),
                               cen = temp_ref)

  # --- Pre-pandemic model ---
  temp_knots_pre <- quantile(df_pre$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb_pre <- crossbasis(
    df_pre$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots_pre),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years_pre <- n_distinct(df_pre$year)
  trend_df_pre <- max(2, round(n_years_pre * 2))
  temp_ref_pre <- median(df_pre$tmax_mean, na.rm = TRUE)
  temp_range_pre <- range(df_pre$tmax_mean, na.rm = TRUE)

  mod_pre <- glm(
    count ~ cb_pre +
      ns(time_index, df = trend_df_pre) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df_pre,
    family = quasipoisson(link = "log")
  )

  pred_pre <- crosspred(cb_pre, mod_pre,
                         at = seq(temp_range_pre[1], temp_range_pre[2], length.out = 80),
                         cen = temp_ref_pre)
  temp_p95_pre <- quantile(df_pre$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99_pre <- quantile(df_pre$tmax_mean, 0.99, na.rm = TRUE)
  pred_pre_pctl <- crosspred(cb_pre, mod_pre,
                              at = c(temp_p95_pre, temp_p99_pre),
                              cen = temp_ref_pre)

  cat("  Full period RR@p95:", round(pred_full_pctl$allRRfit[1], 4), "\n")
  cat("  Pre-pandemic RR@p95:", round(pred_pre_pctl$allRRfit[1], 4), "\n")

  # Store
  prepandemic_results[[grp]] <- data.frame(
    group = grp,
    label = group_labels[grp],
    # Full period
    full_rr_p95    = round(pred_full_pctl$allRRfit[1], 4),
    full_rr_p95_lo = round(pred_full_pctl$allRRlow[1], 4),
    full_rr_p95_hi = round(pred_full_pctl$allRRhigh[1], 4),
    full_rr_p99    = round(pred_full_pctl$allRRfit[2], 4),
    full_rr_p99_lo = round(pred_full_pctl$allRRlow[2], 4),
    full_rr_p99_hi = round(pred_full_pctl$allRRhigh[2], 4),
    # Pre-pandemic
    pre_rr_p95    = round(pred_pre_pctl$allRRfit[1], 4),
    pre_rr_p95_lo = round(pred_pre_pctl$allRRlow[1], 4),
    pre_rr_p95_hi = round(pred_pre_pctl$allRRhigh[1], 4),
    pre_rr_p99    = round(pred_pre_pctl$allRRfit[2], 4),
    pre_rr_p99_lo = round(pred_pre_pctl$allRRlow[2], 4),
    pre_rr_p99_hi = round(pred_pre_pctl$allRRhigh[2], 4),
    n_full = nrow(df_full),
    n_pre  = nrow(df_pre),
    disp_full = round(summary(mod_full)$dispersion, 2),
    disp_pre  = round(summary(mod_pre)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  # Store predictions for ER overlay
  prepandemic_preds[[grp]] <- list(
    full = pred_full, pre = pred_pre,
    temp_ref_full = temp_ref, temp_ref_pre = temp_ref_pre
  )
}


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

int_table <- bind_rows(interaction_results)
write.csv(int_table,
          file.path(tab_dir, "covid_interaction_results.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/covid_interaction_results.csv\n")

pre_table <- bind_rows(prepandemic_results)
write.csv(pre_table,
          file.path(tab_dir, "covid_prepandemic_rr.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/covid_prepandemic_rr.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# --- A. Forest plot: Full vs pre-pandemic RRs ---
forest_rows <- list()
for (i in seq_len(nrow(pre_table))) {
  r <- pre_table[i, ]
  forest_rows[[length(forest_rows) + 1]] <- data.frame(
    label = r$label, period = "Full period (2013-2022)",
    rr = r$full_rr_p95, rr_lo = r$full_rr_p95_lo, rr_hi = r$full_rr_p95_hi)
  forest_rows[[length(forest_rows) + 1]] <- data.frame(
    label = r$label, period = "Pre-pandemic (2013-2019)",
    rr = r$pre_rr_p95, rr_lo = r$pre_rr_p95_lo, rr_hi = r$pre_rr_p95_hi)
}

forest_df <- bind_rows(forest_rows) |>
  mutate(
    period = factor(period, levels = c("Full period (2013-2022)",
                                        "Pre-pandemic (2013-2019)")),
    label_period = paste0(label, "\n", period)
  )

p_forest <- ggplot(forest_df, aes(x = rr, y = label_period, colour = period)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi), height = 0.3, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("Full period (2013-2022)" = "#4575b4",
               "Pre-pandemic (2013-2019)" = "#d73027"),
    name = "Period"
  ) +
  labs(
    title = "Pre-pandemic sensitivity: cumulative RR at 95th percentile",
    subtitle = "Testing whether findings hold without COVID-era data",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "covid_prepandemic_forest.png"),
       p_forest, width = 10, height = 7, dpi = 300)
cat("  -> Saved: covid_prepandemic_forest.png\n")

# --- B. Overlaid ER curves ---
er_plots <- list()
for (grp in primary_groups) {
  preds <- prepandemic_preds[[grp]]
  lbl <- group_labels[grp]

  # Full period ER
  df_full <- data.frame(
    temp  = as.numeric(names(preds$full$allRRfit)),
    rr    = preds$full$allRRfit,
    rr_lo = preds$full$allRRlow,
    rr_hi = preds$full$allRRhigh,
    period = "Full (2013-2022)"
  )

  df_pre <- data.frame(
    temp  = as.numeric(names(preds$pre$allRRfit)),
    rr    = preds$pre$allRRfit,
    rr_lo = preds$pre$allRRlow,
    rr_hi = preds$pre$allRRhigh,
    period = "Pre-pandemic (2013-2019)"
  )

  df_plot <- bind_rows(df_full, df_pre)

  p <- ggplot(df_plot, aes(x = temp, y = rr, colour = period, fill = period)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(
      values = c("Full (2013-2022)" = "#4575b4",
                 "Pre-pandemic (2013-2019)" = "#d73027"),
      name = "Period"
    ) +
    scale_fill_manual(
      values = c("Full (2013-2022)" = "#4575b4",
                 "Pre-pandemic (2013-2019)" = "#d73027"),
      name = "Period"
    ) +
    labs(title = lbl, x = "Weekly mean Tmax (°C)", y = "Cumulative RR") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  er_plots[[grp]] <- p
}

combined <- wrap_plots(er_plots, ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "covid_er_comparison.png"),
       combined, width = 9, height = 14, dpi = 300)
cat("  -> Saved: covid_er_comparison.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 33 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("COVID interaction results:\n")
int_table |>
  select(label, covid_rr, covid_p, int_f_stat, int_p) |>
  print(row.names = FALSE)

cat("\nPre-pandemic vs full period comparison (RR@p95):\n")
pre_table |>
  select(label, full_rr_p95, full_rr_p95_lo, full_rr_p95_hi,
         pre_rr_p95, pre_rr_p95_lo, pre_rr_p95_hi) |>
  print(row.names = FALSE)

cat("\nInterpretation:\n")
cat("  - If pre-pandemic RR similar to full: findings NOT a COVID artifact\n")
cat("  - Significant COVID dummy: dispensing levels shifted during pandemic\n")
cat("  - Significant COVID × temp interaction: heat-dispensing relationship\n")
cat("    changed during pandemic (e.g., due to lockdowns, telehealth)\n")

cat("\nOutputs:\n")
cat("  outputs/tables/covid_interaction_results.csv\n")
cat("  outputs/tables/covid_prepandemic_rr.csv\n")
cat("  outputs/figures/covid_prepandemic_forest.png\n")
cat("  outputs/figures/covid_er_comparison.png\n")
