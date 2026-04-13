# ==============================================================================
# Script 07: Black Summer Bushfire Analysis (Difference-in-Differences)
# ==============================================================================
# Examines whether the 2019-20 Black Summer bushfires produced detectable
# changes in respiratory medication dispensing in smoke-affected SA4 regions.
#
# Design:
#   - Exposure: SA4s classified as high/medium/low smoke based on mean PM2.5
#     during Oct 2019 - Feb 2020 (vs baseline Oct 2018 - Feb 2019)
#   - Outcome: respiratory medication (R03) weekly prescription count
#   - Method: difference-in-differences comparing smoke-exposed vs unexposed
#     areas, pre-bushfire vs during-bushfire periods
#   - Controls: temperature, rainfall, seasonal trends, SA4 fixed effects
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/processed/pm25_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/bushfire_did_results.csv
#   outputs/tables/bushfire_smoke_zones.csv
#   outputs/figures/bushfire_time_series.png
#   outputs/figures/bushfire_did_forest.png
#   outputs/figures/bushfire_parallel_trends.png
#   outputs/tables/bushfire_parallel_trends_wald.csv
#
# Dependencies:
#   install.packages(c("dplyr", "lubridate", "ggplot2", "patchwork",
#                       "fixest", "broom"))
# ==============================================================================

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

cat("=" |> strrep(70), "\n")
cat("Script 07: Black Summer Bushfire DiD Analysis\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA AND CLASSIFY SMOKE EXPOSURE ZONES
# ==============================================================================

cat("Loading data...\n")

# Full panel
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

# PM2.5 data
pm25 <- read.csv(file.path(data_dir, "pm25_weekly_sa4.csv"),
                 stringsAsFactors = FALSE)
pm25$week_start   <- as.Date(pm25$week_start)
pm25$region_code  <- as.character(pm25$region_code)

cat("  Panel rows:", format(nrow(panel), big.mark = ","), "\n")
cat("  PM2.5 rows:", format(nrow(pm25), big.mark = ","), "\n")

# --- Compute smoke exposure during Black Summer ---
cat("\nClassifying smoke exposure zones...\n")

# Black Summer period: Oct 2019 - Feb 2020
bs_start <- as.Date("2019-10-01")
bs_end   <- as.Date("2020-02-28")

# Baseline period: same months, previous year
bl_start <- as.Date("2018-10-01")
bl_end   <- as.Date("2019-02-28")

# Mean weekly PM2.5 during Black Summer
bs_pm25 <- pm25 |>
  filter(week_start >= bs_start, week_start <= bs_end) |>
  group_by(region_code) |>
  summarise(
    bs_pm25_mean = mean(pm25_mean, na.rm = TRUE),
    bs_pm25_max  = max(pm25_max, na.rm = TRUE),
    bs_n_weeks   = n(),
    .groups = "drop"
  )

# Mean weekly PM2.5 during baseline
bl_pm25 <- pm25 |>
  filter(week_start >= bl_start, week_start <= bl_end) |>
  group_by(region_code) |>
  summarise(
    bl_pm25_mean = mean(pm25_mean, na.rm = TRUE),
    .groups = "drop"
  )

# Merge and compute excess PM2.5
smoke <- bs_pm25 |>
  left_join(bl_pm25, by = "region_code") |>
  mutate(
    pm25_excess = bs_pm25_mean - coalesce(bl_pm25_mean, 7),  # default baseline ~7
    pm25_ratio  = bs_pm25_mean / pmax(coalesce(bl_pm25_mean, 7), 1)
  )

# Classify smoke zones using tertiles of Black Summer mean PM2.5
# High: top tertile, Medium: middle, Low: bottom
smoke <- smoke |>
  mutate(
    smoke_zone = case_when(
      bs_pm25_mean >= quantile(bs_pm25_mean, 2/3) ~ "high",
      bs_pm25_mean >= quantile(bs_pm25_mean, 1/3) ~ "medium",
      TRUE ~ "low"
    ),
    smoke_zone = factor(smoke_zone, levels = c("low", "medium", "high"))
  )

cat("  SA4s classified:", nrow(smoke), "\n")
cat("  Smoke zone distribution:\n")
print(table(smoke$smoke_zone))
cat("\n  Mean PM2.5 by zone:\n")
smoke |>
  group_by(smoke_zone) |>
  summarise(
    n = n(),
    pm25_mean = round(mean(bs_pm25_mean), 1),
    pm25_max  = round(max(bs_pm25_max), 0),
    excess    = round(mean(pm25_excess), 1),
    .groups = "drop"
  ) |>
  print()

# Save smoke zone classification
write.csv(smoke,
          file.path(tab_dir, "bushfire_smoke_zones.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/bushfire_smoke_zones.csv\n")


# ==============================================================================
# 2. BUILD DID ANALYSIS DATASET
# ==============================================================================

cat("\nBuilding DiD dataset...\n")

# Study window: 1 year before to end of Black Summer
# Pre-period:    Oct 2018 - Sep 2019 (12 months before)
# Post-period:   Oct 2019 - Feb 2020 (Black Summer)
did_start <- as.Date("2018-10-01")
did_end   <- as.Date("2020-02-28")

# Filter to respiratory medications
resp_groups <- c("resp_total", "resp_relievers", "resp_preventers", "resp_copd_other")

did_data <- panel |>
  filter(
    analysis_group %in% resp_groups,
    week_start >= did_start,
    week_start <= did_end,
    suppressed != "True",
    !is.na(count)
  ) |>
  left_join(smoke |> select(region_code, smoke_zone, bs_pm25_mean),
            by = c("sa4_code" = "region_code")) |>
  filter(!is.na(smoke_zone)) |>
  mutate(
    # DiD indicators
    post         = as.integer(week_start >= bs_start),
    treated      = as.integer(smoke_zone == "high"),
    treated_med  = as.integer(smoke_zone %in% c("high", "medium")),
    did          = post * treated,
    did_med      = post * treated_med,
    # Time controls
    week_of_year = isoweek(week_start),
    time_index   = as.numeric(week_start - did_start) / 7,
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    # Log outcome for semi-elasticity
    log_count    = log(pmax(count, 1))
  )

cat("  DiD dataset rows:", format(nrow(did_data), big.mark = ","), "\n")
cat("  SA4 regions:", n_distinct(did_data$sa4_code), "\n")
cat("  Treated (high smoke):", sum(did_data$treated == 1 & did_data$analysis_group == "resp_total") /
    n_distinct(did_data$week_start[did_data$treated == 1]), "SA4-weeks\n")
cat("  Pre-period weeks:", n_distinct(did_data$week_start[did_data$post == 0]), "\n")
cat("  Post-period weeks:", n_distinct(did_data$week_start[did_data$post == 1]), "\n")


# ==============================================================================
# 3. FIT DID MODELS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting DiD models\n")
cat("=" |> strrep(70), "\n")

resp_labels <- c(
  resp_total      = "All Respiratory",
  resp_relievers  = "Relievers",
  resp_preventers = "Preventers",
  resp_copd_other = "COPD/Other"
)

did_results <- list()

for (grp in resp_groups) {
  cat("\n--- ", resp_labels[grp], " ---\n")

  df <- did_data |> filter(analysis_group == grp)
  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

  if (nrow(df) < 100) {
    cat("  Too few observations, skipping.\n")
    next
  }

  # --- Model 1: High smoke vs low (binary treatment) ---
  cat("  Model 1: High vs Low smoke zones...\n")
  df_hl <- df |> filter(smoke_zone %in% c("high", "low"))

  mod1 <- glm(
    count ~ treated * post +
      sin1 + cos1 +
      tmax_mean + precip_total +
      factor(sa4_code),
    data = df_hl,
    family = quasipoisson(link = "log")
  )

  coef_did <- summary(mod1)$coefficients
  if ("treated:post" %in% rownames(coef_did)) {
    est  <- coef_did["treated:post", "Estimate"]
    se   <- coef_did["treated:post", "Std. Error"]
    pval <- coef_did["treated:post", "Pr(>|t|)"]
    rr   <- exp(est)
    rr_lo <- exp(est - 1.96 * se)
    rr_hi <- exp(est + 1.96 * se)
    cat("  DiD (high vs low): RR =", round(rr, 4),
        "(", round(rr_lo, 4), "-", round(rr_hi, 4), ")",
        "p =", format.pval(pval, digits = 3), "\n")
  }

  # --- Model 2: Three-level smoke zones (low = ref) ---
  cat("  Model 2: Three-level smoke zones...\n")
  df$smoke_zone <- relevel(df$smoke_zone, ref = "low")

  mod2 <- glm(
    count ~ smoke_zone * post +
      sin1 + cos1 +
      tmax_mean + precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  coef2 <- summary(mod2)$coefficients

  # Extract interaction terms
  results_row <- data.frame(
    group = grp,
    label = resp_labels[grp],
    stringsAsFactors = FALSE
  )

  for (zone in c("medium", "high")) {
    term <- paste0("smoke_zone", zone, ":post")
    if (term %in% rownames(coef2)) {
      est  <- coef2[term, "Estimate"]
      se   <- coef2[term, "Std. Error"]
      pval <- coef2[term, "Pr(>|t|)"]
      results_row[[paste0("rr_", zone)]]    <- round(exp(est), 4)
      results_row[[paste0("rr_", zone, "_lo")]] <- round(exp(est - 1.96 * se), 4)
      results_row[[paste0("rr_", zone, "_hi")]] <- round(exp(est + 1.96 * se), 4)
      results_row[[paste0("p_", zone)]]     <- pval
      cat("    ", zone, " vs low: RR =", round(exp(est), 4),
          "(", round(exp(est - 1.96 * se), 4), "-", round(exp(est + 1.96 * se), 4), ")",
          "p =", format.pval(pval, digits = 3), "\n")
    }
  }

  results_row$n_obs      <- nrow(df)
  results_row$n_sa4      <- n_distinct(df$sa4_code)
  results_row$dispersion <- round(summary(mod2)$dispersion, 2)

  did_results[[grp]] <- results_row
}

# Save results
did_table <- bind_rows(did_results)
write.csv(did_table,
          file.path(tab_dir, "bushfire_did_results.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/bushfire_did_results.csv\n")


# ==============================================================================
# 4. PARALLEL TRENDS CHECK
# ==============================================================================

cat("\nChecking parallel trends assumption...\n")

# Compute weekly mean prescriptions by smoke zone (resp_total only)
trends_data <- did_data |>
  filter(analysis_group == "resp_total") |>
  group_by(smoke_zone, week_start) |>
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    .groups = "drop"
  )

p_trends <- ggplot(trends_data, aes(x = week_start, y = mean_count,
                                     colour = smoke_zone)) +
  geom_line(alpha = 0.6, linewidth = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1) +
  geom_vline(xintercept = bs_start, linetype = "dashed", colour = "red") +
  annotate("text", x = bs_start + 7, y = Inf, label = "Black Summer onset",
           hjust = 0, vjust = 1.5, colour = "red", size = 3.5) +
  annotate("rect", xmin = bs_start, xmax = bs_end,
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.05) +
  scale_colour_manual(
    values = c(low = "#4575b4", medium = "#fee090", high = "#d73027"),
    labels = c(low = "Low smoke", medium = "Medium smoke", high = "High smoke"),
    name = "Smoke Zone"
  ) +
  labs(
    title = "Respiratory prescriptions by smoke exposure zone",
    subtitle = "Weekly mean prescription count per SA4",
    x = NULL,
    y = "Mean weekly prescriptions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "bushfire_parallel_trends.png"),
       p_trends, width = 10, height = 6, dpi = 300)
cat("  -> Saved: bushfire_parallel_trends.png\n")


# ==============================================================================
# 5. INDEXED TIME SERIES (normalised to pre-period)
# ==============================================================================

cat("Generating indexed time series...\n")

# Normalise each zone to its pre-period mean
pre_means <- trends_data |>
  filter(week_start < bs_start) |>
  group_by(smoke_zone) |>
  summarise(pre_mean = mean(mean_count), .groups = "drop")

trends_indexed <- trends_data |>
  left_join(pre_means, by = "smoke_zone") |>
  mutate(indexed = mean_count / pre_mean * 100)

p_indexed <- ggplot(trends_indexed, aes(x = week_start, y = indexed,
                                         colour = smoke_zone)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 100, linetype = "dotted", colour = "grey50") +
  geom_vline(xintercept = bs_start, linetype = "dashed", colour = "red") +
  annotate("rect", xmin = bs_start, xmax = bs_end,
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.05) +
  scale_colour_manual(
    values = c(low = "#4575b4", medium = "#fee090", high = "#d73027"),
    labels = c(low = "Low smoke", medium = "Medium smoke", high = "High smoke"),
    name = "Smoke Zone"
  ) +
  labs(
    title = "Respiratory prescriptions indexed to pre-bushfire baseline",
    subtitle = "Pre-period mean = 100",
    x = NULL,
    y = "Index (pre-period = 100)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "bushfire_time_series.png"),
       p_indexed, width = 10, height = 6, dpi = 300)
cat("  -> Saved: bushfire_time_series.png\n")


# ==============================================================================
# 6. FOREST PLOT OF DID EFFECTS
# ==============================================================================

cat("Generating DiD forest plot...\n")

# Reshape results for forest plot
forest_rows <- list()
for (i in seq_len(nrow(did_table))) {
  r <- did_table[i, ]
  for (zone in c("medium", "high")) {
    rr_col <- paste0("rr_", zone)
    if (rr_col %in% names(r) && !is.na(r[[rr_col]])) {
      forest_rows[[length(forest_rows) + 1]] <- data.frame(
        label     = r$label,
        zone      = zone,
        rr        = r[[paste0("rr_", zone)]],
        rr_lo     = r[[paste0("rr_", zone, "_lo")]],
        rr_hi     = r[[paste0("rr_", zone, "_hi")]],
        p_val     = r[[paste0("p_", zone)]],
        stringsAsFactors = FALSE
      )
    }
  }
}

forest_df <- bind_rows(forest_rows) |>
  mutate(
    zone_label = ifelse(zone == "high", "High smoke", "Medium smoke"),
    group_zone = paste0(label, " — ", zone_label),
    sig = ifelse(p_val < 0.05, "*", "")
  )

p_forest <- ggplot(forest_df, aes(x = rr, y = reorder(group_zone, rr),
                                    colour = zone)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi), height = 0.25, linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sig), hjust = -1, size = 5, show.legend = FALSE) +
  scale_colour_manual(
    values = c(medium = "#fc8d59", high = "#d73027"),
    labels = c(medium = "Medium smoke", high = "High smoke"),
    name = "Smoke Zone"
  ) +
  labs(
    title = "Black Summer bushfire effect on respiratory prescriptions",
    subtitle = "Difference-in-differences RR vs low-smoke SA4s (reference)",
    x = "Rate Ratio (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "bushfire_did_forest.png"),
       p_forest, width = 9, height = 6, dpi = 300)
cat("  -> Saved: bushfire_did_forest.png\n")


# ==============================================================================
# 7. SENSITIVITY: Event-study specification
# ==============================================================================

cat("\nRunning event-study sensitivity analysis...\n")

# Create monthly periods relative to Black Summer onset
es_data <- did_data |>
  filter(analysis_group == "resp_total") |>
  mutate(
    months_since = as.numeric(interval(bs_start, week_start) %/% months(1)),
    # Bin into monthly periods
    period = case_when(
      months_since < -9  ~ -10,
      months_since > 4   ~ 5,
      TRUE               ~ months_since
    ),
    period_f = factor(period),
    treated_f = factor(treated)
  )

# Relevel: period -1 (month before onset) as reference
es_data$period_f <- relevel(es_data$period_f, ref = "-1")

mod_es <- glm(
  count ~ treated_f * period_f +
    sin1 + cos1 +
    tmax_mean + precip_total +
    factor(sa4_code),
  data = es_data,
  family = quasipoisson(link = "log")
)

# Extract interaction coefficients
es_coef <- summary(mod_es)$coefficients
int_terms <- grep("^treated_f1:period_f", rownames(es_coef), value = TRUE)

if (length(int_terms) > 0) {
  es_results <- data.frame(
    term   = int_terms,
    period = as.numeric(gsub("treated_f1:period_f", "", int_terms)),
    est    = es_coef[int_terms, "Estimate"],
    se     = es_coef[int_terms, "Std. Error"]
  ) |>
    mutate(
      rr    = exp(est),
      rr_lo = exp(est - 1.96 * se),
      rr_hi = exp(est + 1.96 * se)
    )

  # Add reference period
  es_results <- bind_rows(
    es_results,
    data.frame(term = "ref", period = -1, est = 0, se = 0,
               rr = 1, rr_lo = 1, rr_hi = 1)
  ) |>
    arrange(period)

  p_es <- ggplot(es_results, aes(x = period, y = rr)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "#d73027", alpha = 0.15) +
    geom_line(colour = "#d73027", linewidth = 0.8) +
    geom_point(colour = "#d73027", size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red", alpha = 0.5) +
    annotate("text", x = 0.2, y = Inf, label = "Black Summer\nonset",
             hjust = 0, vjust = 1.2, colour = "red", size = 3) +
    scale_x_continuous(breaks = seq(-10, 5, 2)) +
    labs(
      title = "Event-study: High-smoke SA4s vs low-smoke (respiratory prescriptions)",
      subtitle = "Monthly interaction coefficients, reference = 1 month pre-onset",
      x = "Months relative to Black Summer onset (Oct 2019)",
      y = "Rate Ratio"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(fig_dir, "bushfire_event_study.png"),
         p_es, width = 10, height = 6, dpi = 300)
  cat("  -> Saved: bushfire_event_study.png\n")

  # --- Joint Wald test of pre-period coefficients (parallel trends) ---
  # Under the parallel trends assumption, all pre-period interaction
  # coefficients should be jointly zero. We use an F-test (not chi-square)
  # because the quasi-Poisson dispersion parameter must be accounted for.

  pre_terms <- es_results |>
    filter(period < -1, term != "ref") |>
    pull(term)

  if (length(pre_terms) > 0) {
    # Extract coefficient vector and variance-covariance submatrix
    beta_pre <- coef(mod_es)[pre_terms]
    vcov_pre <- vcov(mod_es)[pre_terms, pre_terms]
    k <- length(pre_terms)
    n_obs <- nrow(es_data)
    n_params <- length(coef(mod_es))

    # F-statistic: (beta' %*% V^{-1} %*% beta) / k
    F_stat <- as.numeric(t(beta_pre) %*% solve(vcov_pre) %*% beta_pre) / k
    df_resid <- n_obs - n_params
    p_wald <- pf(F_stat, df1 = k, df2 = df_resid, lower.tail = FALSE)

    cat("\n  Joint Wald test of pre-period parallel trends:\n")
    cat("    Pre-period terms tested:", k, "\n")
    cat("    F-statistic:", round(F_stat, 3), "\n")
    cat("    df1 =", k, ", df2 =", df_resid, "\n")
    cat("    p-value:", format.pval(p_wald, digits = 3), "\n")
    cat("    Interpretation:",
        ifelse(p_wald > 0.05,
               "No evidence against parallel trends (fail to reject H0)\n",
               "Evidence against parallel trends (reject H0)\n"))

    # Save results
    wald_result <- data.frame(
      test = "Joint Wald F-test of pre-period interaction coefficients",
      n_pre_terms = k,
      F_stat = round(F_stat, 3),
      df1 = k,
      df2 = df_resid,
      p_value = p_wald
    )
    write.csv(wald_result,
              file.path(tab_dir, "bushfire_parallel_trends_wald.csv"),
              row.names = FALSE)
    cat("  -> Saved: bushfire_parallel_trends_wald.csv\n")
  }
}


# ==============================================================================
# 8. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 07 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("DiD results (vs low-smoke SA4s):\n")
did_table |>
  select(label, rr_medium, rr_medium_lo, rr_medium_hi, p_medium,
         rr_high, rr_high_lo, rr_high_hi, p_high) |>
  print(row.names = FALSE)

cat("\nSmoke zone summary:\n")
smoke |>
  group_by(smoke_zone) |>
  summarise(n = n(), mean_pm25 = round(mean(bs_pm25_mean), 1), .groups = "drop") |>
  print()

cat("\nOutputs:\n")
cat("  outputs/tables/bushfire_did_results.csv\n")
cat("  outputs/tables/bushfire_smoke_zones.csv\n")
cat("  outputs/tables/bushfire_parallel_trends_wald.csv\n")
cat("  outputs/figures/bushfire_parallel_trends.png\n")
cat("  outputs/figures/bushfire_time_series.png\n")
cat("  outputs/figures/bushfire_did_forest.png\n")
cat("  outputs/figures/bushfire_event_study.png\n")

cat("\nNext step: Script 08 — Publication visualisation and maps\n")
