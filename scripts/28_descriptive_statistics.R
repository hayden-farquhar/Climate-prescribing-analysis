# ==============================================================================
# Script 28: Descriptive Statistics (Table 1)
# ==============================================================================
# Produces publication-ready descriptive statistics for the study population,
# exposures, and outcomes. Generates Table 1 for the manuscript.
#
# Sections:
#   A. Study period and spatial coverage
#   B. Outcome distributions by medication group
#   C. Temperature and climate exposure distributions
#   D. Heatwave frequency and duration
#   E. Socioeconomic and remoteness distributions
#   F. Bushfire period PM2.5 summary
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv
#   data/reference/remoteness_sa4.csv
#
# Outputs:
#   outputs/tables/table1_descriptive.csv
#   outputs/tables/descriptive_by_group.csv
#   outputs/tables/descriptive_climate.csv
#   outputs/figures/descriptive_outcome_distributions.png
#   outputs/figures/descriptive_temperature_distribution.png
#
# Dependencies:
#   install.packages(c("dplyr", "lubridate", "tidyr", "ggplot2", "patchwork"))
# ==============================================================================

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(patchwork)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("=" |> strrep(70), "\n")
cat("Script 28: Descriptive Statistics (Table 1)\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)
panel$year       <- year(panel$week_start)

# SEIFA quintiles
seifa_file <- file.path(ref_dir, "seifa_sa4_quintiles.csv")
if (file.exists(seifa_file)) {
  seifa <- read.csv(seifa_file, stringsAsFactors = FALSE)
  seifa$sa4_code <- as.character(seifa$sa4_code)
} else {
  seifa <- NULL
  cat("  SEIFA file not found.\n")
}

# Remoteness
remote_file <- file.path(ref_dir, "remoteness_sa4.csv")
if (file.exists(remote_file)) {
  remote <- read.csv(remote_file, stringsAsFactors = FALSE)
  remote$sa4_code <- as.character(remote$sa4_code)
} else {
  remote <- NULL
  cat("  Remoteness file not found.\n")
}


# ==============================================================================
# 2. STUDY OVERVIEW
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("A. Study Overview\n")
cat("=" |> strrep(50), "\n")

cat("  Study period:", as.character(min(panel$week_start)), "to",
    as.character(max(panel$week_start)), "\n")
cat("  Duration:", round(as.numeric(max(panel$week_start) - min(panel$week_start)) / 365.25, 1),
    "years\n")
cat("  Total weeks:", n_distinct(panel$week_start), "\n")
cat("  SA4 regions:", n_distinct(panel$sa4_code), "\n")
cat("  Jurisdictions:", paste(sort(unique(panel$jurisdiction)), collapse = ", "), "\n")
cat("  Total panel rows:", format(nrow(panel), big.mark = ","), "\n")
cat("  Medication groups:", n_distinct(panel$analysis_group), "\n")
cat("  Groups:", paste(sort(unique(panel$analysis_group)), collapse = ", "), "\n")


# ==============================================================================
# 3. OUTCOME DISTRIBUTIONS BY MEDICATION GROUP
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("B. Outcome Distributions\n")
cat("=" |> strrep(50), "\n")

group_labels <- c(
  cvd_total              = "Cardiovascular (total)",
  cvd_antithrombotic     = "CVD: Antithrombotic",
  cvd_bp_lowering        = "CVD: BP Lowering",
  cvd_lipid_modifying    = "CVD: Lipid Modifying",
  cvd_other              = "CVD: Other",
  mh_total               = "Mental Health (total)",
  mh_antidepressants     = "MH: Antidepressants",
  mh_anxiolytics         = "MH: Anxiolytics",
  mh_other               = "MH: Other",
  resp_total             = "Respiratory (total)",
  resp_relievers         = "Resp: Relievers",
  resp_preventers        = "Resp: Preventers",
  resp_copd_other        = "Resp: COPD/Other",
  resp_oral_corticosteroids = "Resp: Oral Corticosteroids"
)

outcome_stats <- panel |>
  filter(suppressed != "True", !is.na(count)) |>
  group_by(analysis_group) |>
  summarise(
    n_obs          = n(),
    n_sa4          = n_distinct(sa4_code),
    n_weeks        = n_distinct(week_start),
    total_rx       = sum(count, na.rm = TRUE),
    mean_count     = round(mean(count, na.rm = TRUE), 1),
    sd_count       = round(sd(count, na.rm = TRUE), 1),
    median_count   = median(count, na.rm = TRUE),
    iqr_lo         = quantile(count, 0.25, na.rm = TRUE),
    iqr_hi         = quantile(count, 0.75, na.rm = TRUE),
    min_count      = min(count, na.rm = TRUE),
    max_count      = max(count, na.rm = TRUE),
    pct_zero       = round(100 * mean(count == 0, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  mutate(label = group_labels[analysis_group]) |>
  arrange(match(analysis_group, names(group_labels)))

cat("\n")
print(as.data.frame(
  outcome_stats |>
    select(label, n_obs, total_rx, mean_count, sd_count, median_count, iqr_lo, iqr_hi, pct_zero)
), row.names = FALSE)


# ==============================================================================
# 4. CLIMATE EXPOSURE DISTRIBUTIONS
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("C. Climate Exposure Distributions\n")
cat("=" |> strrep(50), "\n")

# One row per SA4-week (use a single analysis group to avoid duplication)
climate_df <- panel |>
  filter(analysis_group == "cvd_total",
         !is.na(tmax_mean))

climate_stats <- data.frame(
  variable = c("Weekly mean Tmax (°C)", "Weekly max Tmax (°C)",
               "Weekly mean Tmin (°C)", "Weekly total precip (mm)",
               "Days above 90th pctl", "Days above 95th pctl",
               "Days above 99th pctl", "Heatwave days/week"),
  n        = c(sum(!is.na(climate_df$tmax_mean)),
               sum(!is.na(climate_df$tmax_max)),
               sum(!is.na(climate_df$tmin_mean)),
               sum(!is.na(climate_df$precip_total)),
               sum(!is.na(climate_df$days_above_p90)),
               sum(!is.na(climate_df$days_above_p95)),
               sum(!is.na(climate_df$days_above_p99)),
               sum(!is.na(climate_df$heatwave_days))),
  mean     = round(c(mean(climate_df$tmax_mean, na.rm = TRUE),
                      mean(climate_df$tmax_max, na.rm = TRUE),
                      mean(climate_df$tmin_mean, na.rm = TRUE),
                      mean(climate_df$precip_total, na.rm = TRUE),
                      mean(climate_df$days_above_p90, na.rm = TRUE),
                      mean(climate_df$days_above_p95, na.rm = TRUE),
                      mean(climate_df$days_above_p99, na.rm = TRUE),
                      mean(climate_df$heatwave_days, na.rm = TRUE)), 2),
  sd       = round(c(sd(climate_df$tmax_mean, na.rm = TRUE),
                      sd(climate_df$tmax_max, na.rm = TRUE),
                      sd(climate_df$tmin_mean, na.rm = TRUE),
                      sd(climate_df$precip_total, na.rm = TRUE),
                      sd(climate_df$days_above_p90, na.rm = TRUE),
                      sd(climate_df$days_above_p95, na.rm = TRUE),
                      sd(climate_df$days_above_p99, na.rm = TRUE),
                      sd(climate_df$heatwave_days, na.rm = TRUE)), 2),
  p05      = round(c(quantile(climate_df$tmax_mean, 0.05, na.rm = TRUE),
                      quantile(climate_df$tmax_max, 0.05, na.rm = TRUE),
                      quantile(climate_df$tmin_mean, 0.05, na.rm = TRUE),
                      quantile(climate_df$precip_total, 0.05, na.rm = TRUE),
                      quantile(climate_df$days_above_p90, 0.05, na.rm = TRUE),
                      quantile(climate_df$days_above_p95, 0.05, na.rm = TRUE),
                      quantile(climate_df$days_above_p99, 0.05, na.rm = TRUE),
                      quantile(climate_df$heatwave_days, 0.05, na.rm = TRUE)), 2),
  median   = round(c(median(climate_df$tmax_mean, na.rm = TRUE),
                      median(climate_df$tmax_max, na.rm = TRUE),
                      median(climate_df$tmin_mean, na.rm = TRUE),
                      median(climate_df$precip_total, na.rm = TRUE),
                      median(climate_df$days_above_p90, na.rm = TRUE),
                      median(climate_df$days_above_p95, na.rm = TRUE),
                      median(climate_df$days_above_p99, na.rm = TRUE),
                      median(climate_df$heatwave_days, na.rm = TRUE)), 2),
  p95      = round(c(quantile(climate_df$tmax_mean, 0.95, na.rm = TRUE),
                      quantile(climate_df$tmax_max, 0.95, na.rm = TRUE),
                      quantile(climate_df$tmin_mean, 0.95, na.rm = TRUE),
                      quantile(climate_df$precip_total, 0.95, na.rm = TRUE),
                      quantile(climate_df$days_above_p90, 0.95, na.rm = TRUE),
                      quantile(climate_df$days_above_p95, 0.95, na.rm = TRUE),
                      quantile(climate_df$days_above_p99, 0.95, na.rm = TRUE),
                      quantile(climate_df$heatwave_days, 0.95, na.rm = TRUE)), 2),
  stringsAsFactors = FALSE
)

print(as.data.frame(climate_stats), row.names = FALSE)


# ==============================================================================
# 5. SOCIOECONOMIC AND REMOTENESS DISTRIBUTIONS
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("D. SA4 Characteristics\n")
cat("=" |> strrep(50), "\n")

sa4_chars <- list()

if (!is.null(seifa)) {
  seifa_dist <- seifa |>
    group_by(irsd_quintile) |>
    summarise(
      n_sa4 = n(),
      mean_irsd = round(mean(irsd_score, na.rm = TRUE), 1),
      total_pop = sum(population, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(pct_pop = round(100 * total_pop / sum(total_pop), 1))

  cat("\nSEIFA IRSD quintile distribution:\n")
  print(as.data.frame(seifa_dist), row.names = FALSE)
  sa4_chars[["seifa"]] <- seifa_dist
}

if (!is.null(remote)) {
  remote_dist <- remote |>
    group_by(remoteness) |>
    summarise(n_sa4 = n(), .groups = "drop") |>
    mutate(pct = round(100 * n_sa4 / sum(n_sa4), 1))

  cat("\nRemoteness distribution:\n")
  print(as.data.frame(remote_dist), row.names = FALSE)
  sa4_chars[["remoteness"]] <- remote_dist
}

# Jurisdiction distribution
juris_dist <- panel |>
  filter(analysis_group == "cvd_total") |>
  distinct(sa4_code, jurisdiction) |>
  group_by(jurisdiction) |>
  summarise(n_sa4 = n(), .groups = "drop") |>
  mutate(pct = round(100 * n_sa4 / sum(n_sa4), 1)) |>
  arrange(desc(n_sa4))

cat("\nJurisdiction distribution:\n")
print(as.data.frame(juris_dist), row.names = FALSE)


# ==============================================================================
# 6. TEMPORAL PATTERNS
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("E. Temporal Patterns\n")
cat("=" |> strrep(50), "\n")

yearly_stats <- panel |>
  filter(analysis_group == "resp_total",
         suppressed != "True", !is.na(count)) |>
  group_by(year) |>
  summarise(
    n_obs     = n(),
    total_rx  = sum(count),
    mean_count = round(mean(count), 1),
    mean_tmax = round(mean(tmax_mean, na.rm = TRUE), 1),
    .groups = "drop"
  )

cat("\nYearly summary (Respiratory):\n")
print(as.data.frame(yearly_stats), row.names = FALSE)


# ==============================================================================
# 7. PM2.5 SUMMARY (BUSHFIRE PERIOD)
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("F. PM2.5 / Bushfire Period\n")
cat("=" |> strrep(50), "\n")

if ("pm25_mean" %in% names(panel)) {
  pm25_avail <- panel |>
    filter(!is.na(pm25_mean), analysis_group == "cvd_total")

  cat("  PM2.5 coverage:", nrow(pm25_avail), "SA4-weeks\n")

  if (nrow(pm25_avail) > 0) {
    # Overall
    cat("  Overall PM2.5: mean =", round(mean(pm25_avail$pm25_mean, na.rm = TRUE), 1),
        ", median =", round(median(pm25_avail$pm25_mean, na.rm = TRUE), 1),
        ", p95 =", round(quantile(pm25_avail$pm25_mean, 0.95, na.rm = TRUE), 1), "\n")

    # Black Summer period
    bs_start <- as.Date("2019-10-01")
    bs_end   <- as.Date("2020-02-28")
    bs <- pm25_avail |> filter(week_start >= bs_start, week_start <= bs_end)
    if (nrow(bs) > 0) {
      cat("  Black Summer PM2.5: mean =", round(mean(bs$pm25_mean, na.rm = TRUE), 1),
          ", max =", round(max(bs$pm25_max, na.rm = TRUE), 1), "\n")
    }
  }
} else {
  cat("  No PM2.5 column in panel.\n")
}


# ==============================================================================
# 8. SUPPRESSION SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(50), "\n")
cat("G. Data Suppression\n")
cat("=" |> strrep(50), "\n")

supp_stats <- panel |>
  group_by(analysis_group) |>
  summarise(
    total = n(),
    suppressed = sum(suppressed == "True", na.rm = TRUE),
    pct_suppressed = round(100 * suppressed / total, 1),
    .groups = "drop"
  ) |>
  mutate(label = group_labels[analysis_group]) |>
  filter(!is.na(label)) |>
  arrange(match(analysis_group, names(group_labels)))

print(as.data.frame(
  supp_stats |> select(label, total, suppressed, pct_suppressed)
), row.names = FALSE)


# ==============================================================================
# 9. SAVE TABLES
# ==============================================================================

# Table 1: combined
table1 <- list(
  study_period = data.frame(
    metric = c("Study period", "Duration (years)", "Total weeks",
               "SA4 regions", "Panel rows"),
    value = c(
      paste(min(panel$week_start), "to", max(panel$week_start)),
      round(as.numeric(max(panel$week_start) - min(panel$week_start)) / 365.25, 1),
      n_distinct(panel$week_start),
      n_distinct(panel$sa4_code),
      format(nrow(panel), big.mark = ",")
    ),
    stringsAsFactors = FALSE
  )
)

write.csv(outcome_stats,
          file.path(tab_dir, "descriptive_by_group.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/descriptive_by_group.csv\n")

write.csv(climate_stats,
          file.path(tab_dir, "descriptive_climate.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/descriptive_climate.csv\n")

# Combined Table 1
t1_rows <- data.frame(
  section  = c(rep("Study", 4),
               rep("Outcome", nrow(outcome_stats)),
               rep("Climate", nrow(climate_stats))),
  variable = c("Study period", "SA4 regions", "Total weeks", "Jurisdictions",
               outcome_stats$label,
               climate_stats$variable),
  summary  = c(
    paste(min(panel$week_start), "to", max(panel$week_start)),
    as.character(n_distinct(panel$sa4_code)),
    as.character(n_distinct(panel$week_start)),
    as.character(n_distinct(panel$jurisdiction)),
    paste0(format(outcome_stats$total_rx, big.mark = ","),
           " (mean ", outcome_stats$mean_count,
           ", SD ", outcome_stats$sd_count, ")"),
    paste0(climate_stats$mean, " (", climate_stats$sd, ")")
  ),
  stringsAsFactors = FALSE
)

write.csv(t1_rows,
          file.path(tab_dir, "table1_descriptive.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/table1_descriptive.csv\n")


# ==============================================================================
# 10. FIGURES
# ==============================================================================

cat("\nGenerating descriptive figures...\n")

# --- Outcome distribution by group ---
plot_df <- panel |>
  filter(suppressed != "True", !is.na(count),
         analysis_group %in% names(group_labels)) |>
  mutate(label = group_labels[analysis_group])

p_outcome <- ggplot(plot_df, aes(x = count)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, colour = "white") +
  facet_wrap(~ label, scales = "free", ncol = 2) +
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = "Distribution of weekly prescription counts by medication group",
    x = "Weekly count per SA4",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "descriptive_outcome_distributions.png"),
       p_outcome, width = 12, height = 10, dpi = 300)
cat("  -> Saved: descriptive_outcome_distributions.png\n")

# --- Temperature distribution ---
p_temp <- ggplot(climate_df, aes(x = tmax_mean)) +
  geom_histogram(bins = 60, fill = "#d73027", alpha = 0.6, colour = "white") +
  geom_vline(xintercept = quantile(climate_df$tmax_mean, c(0.05, 0.50, 0.95), na.rm = TRUE),
             linetype = c("dashed", "solid", "dashed"),
             colour = c("blue", "black", "red"), linewidth = 0.7) +
  annotate("text", x = quantile(climate_df$tmax_mean, 0.05, na.rm = TRUE), y = Inf,
           label = "p05", vjust = 2, hjust = 1.2, size = 3.5, colour = "blue") +
  annotate("text", x = quantile(climate_df$tmax_mean, 0.50, na.rm = TRUE), y = Inf,
           label = "Median", vjust = 2, hjust = -0.2, size = 3.5) +
  annotate("text", x = quantile(climate_df$tmax_mean, 0.95, na.rm = TRUE), y = Inf,
           label = "p95", vjust = 2, hjust = -0.2, size = 3.5, colour = "red") +
  labs(
    title = "Distribution of weekly mean maximum temperature across SA4 regions",
    x = "Weekly mean Tmax (\u00b0C)",
    y = "Frequency (SA4-weeks)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "descriptive_temperature_distribution.png"),
       p_temp, width = 9, height = 5, dpi = 300)
cat("  -> Saved: descriptive_temperature_distribution.png\n")

# --- Time series of mean prescribing ---
ts_df <- panel |>
  filter(suppressed != "True", !is.na(count),
         analysis_group %in% c("cvd_total", "mh_total", "resp_total")) |>
  mutate(label = group_labels[analysis_group]) |>
  group_by(label, week_start) |>
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

p_ts <- ggplot(ts_df, aes(x = week_start, y = mean_count, colour = label)) +
  geom_line(alpha = 0.5, linewidth = 0.3) +
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 0.8) +
  facet_wrap(~ label, scales = "free_y", ncol = 1) +
  labs(
    title = "Weekly mean prescription counts over study period",
    x = NULL, y = "Mean weekly count per SA4"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

ggsave(file.path(fig_dir, "descriptive_time_series.png"),
       p_ts, width = 12, height = 8, dpi = 300)
cat("  -> Saved: descriptive_time_series.png\n")


# ==============================================================================
# 11. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 28 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Outputs:\n")
cat("  outputs/tables/table1_descriptive.csv\n")
cat("  outputs/tables/descriptive_by_group.csv\n")
cat("  outputs/tables/descriptive_climate.csv\n")
cat("  outputs/figures/descriptive_outcome_distributions.png\n")
cat("  outputs/figures/descriptive_temperature_distribution.png\n")
cat("  outputs/figures/descriptive_time_series.png\n")
