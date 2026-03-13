# ==============================================================================
# Script 10: Attributable Burden Estimation
# ==============================================================================
# Computes the fraction and number of respiratory prescriptions attributable
# to heat exposure, using the attrdl() function from the dlnm package.
#
# Also estimates attributable burden by SEIFA quintile to quantify the equity
# dimension in absolute terms (not just relative risks).
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv
#
# Outputs:
#   outputs/tables/attributable_burden_overall.csv
#   outputs/tables/attributable_burden_by_seifa.csv
#   outputs/figures/attributable_burden_bar.png
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

# Source the attrdl() function (provided by dlnm authors, not in the CRAN package)
# See: https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata
attrdl_url <- "https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/attrdl.R"
attrdl_file <- file.path(tempdir(), "attrdl.R")
if (!exists("attrdl", mode = "function")) {
  cat("Downloading attrdl() function...\n")
  download.file(attrdl_url, attrdl_file, quiet = TRUE)
  source(attrdl_file)
  cat("  -> attrdl() loaded.\n")
}

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

# Focus on medication groups with significant heat effects
outcome_groups <- c(
  "resp_total",
  "resp_relievers",
  "resp_preventers"
)

group_labels <- c(
  resp_total      = "Respiratory (total)",
  resp_relievers  = "Resp. Relievers",
  resp_preventers = "Resp. Preventers"
)

cat("=" |> strrep(70), "\n")
cat("Script 10: Attributable Burden Estimation\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("Loading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

# Load SEIFA
seifa <- read.csv(file.path(ref_dir, "seifa_sa4_quintiles.csv"),
                  stringsAsFactors = FALSE)
seifa$sa4_code <- as.character(seifa$sa4_code)

panel <- panel |>
  left_join(seifa |> select(sa4_code, irsd_quintile), by = "sa4_code")

# Time variables
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
# 2. FIT DLNM AND COMPUTE ATTRIBUTABLE BURDEN (OVERALL)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Computing attributable burden (overall)\n")
cat("=" |> strrep(70), "\n")

compute_attributable <- function(data, group_name, label_suffix = "") {
  df <- data |>
    filter(analysis_group == group_name,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 500) {
    cat("  Too few observations for", group_name, label_suffix, "\n")
    return(NULL)
  }

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

  temp_ref <- median(df$tmax_mean, na.rm = TRUE)

  # --- Attributable burden: heat (above median) ---
  # attrdl computes attributable number for the entire series
  cat("  Computing attributable numbers for heat (above median)...\n")

  # Forward perspective: burden attributable to temperatures above the reference
  an_heat <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(temp_ref, max(df$tmax_mean, na.rm = TRUE)),
    type = "an",
    tot  = TRUE
  )

  af_heat <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(temp_ref, max(df$tmax_mean, na.rm = TRUE)),
    type = "af",
    tot  = TRUE
  )

  # Also compute for extreme heat only (above 90th percentile)
  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)

  an_extreme <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(temp_p90, max(df$tmax_mean, na.rm = TRUE)),
    type = "an",
    tot  = TRUE
  )

  af_extreme <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(temp_p90, max(df$tmax_mean, na.rm = TRUE)),
    type = "af",
    tot  = TRUE
  )

  # Cold burden (below median)
  an_cold <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(min(df$tmax_mean, na.rm = TRUE), temp_ref),
    type = "an",
    tot  = TRUE
  )

  af_cold <- attrdl(
    x    = df$tmax_mean,
    basis = cb,
    cases = df$count,
    model = model,
    cen  = temp_ref,
    range = c(min(df$tmax_mean, na.rm = TRUE), temp_ref),
    type = "af",
    tot  = TRUE
  )

  total_rx <- sum(df$count, na.rm = TRUE)

  result <- data.frame(
    group         = group_name,
    label         = group_labels[group_name],
    suffix        = label_suffix,
    total_rx      = total_rx,
    temp_ref      = round(temp_ref, 1),
    # Heat (above median)
    an_heat       = round(an_heat),
    af_heat_pct   = round(af_heat * 100, 2),
    # Extreme heat (above p90)
    an_extreme    = round(an_extreme),
    af_extreme_pct = round(af_extreme * 100, 2),
    # Cold (below median)
    an_cold       = round(an_cold),
    af_cold_pct   = round(af_cold * 100, 2),
    n_obs         = nrow(df),
    n_sa4         = n_distinct(df$sa4_code),
    stringsAsFactors = FALSE
  )

  cat("  ", group_labels[group_name], label_suffix, ":\n")
  cat("    Total prescriptions:", format(total_rx, big.mark = ","), "\n")
  cat("    Heat (>median):  AN =", format(round(an_heat), big.mark = ","),
      " AF =", round(af_heat * 100, 2), "%\n")
  cat("    Extreme (>p90):  AN =", format(round(an_extreme), big.mark = ","),
      " AF =", round(af_extreme * 100, 2), "%\n")
  cat("    Cold (<median):  AN =", format(round(an_cold), big.mark = ","),
      " AF =", round(af_cold * 100, 2), "%\n")

  return(result)
}


# Overall attributable burden
overall_results <- list()
for (grp in outcome_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  res <- compute_attributable(panel, grp)
  if (!is.null(res)) overall_results[[grp]] <- res
}

overall_table <- bind_rows(overall_results)
write.csv(overall_table,
          file.path(tab_dir, "attributable_burden_overall.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/attributable_burden_overall.csv\n")


# ==============================================================================
# 3. ATTRIBUTABLE BURDEN BY SEIFA QUINTILE
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Computing attributable burden by SEIFA quintile (respiratory only)\n")
cat("=" |> strrep(70), "\n")

seifa_results <- list()

for (q in 1:5) {
  cat("\n  Quintile ", q, ":\n")
  panel_q <- panel |> filter(irsd_quintile == q)

  for (grp in outcome_groups) {
    res <- compute_attributable(panel_q, grp,
                                label_suffix = paste0("Q", q))
    if (!is.null(res)) {
      res$quintile <- q
      seifa_results[[paste0(grp, "_q", q)]] <- res
    }
  }
}

seifa_table <- bind_rows(seifa_results)
write.csv(seifa_table,
          file.path(tab_dir, "attributable_burden_by_seifa.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/attributable_burden_by_seifa.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating plots...\n")

# --- Bar chart: AF by medication group (overall) ---
p_overall <- ggplot(overall_table,
                     aes(x = label, y = af_heat_pct, fill = label)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(af_heat_pct, "%")),
            vjust = ifelse(overall_table$af_heat_pct >= 0, -0.5, 1.5),
            size = 4) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(
    title = "Attributable fraction: heat exposure (above median temperature)",
    subtitle = "Percentage of total prescriptions attributable to heat",
    x = NULL,
    y = "Attributable fraction (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13))

# --- Bar chart: AF by SEIFA quintile (resp_total) ---
seifa_resp <- seifa_table |> filter(group == "resp_total")

if (nrow(seifa_resp) > 0) {
  quintile_labels <- c("Q1\n(Most\ndisadvantaged)", "Q2", "Q3", "Q4",
                       "Q5\n(Least\ndisadvantaged)")
  quintile_colours <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")

  p_seifa <- ggplot(seifa_resp,
                     aes(x = factor(quintile), y = af_heat_pct,
                         fill = factor(quintile))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = paste0(af_heat_pct, "%")),
              vjust = ifelse(seifa_resp$af_heat_pct >= 0, -0.5, 1.5),
              size = 4) +
    geom_hline(yintercept = 0, colour = "grey40") +
    scale_fill_manual(values = quintile_colours, guide = "none") +
    scale_x_discrete(labels = quintile_labels) +
    labs(
      title = "Heat-attributable fraction for respiratory prescriptions by SEIFA quintile",
      x = "SEIFA IRSD Quintile",
      y = "Attributable fraction (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  combined <- p_overall / p_seifa + plot_layout(heights = c(1, 1))
  ggsave(file.path(fig_dir, "attributable_burden_bar.png"),
         combined, width = 9, height = 10, dpi = 300)
  cat("  -> Saved: attributable_burden_bar.png\n")
} else {
  ggsave(file.path(fig_dir, "attributable_burden_bar.png"),
         p_overall, width = 8, height = 5, dpi = 300)
  cat("  -> Saved: attributable_burden_bar.png\n")
}


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 10 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Overall attributable burden:\n")
overall_table |>
  select(label, total_rx, an_heat, af_heat_pct, an_extreme, af_extreme_pct) |>
  print(row.names = FALSE)

if (nrow(seifa_resp) > 0) {
  cat("\nRespiratory AF by SEIFA quintile:\n")
  seifa_resp |>
    select(quintile, total_rx, an_heat, af_heat_pct, an_extreme, af_extreme_pct) |>
    print(row.names = FALSE)
}

cat("\nOutputs:\n")
cat("  outputs/tables/attributable_burden_overall.csv\n")
cat("  outputs/tables/attributable_burden_by_seifa.csv\n")
cat("  outputs/figures/attributable_burden_bar.png\n")
