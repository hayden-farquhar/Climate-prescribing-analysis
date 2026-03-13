# ==============================================================================
# Script 16: Leave-One-Out Jurisdiction Sensitivity Analysis
# ==============================================================================
# Tests whether the primary DLNM results are driven by any single state or
# territory by re-fitting models after dropping each jurisdiction one at a time.
#
# Also includes the SEIFA x remoteness cross-tabulation analysis.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/reference/seifa_sa4_quintiles.csv
#   data/reference/remoteness_sa4.csv  (from script 09)
#
# Outputs:
#   outputs/tables/loo_jurisdiction_rr.csv
#   outputs/tables/seifa_remoteness_crosstab.csv
#   outputs/figures/loo_jurisdiction_forest.png
#   outputs/figures/seifa_remoteness_heatmap.png
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

primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

cat("=" |> strrep(70), "\n")
cat("Script 16: Leave-One-Out Jurisdiction Sensitivity\n")
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

jurisdictions <- sort(unique(panel$jurisdiction))
cat("  Jurisdictions:", paste(jurisdictions, collapse = ", "), "\n")
cat("  Total:", length(jurisdictions), "states/territories\n")


# ==============================================================================
# 2. LEAVE-ONE-OUT DLNM
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Leave-one-out analysis:", length(jurisdictions), "jurisdictions x",
    length(primary_groups), "groups =",
    length(jurisdictions) * length(primary_groups), "models\n")
cat("=" |> strrep(70), "\n")

fit_dlnm_loo <- function(data, group_name, drop_jur) {
  df <- data |>
    filter(analysis_group == group_name,
           jurisdiction != drop_jur,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 500) return(NULL)

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
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  pred_pctl <- crosspred(cb, model, at = c(temp_p95, temp_p99), cen = temp_ref)

  data.frame(
    group       = group_name,
    dropped     = drop_jur,
    temp_ref    = round(temp_ref, 1),
    temp_p95    = round(temp_p95, 1),
    rr_p95      = round(pred_pctl$allRRfit[1], 4),
    rr_p95_lo   = round(pred_pctl$allRRlow[1], 4),
    rr_p95_hi   = round(pred_pctl$allRRhigh[1], 4),
    temp_p99    = round(temp_p99, 1),
    rr_p99      = round(pred_pctl$allRRfit[2], 4),
    rr_p99_lo   = round(pred_pctl$allRRlow[2], 4),
    rr_p99_hi   = round(pred_pctl$allRRhigh[2], 4),
    n_obs       = nrow(df),
    n_sa4       = n_distinct(df$sa4_code),
    stringsAsFactors = FALSE
  )
}

loo_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  for (jur in jurisdictions) {
    cat("  Dropping ", jur, "... ", sep = "")
    res <- fit_dlnm_loo(panel, grp, jur)
    if (!is.null(res)) {
      loo_results[[paste0(grp, "_", jur)]] <- res
      cat("RR@p95 =", res$rr_p95,
          "(", res$rr_p95_lo, "-", res$rr_p95_hi, ")\n")
    } else {
      cat("skipped\n")
    }
  }
}

loo_table <- bind_rows(loo_results)
loo_table$label <- group_labels[loo_table$group]

write.csv(loo_table,
          file.path(tab_dir, "loo_jurisdiction_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/loo_jurisdiction_rr.csv\n")


# ==============================================================================
# 3. VISUALISATION: LOO FOREST PLOT
# ==============================================================================

cat("\nGenerating leave-one-out forest plots...\n")

# Load full-sample results for reference line
full_file <- file.path(tab_dir, "dlnm_results_summary.csv")
full_rr <- NULL
if (file.exists(full_file)) {
  full_results <- read.csv(full_file, stringsAsFactors = FALSE)
  full_rr <- full_results |>
    filter(group %in% primary_groups) |>
    select(group, full_rr = rr_p95)
}

for (grp in primary_groups) {
  df <- loo_table |> filter(group == grp)

  p <- ggplot(df, aes(x = rr_p95, y = reorder(dropped, rr_p95))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5, colour = "steelblue") +
    geom_point(colour = "steelblue", size = 2.5)

  # Add full-sample reference
  if (!is.null(full_rr)) {
    ref_val <- full_rr$full_rr[full_rr$group == grp]
    if (length(ref_val) == 1) {
      p <- p + geom_vline(xintercept = ref_val, linetype = "dotted",
                           colour = "#d73027", linewidth = 0.7) +
        annotate("text", x = ref_val, y = 0.5, label = "Full sample",
                 colour = "#d73027", size = 3, hjust = -0.1)
    }
  }

  p <- p +
    labs(
      title = paste0(group_labels[grp], " — Leave-one-out sensitivity"),
      subtitle = "Each row drops one state/territory; red line = full-sample estimate",
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = "Dropped jurisdiction"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13))

  ggsave(file.path(fig_dir, paste0("loo_jurisdiction_", grp, ".png")),
         p, width = 8, height = 5, dpi = 300)
}
cat("  -> Saved: loo_jurisdiction_*.png\n")

# Combined panel
plots <- list()
for (grp in primary_groups) {
  df <- loo_table |> filter(group == grp)
  p <- ggplot(df, aes(x = rr_p95, y = reorder(dropped, rr_p95))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5, colour = "steelblue") +
    geom_point(colour = "steelblue", size = 2) +
    labs(title = group_labels[grp], x = "Cumulative RR (95% CI)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
  plots[[grp]] <- p
}

combined <- wrap_plots(plots, ncol = 3)
ggsave(file.path(fig_dir, "loo_jurisdiction_panel.png"),
       combined, width = 16, height = 5, dpi = 300)
cat("  -> Saved: loo_jurisdiction_panel.png\n")


# ==============================================================================
# 4. SEIFA x REMOTENESS CROSS-TABULATION
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("SEIFA x Remoteness cross-tabulation\n")
cat("=" |> strrep(70), "\n")

# Load SEIFA and remoteness
seifa <- read.csv(file.path(ref_dir, "seifa_sa4_quintiles.csv"),
                  stringsAsFactors = FALSE)
seifa$sa4_code <- as.character(seifa$sa4_code)

remoteness_file <- file.path(ref_dir, "remoteness_sa4.csv")
if (!file.exists(remoteness_file)) {
  cat("  remoteness_sa4.csv not found — run script 09 first.\n")
  cat("  Skipping cross-tabulation.\n")
} else {
  remote <- read.csv(remoteness_file, stringsAsFactors = FALSE)
  remote$sa4_code <- as.character(remote$sa4_code)

  # Merge
  crosstab_data <- seifa |>
    select(sa4_code, irsd_quintile) |>
    left_join(remote |> select(sa4_code, remoteness), by = "sa4_code") |>
    filter(!is.na(remoteness))

  cat("  SA4s with both classifications:", nrow(crosstab_data), "\n")

  # Cross-tabulation
  ct <- table(
    SEIFA = paste0("Q", crosstab_data$irsd_quintile),
    Remoteness = crosstab_data$remoteness
  )
  cat("\n  Cross-tabulation (SA4 count):\n")
  print(ct)

  # Save
  ct_df <- as.data.frame(ct)
  names(ct_df) <- c("seifa_quintile", "remoteness", "n_sa4")
  write.csv(ct_df,
            file.path(tab_dir, "seifa_remoteness_crosstab.csv"),
            row.names = FALSE)
  cat("\n  -> Saved: outputs/tables/seifa_remoteness_crosstab.csv\n")

  # Heatmap
  p_heatmap <- ggplot(ct_df, aes(x = remoteness, y = seifa_quintile, fill = n_sa4)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = n_sa4), colour = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "#d73027", name = "n SA4s") +
    labs(
      title = "Distribution of SA4 regions: SEIFA x Remoteness",
      x = "Remoteness classification",
      y = "SEIFA IRSD Quintile"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )

  ggsave(file.path(fig_dir, "seifa_remoteness_heatmap.png"),
         p_heatmap, width = 8, height = 5, dpi = 300)
  cat("  -> Saved: seifa_remoteness_heatmap.png\n")

  # --- Cross-stratified analysis (if cell sizes allow) ---
  # Only for respiratory (strongest effect) — combine where needed
  cat("\n  Attempting SEIFA x remoteness stratified analysis (respiratory)...\n")

  panel_merged <- panel |>
    left_join(seifa |> select(sa4_code, irsd_quintile), by = "sa4_code") |>
    left_join(remote |> select(sa4_code, remoteness), by = "sa4_code") |>
    filter(!is.na(irsd_quintile), !is.na(remoteness))

  # Simplify SEIFA to 2 groups for adequate cell sizes
  panel_merged <- panel_merged |>
    mutate(
      seifa_group = ifelse(irsd_quintile <= 2, "Disadvantaged (Q1-2)",
                           "Advantaged (Q3-5)"),
      stratum = paste0(seifa_group, " / ", remoteness)
    )

  strata <- unique(panel_merged$stratum)
  cat("  Strata:", paste(strata, collapse = "; "), "\n")

  cross_results <- list()
  for (s in strata) {
    df_s <- panel_merged |>
      filter(analysis_group == "resp_total",
             stratum == s,
             suppressed != "True",
             !is.na(count),
             !is.na(tmax_mean)) |>
      arrange(sa4_code, week_start)

    cat("    ", s, ": n =", nrow(df_s), ", SA4s =", n_distinct(df_s$sa4_code), "... ")

    if (nrow(df_s) < 300 || n_distinct(df_s$sa4_code) < 3) {
      cat("too few, skipping\n")
      next
    }

    temp_knots <- quantile(df_s$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
    cb <- crossbasis(
      df_s$tmax_mean,
      lag = max_lag,
      argvar = list(fun = "ns", knots = temp_knots),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    )

    n_years  <- n_distinct(year(df_s$week_start))
    trend_df <- max(2, round(n_years * 2))

    model <- tryCatch(
      glm(
        count ~ cb +
          ns(time_index, df = trend_df) +
          sin1 + cos1 + sin2 + cos2 +
          precip_total +
          factor(sa4_code),
        data = df_s,
        family = quasipoisson(link = "log")
      ),
      error = function(e) NULL
    )

    if (is.null(model)) {
      cat("model failed\n")
      next
    }

    temp_ref <- median(df_s$tmax_mean, na.rm = TRUE)
    temp_p95 <- quantile(df_s$tmax_mean, 0.95, na.rm = TRUE)
    pred_pctl <- crosspred(cb, model, at = temp_p95, cen = temp_ref)

    cross_results[[s]] <- data.frame(
      stratum   = s,
      rr_p95    = round(pred_pctl$allRRfit[1], 4),
      rr_p95_lo = round(pred_pctl$allRRlow[1], 4),
      rr_p95_hi = round(pred_pctl$allRRhigh[1], 4),
      n_obs     = nrow(df_s),
      n_sa4     = n_distinct(df_s$sa4_code),
      stringsAsFactors = FALSE
    )

    cat("RR@p95 =", round(pred_pctl$allRRfit[1], 4), "\n")
  }

  if (length(cross_results) > 0) {
    cross_table <- bind_rows(cross_results)
    write.csv(cross_table,
              file.path(tab_dir, "seifa_remoteness_stratified_rr.csv"),
              row.names = FALSE)
    cat("  -> Saved: seifa_remoteness_stratified_rr.csv\n")

    # Forest plot
    p_cross <- ggplot(cross_table,
                       aes(x = rr_p95, y = reorder(stratum, rr_p95))) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                     height = 0.3, linewidth = 0.6, colour = "#d73027") +
      geom_point(colour = "#d73027", size = 2.5) +
      labs(
        title = "Respiratory — heat effect by SEIFA x remoteness",
        x = "Cumulative RR at 95th percentile (95% CI)",
        y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold", size = 13))

    ggsave(file.path(fig_dir, "seifa_remoteness_forest.png"),
           p_cross, width = 9, height = 5, dpi = 300)
    cat("  -> Saved: seifa_remoteness_forest.png\n")
  }
}


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 16 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Leave-one-out sensitivity (RR at p95):\n")
print(
  loo_table |>
    select(label, dropped, rr_p95, rr_p95_lo, rr_p95_hi) |>
    as.data.frame(),
  row.names = FALSE
)

cat("\nInterpretation:\n")
cat("  If all RRs remain similar after dropping each jurisdiction,\n")
cat("  the results are not driven by any single state/territory.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/loo_jurisdiction_rr.csv\n")
cat("  outputs/figures/loo_jurisdiction_*.png\n")
if (file.exists(remoteness_file)) {
  cat("  outputs/tables/seifa_remoteness_crosstab.csv\n")
  cat("  outputs/tables/seifa_remoteness_stratified_rr.csv\n")
  cat("  outputs/figures/seifa_remoteness_heatmap.png\n")
  cat("  outputs/figures/seifa_remoteness_forest.png\n")
}
