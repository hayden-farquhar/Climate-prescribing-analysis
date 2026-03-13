# ==============================================================================
# Script 09: Remoteness Stratification
# ==============================================================================
# Stratifies DLNM heat-prescribing models by ABS Remoteness Area classification
# to test whether rural/remote areas show greater sensitivity to heat exposure.
#
# Approach:
#   1. Spatial overlay: SA4 polygons x Remoteness Areas -> classify each SA4
#   2. Fit DLNM separately for each remoteness category (3 primary groups)
#   3. Fit pooled model with heat x remoteness interaction (F-test)
#   4. Generate comparative exposure-response curves and forest plot
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/spatial/SA4_2021/SA4_2021_AUST_GDA2020.shp
#   data/spatial/RA_2021/RA_2021_AUST_GDA2020.shp
#
# Outputs:
#   data/reference/remoteness_sa4.csv
#   outputs/tables/remoteness_stratified_rr.csv
#   outputs/tables/remoteness_interaction_tests.csv
#   outputs/figures/remoteness_exposure_response_*.png
#   outputs/figures/remoteness_forest_plot.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "sf", "dplyr", "lubridate",
#                       "ggplot2", "patchwork", "tidyr"))
# ==============================================================================

library(dlnm)
library(splines)
library(sf)
sf_use_s2(FALSE)  # Disable s2 — ABS shapefiles have degenerate edges that s2 rejects
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)
library(tidyr)

# --- Project paths ---
project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
ref_dir     <- file.path(project_dir, "data", "reference")
spatial_dir <- file.path(project_dir, "data", "spatial")
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
cat("Script 09: Remoteness Stratification\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. SPATIAL OVERLAY: SA4 -> REMOTENESS CLASSIFICATION
# ==============================================================================

cat("Classifying SA4 regions by remoteness...\n")

# Load shapefiles
sa4_shp <- st_read(file.path(spatial_dir, "SA4_2021", "SA4_2021_AUST_GDA2020.shp"),
                   quiet = TRUE)
ra_shp  <- st_read(file.path(spatial_dir, "RA_2021", "RA_2021_AUST_GDA2020.shp"),
                   quiet = TRUE)

# Ensure matching CRS
ra_shp <- st_transform(ra_shp, st_crs(sa4_shp))

cat("  SA4 polygons:", nrow(sa4_shp), "\n")
cat("  RA polygons:", nrow(ra_shp), "\n")
cat("  RA categories:", paste(sort(unique(ra_shp$RA_NAME21)), collapse = ", "), "\n")

# Area-weighted overlay: for each SA4, find the dominant remoteness category
# Compute intersection areas
sa4_simple <- sa4_shp |> select(SA4_CODE21, SA4_NAME21)
ra_simple  <- ra_shp |> select(RA_CODE21, RA_NAME21)

cat("  Computing spatial intersection (may take a moment)...\n")
intersection <- st_intersection(sa4_simple, ra_simple)
intersection$area <- as.numeric(st_area(intersection))

# For each SA4, find the RA category covering the largest area
sa4_remoteness <- intersection |>
  st_drop_geometry() |>
  group_by(SA4_CODE21, SA4_NAME21) |>
  slice_max(area, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(sa4_code = SA4_CODE21, sa4_name = SA4_NAME21,
         ra_code = RA_CODE21, ra_name = RA_NAME21, ra_area = area)

# Simplify remoteness into 3 categories for adequate sample size:
# Major Cities, Inner Regional, Outer Regional + Remote + Very Remote
sa4_remoteness <- sa4_remoteness |>
  mutate(
    remoteness = case_when(
      grepl("Major Cities", ra_name) ~ "Major Cities",
      grepl("Inner Regional", ra_name) ~ "Inner Regional",
      TRUE ~ "Outer Regional/Remote"
    ),
    remoteness = factor(remoteness,
                        levels = c("Major Cities", "Inner Regional", "Outer Regional/Remote"))
  )

cat("\n  SA4 remoteness classification:\n")
print(table(sa4_remoteness$remoteness))

# Save
write.csv(sa4_remoteness,
          file.path(ref_dir, "remoteness_sa4.csv"),
          row.names = FALSE)
cat("  -> Saved: data/reference/remoteness_sa4.csv\n")


# ==============================================================================
# 2. LOAD PANEL AND MERGE REMOTENESS
# ==============================================================================

cat("\nLoading panel data...\n")
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

panel <- panel |>
  left_join(sa4_remoteness |> select(sa4_code, remoteness),
            by = "sa4_code")

n_matched <- sum(!is.na(panel$remoteness))
cat("  Panel rows with remoteness match:", format(n_matched, big.mark = ","),
    "/", format(nrow(panel), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * n_matched / nrow(panel)))

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
# 3. STRATIFIED DLNM BY REMOTENESS CATEGORY
# ==============================================================================

remoteness_levels <- levels(sa4_remoteness$remoteness)

remoteness_colours <- c(
  "Major Cities"           = "#4575b4",
  "Inner Regional"         = "#fee090",
  "Outer Regional/Remote"  = "#d73027"
)

fit_dlnm_remoteness <- function(data, group_name, remote_cat) {
  df <- data |>
    filter(analysis_group == group_name,
           remoteness == remote_cat,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  if (nrow(df) < 200) {
    cat("    ", remote_cat, ": too few obs (", nrow(df), "), skipping\n")
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
    group      = group_name,
    remoteness = remote_cat,
    temp_ref   = round(temp_ref, 1),
    temp_p95   = round(temp_p95, 1),
    rr_p95     = round(pred_pctl$allRRfit[1], 4),
    rr_p95_lo  = round(pred_pctl$allRRlow[1], 4),
    rr_p95_hi  = round(pred_pctl$allRRhigh[1], 4),
    temp_p99   = round(temp_p99, 1),
    rr_p99     = round(pred_pctl$allRRfit[2], 4),
    rr_p99_lo  = round(pred_pctl$allRRlow[2], 4),
    rr_p99_hi  = round(pred_pctl$allRRhigh[2], 4),
    n_obs      = nrow(df),
    n_sa4      = n_distinct(df$sa4_code),
    dispersion = round(summary(model)$dispersion, 2),
    stringsAsFactors = FALSE
  )

  return(list(pred = pred, summary = rr_row, model = model,
              cb = cb, data = df, temp_ref = temp_ref))
}


cat("\n", "=" |> strrep(70), "\n")
cat("Fitting stratified DLNM models (3 groups x 3 remoteness = 9 models)\n")
cat("=" |> strrep(70), "\n")

all_results   <- list()
all_summaries <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")
  for (rc in remoteness_levels) {
    cat("  ", rc, "... ", sep = "")
    res <- fit_dlnm_remoteness(panel, grp, rc)
    if (!is.null(res)) {
      key <- paste0(grp, "_", gsub("[/ ]", "_", rc))
      all_results[[key]] <- res
      all_summaries[[key]] <- res$summary
      cat("RR@p95 =", res$summary$rr_p95,
          "(", res$summary$rr_p95_lo, "-", res$summary$rr_p95_hi, ")",
          "n =", format(res$summary$n_obs, big.mark = ","), "\n")
    }
  }
}


# ==============================================================================
# 4. INTERACTION TEST (pooled model with remoteness x temperature)
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Testing heat x remoteness interaction\n")
cat("=" |> strrep(70), "\n")

interaction_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(remoteness)) |>
    mutate(remote_f = factor(remoteness,
                             levels = c("Major Cities", "Inner Regional",
                                        "Outer Regional/Remote"))) |>
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

  cat("  Fitting null model (no interaction)...\n")
  mod_null <- glm(
    count ~ cb + remote_f +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  cat("  Fitting interaction model (tmax_mean x remoteness)...\n")
  mod_int <- glm(
    count ~ cb + remote_f + tmax_mean:remote_f +
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
# 5. SAVE RESULTS TABLES
# ==============================================================================

strat_table <- bind_rows(all_summaries)
strat_table$label <- group_labels[strat_table$group]

write.csv(strat_table,
          file.path(tab_dir, "remoteness_stratified_rr.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/remoteness_stratified_rr.csv\n")

int_table <- bind_rows(interaction_results)
write.csv(int_table,
          file.path(tab_dir, "remoteness_interaction_tests.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/remoteness_interaction_tests.csv\n")


# ==============================================================================
# 6. VISUALISATION
# ==============================================================================

cat("\nGenerating remoteness plots...\n")

# --- Exposure-response by remoteness ---
plot_remoteness_er <- function(grp) {
  lbl <- group_labels[grp]

  plot_data <- list()
  for (rc in remoteness_levels) {
    key <- paste0(grp, "_", gsub("[/ ]", "_", rc))
    if (!key %in% names(all_results)) next
    res <- all_results[[key]]
    pred <- res$pred

    df_q <- data.frame(
      temp       = as.numeric(names(pred$allRRfit)),
      rr         = pred$allRRfit,
      rr_lo      = pred$allRRlow,
      rr_hi      = pred$allRRhigh,
      remoteness = rc
    )
    plot_data[[rc]] <- df_q
  }

  if (length(plot_data) == 0) return(NULL)
  df_plot <- bind_rows(plot_data)

  p <- ggplot(df_plot, aes(x = temp, y = rr, colour = remoteness, fill = remoteness)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = remoteness_colours, name = "Remoteness") +
    scale_fill_manual(values = remoteness_colours, name = "Remoteness") +
    labs(
      title = paste0(lbl, " — Exposure-response by remoteness"),
      x = "Weekly mean Tmax (°C)",
      y = "Cumulative RR (vs median)"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))

  return(p)
}

for (grp in primary_groups) {
  p <- plot_remoteness_er(grp)
  if (!is.null(p)) {
    ggsave(file.path(fig_dir, paste0("remoteness_exposure_response_", grp, ".png")),
           p, width = 8, height = 6, dpi = 300)
  }
}
cat("  -> Saved: remoteness_exposure_response_*.png\n")

# Combined panel
plots <- lapply(primary_groups, plot_remoteness_er)
plots <- plots[!sapply(plots, is.null)]
if (length(plots) >= 2) {
  combined <- wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "remoteness_exposure_response_panel.png"),
         combined, width = 9, height = 14, dpi = 300)
  cat("  -> Saved: remoteness_exposure_response_panel.png\n")
}

# --- Forest plot ---
cat("Generating forest plot...\n")
forest_data <- strat_table |>
  mutate(
    remoteness_f = factor(remoteness,
                          levels = rev(c("Major Cities", "Inner Regional",
                                         "Outer Regional/Remote"))),
    label_f = factor(label, levels = c("Cardiovascular", "Mental Health", "Respiratory"))
  )

p_forest <- ggplot(forest_data,
                    aes(x = rr_p95, y = remoteness_f, colour = remoteness)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = remoteness_colours, guide = "none") +
  facet_wrap(~ label_f, scales = "free_x", ncol = 3) +
  labs(
    title = "Cumulative RR at 95th percentile temperature by remoteness",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "remoteness_forest_plot.png"),
       p_forest, width = 12, height = 5, dpi = 300)
cat("  -> Saved: remoteness_forest_plot.png\n")


# ==============================================================================
# 7. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 09 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Stratified RRs at 95th percentile:\n")
print(
  strat_table |>
    select(label, remoteness, rr_p95, rr_p95_lo, rr_p95_hi, n_obs, n_sa4) |>
    as.data.frame(),
  row.names = FALSE
)

cat("\nInteraction tests (heat x remoteness):\n")
print(
  int_table |>
    select(label, f_stat, p_value, significant) |>
    as.data.frame(),
  row.names = FALSE
)

cat("\nOutputs:\n")
cat("  data/reference/remoteness_sa4.csv\n")
cat("  outputs/tables/remoteness_stratified_rr.csv\n")
cat("  outputs/tables/remoteness_interaction_tests.csv\n")
cat("  outputs/figures/remoteness_exposure_response_*.png\n")
cat("  outputs/figures/remoteness_forest_plot.png\n")
