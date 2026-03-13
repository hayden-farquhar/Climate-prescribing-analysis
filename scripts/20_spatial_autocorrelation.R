# ==============================================================================
# Script 20: Spatial Autocorrelation of Model Residuals
# ==============================================================================
# Tests for spatial autocorrelation in DLNM residuals using Moran's I.
# If significant spatial correlation exists, standard errors may be too narrow.
#
# Approach:
#   1. Fit primary DLNM models
#   2. Aggregate deviance residuals to SA4 level (mean over time)
#   3. Compute spatial weights matrix from SA4 polygon contiguity
#   4. Moran's I test on mean residuals
#   5. Also test at weekly level (sample of weeks) to check temporal variation
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/spatial/SA4_2021/SA4_2021_AUST_GDA2020.shp
#
# Outputs:
#   outputs/tables/spatial_autocorrelation_tests.csv
#   outputs/figures/spatial_moran_*.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "sf", "spdep", "dplyr",
#                       "lubridate", "ggplot2"))
# ==============================================================================

library(dlnm)
library(splines)
library(sf)
sf_use_s2(FALSE)
library(spdep)
library(dplyr)
library(lubridate)
library(ggplot2)
library(patchwork)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
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
cat("Script 20: Spatial Autocorrelation\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. BUILD SPATIAL WEIGHTS MATRIX
# ==============================================================================

cat("Building spatial weights matrix...\n")

sa4_shp <- st_read(file.path(spatial_dir, "SA4_2021", "SA4_2021_AUST_GDA2020.shp"),
                   quiet = TRUE)

# Load panel to get the SA4s actually in the analysis
panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

analysis_sa4s <- unique(panel$sa4_code)

# Filter shapefile to analysis SA4s
sa4_shp <- sa4_shp |>
  filter(SA4_CODE21 %in% analysis_sa4s)

cat("  SA4s in analysis:", length(analysis_sa4s), "\n")
cat("  SA4s in shapefile:", nrow(sa4_shp), "\n")

# Queen contiguity neighbours
nb <- poly2nb(sa4_shp, queen = TRUE)
# Some SA4s may be islands (no neighbours) — use k-nearest as fallback
n_islands <- sum(card(nb) == 0)
cat("  Island SA4s (no contiguous neighbours):", n_islands, "\n")

if (n_islands > 0) {
  # For island SA4s, add nearest neighbour
  coords <- st_centroid(sa4_shp) |> st_coordinates()
  knn <- knearneigh(coords, k = 1)
  nb_knn <- knn2nb(knn)
  # Combine: use contiguity where available, knn for islands
  for (i in which(card(nb) == 0)) {
    nb[[i]] <- nb_knn[[i]]
  }
  cat("  Added k-nearest neighbours for islands.\n")
}

listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
cat("  Spatial weights matrix: W-style, queen contiguity\n")


# ==============================================================================
# 2. FIT MODELS AND COMPUTE MORAN'S I
# ==============================================================================

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

all_results <- list()

for (grp in primary_groups) {
  cat("\n", "=" |> strrep(50), "\n")
  cat("Moran's I for:", group_labels[grp], "\n")
  cat("=" |> strrep(50), "\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  # Fit model
  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    df$tmax_mean,
    lag = max_lag,
    argvar = list(fun = "ns", knots = temp_knots),
    arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
  )

  n_years  <- n_distinct(df$year)
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

  # Subset df to rows actually used by glm (drops NA from crossbasis lags etc.)
  if (!is.null(model$na.action)) {
    df <- df[-model$na.action, ]
  }
  df$dev_resid <- as.numeric(residuals(model, type = "deviance"))

  # --- A. Mean residuals by SA4 (time-averaged) ---
  mean_resid <- df |>
    group_by(sa4_code) |>
    summarise(mean_resid = mean(dev_resid, na.rm = TRUE), .groups = "drop")

  # Match order to shapefile
  sa4_order <- match(sa4_shp$SA4_CODE21, mean_resid$sa4_code)
  resid_vec <- mean_resid$mean_resid[sa4_order]

  # Handle NAs (SA4s in shapefile but not in data)
  valid <- !is.na(resid_vec)

  if (sum(valid) >= 10) {
    # Subset listw to valid SA4s
    valid_idx <- which(valid)
    nb_sub <- subset(nb, valid)
    listw_sub <- nb2listw(nb_sub, style = "W", zero.policy = TRUE)
    resid_sub <- resid_vec[valid_idx]

    moran_overall <- moran.test(resid_sub, listw_sub, zero.policy = TRUE)
    cat("  Overall Moran's I:", round(moran_overall$estimate[1], 4),
        "  p =", format.pval(moran_overall$p.value, digits = 3), "\n")
  } else {
    moran_overall <- list(estimate = c(NA, NA), p.value = NA)
    cat("  Too few valid SA4s for Moran's I\n")
  }

  # --- B. Weekly Moran's I (sample of weeks) ---
  cat("  Computing weekly Moran's I (sample of 50 weeks)...\n")
  set.seed(42)
  sample_weeks <- sample(unique(df$week_start), min(50, n_distinct(df$week_start)))

  weekly_morans <- lapply(sample_weeks, function(wk) {
    resid_wk <- df |>
      filter(week_start == wk) |>
      group_by(sa4_code) |>
      summarise(resid = mean(dev_resid), .groups = "drop")

    resid_vec_wk <- resid_wk$resid[match(sa4_shp$SA4_CODE21[valid_idx],
                                          resid_wk$sa4_code)]
    v <- !is.na(resid_vec_wk)
    if (sum(v) < 10) return(c(I = NA_real_, p = NA_real_))

    nb_v <- tryCatch(subset(nb_sub, v), error = function(e) NULL)
    if (is.null(nb_v)) return(c(I = NA_real_, p = NA_real_))
    lw_v <- tryCatch(nb2listw(nb_v, style = "W", zero.policy = TRUE),
                     error = function(e) NULL)
    if (is.null(lw_v)) return(c(I = NA_real_, p = NA_real_))

    mt <- tryCatch(moran.test(resid_vec_wk[v], lw_v, zero.policy = TRUE),
                   error = function(e) NULL)
    if (is.null(mt)) return(c(I = NA_real_, p = NA_real_))
    c(I = unname(mt$estimate[1]), p = mt$p.value)
  })

  weekly_I <- sapply(weekly_morans, function(x) x["I"])
  weekly_p <- sapply(weekly_morans, function(x) x["p"])

  cat("  Weekly Moran's I — mean:", round(mean(weekly_I, na.rm = TRUE), 4),
      " median:", round(median(weekly_I, na.rm = TRUE), 4),
      " % significant (p<0.05):", round(100 * mean(weekly_p < 0.05, na.rm = TRUE), 1), "%\n")

  all_results[[grp]] <- data.frame(
    group                = grp,
    label                = group_labels[grp],
    moran_I_overall      = round(moran_overall$estimate[1], 4),
    moran_p_overall      = moran_overall$p.value,
    moran_I_weekly_mean  = round(mean(weekly_I, na.rm = TRUE), 4),
    moran_I_weekly_median = round(median(weekly_I, na.rm = TRUE), 4),
    pct_weeks_significant = round(100 * mean(weekly_p < 0.05, na.rm = TRUE), 1),
    n_weeks_tested       = sum(!is.na(weekly_I)),
    stringsAsFactors     = FALSE
  )

  # --- C. Moran scatter plot ---
  if (sum(valid) >= 10) {
    # Spatial lag of residuals
    lag_resid <- lag.listw(listw_sub, resid_sub, zero.policy = TRUE)

    moran_df <- data.frame(
      resid     = resid_sub,
      lag_resid = lag_resid,
      sa4       = sa4_shp$SA4_NAME21[valid_idx]
    )

    p_moran <- ggplot(moran_df, aes(x = resid, y = lag_resid)) +
      geom_point(alpha = 0.6, size = 2, colour = "steelblue") +
      geom_smooth(method = "lm", se = FALSE, colour = "#d73027", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(
        title = paste0(group_labels[grp], " — Moran scatter plot"),
        subtitle = paste0("Moran's I = ",
                          round(moran_overall$estimate[1], 4),
                          ", p = ", format.pval(moran_overall$p.value, digits = 3)),
        x = "Mean deviance residual",
        y = "Spatial lag of residual"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(fig_dir, paste0("spatial_moran_", grp, ".png")),
           p_moran, width = 7, height = 6, dpi = 300)
    cat("  -> Saved: spatial_moran_", grp, ".png\n")
  }
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

moran_table <- bind_rows(all_results)
write.csv(moran_table,
          file.path(tab_dir, "spatial_autocorrelation_tests.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/spatial_autocorrelation_tests.csv\n")


# ==============================================================================
# 4. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 20 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Spatial autocorrelation results:\n")
print(as.data.frame(moran_table), row.names = FALSE)

cat("\nInterpretation guide:\n")
cat("  - Moran's I near 0: no spatial autocorrelation\n")
cat("  - Moran's I > 0: positive spatial clustering of residuals\n")
cat("  - If I is small (<0.1) even if significant, practical impact is minimal\n")
cat("  - SA4 fixed effects already absorb spatial heterogeneity in levels\n")
cat("  - If substantial: consider spatial error model or cluster-robust SEs\n")

cat("\nOutputs:\n")
cat("  outputs/tables/spatial_autocorrelation_tests.csv\n")
cat("  outputs/figures/spatial_moran_*.png\n")
