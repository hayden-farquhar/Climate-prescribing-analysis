# ==============================================================================
# Script 31: Spatial CAR Model Sensitivity Analysis
# ==============================================================================
# Addresses spatial autocorrelation (Moran's I = 0.28–0.35) by fitting a
# Conditional Autoregressive (CAR) spatial random effects model alongside the
# standard DLNM. Compares estimates from:
#   1. Original quasi-Poisson GLM with SA4 fixed effects
#   2. Spatial GLMM with BYM2 (Besag-York-Mollié) random effects via INLA
#
# The BYM2 model explicitly accounts for spatial dependency between neighboring
# SA4 regions, producing confidence intervals that are not artificially narrow.
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#   data/spatial/SA4_2021/SA4_2021_AUST_GDA2020.shp
#
# Outputs:
#   outputs/tables/spatial_car_comparison.csv
#   outputs/tables/spatial_car_moran_after.csv
#   outputs/figures/spatial_car_forest.png
#   outputs/figures/spatial_car_er_comparison.png
#
# Dependencies:
#   install.packages("INLA",
#     repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#     dep = TRUE)
#   install.packages(c("dlnm", "splines", "sf", "spdep", "dplyr",
#                       "lubridate", "ggplot2", "patchwork"))
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
cat("Script 31: Spatial CAR Model Sensitivity Analysis\n")
cat("=" |> strrep(70), "\n\n")

# Check for INLA
if (!requireNamespace("INLA", quietly = TRUE)) {
  cat("ERROR: INLA package not installed.\n")
  cat("Install with:\n")
  cat('  install.packages("INLA",\n')
  cat('    repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),\n')
  cat("    dep = TRUE)\n")
  stop("INLA required for BYM2 spatial model.")
}
library(INLA)


# ==============================================================================
# 1. BUILD SPATIAL ADJACENCY GRAPH
# ==============================================================================

cat("Building spatial adjacency structure...\n")

sa4_shp <- st_read(file.path(spatial_dir, "SA4_2021", "SA4_2021_AUST_GDA2020.shp"),
                   quiet = TRUE)

panel <- read.csv(file.path(data_dir, "panel_weekly_sa4.csv"),
                  stringsAsFactors = FALSE)
panel$week_start <- as.Date(panel$week_start)
panel$sa4_code   <- as.character(panel$sa4_code)

analysis_sa4s <- unique(panel$sa4_code)

sa4_shp <- sa4_shp |>
  filter(SA4_CODE21 %in% analysis_sa4s) |>
  arrange(SA4_CODE21)

cat("  SA4s in analysis:", nrow(sa4_shp), "\n")

# Queen contiguity neighbours
nb <- poly2nb(sa4_shp, queen = TRUE)

# Fix islands (no contiguous neighbours)
n_islands <- sum(card(nb) == 0)
if (n_islands > 0) {
  coords <- st_centroid(sa4_shp) |> st_coordinates()
  knn <- knearneigh(coords, k = 1)
  nb_knn <- knn2nb(knn)
  for (i in which(card(nb) == 0)) {
    nb[[i]] <- nb_knn[[i]]
  }
  cat("  Fixed", n_islands, "island SA4s with k-nearest neighbour.\n")
}

# Create adjacency graph for INLA (save temp file)
graph_file <- tempfile(fileext = ".adj")
nb2INLA(graph_file, nb)
cat("  Adjacency graph created for INLA.\n")

# Create SA4 index mapping (for INLA spatial random effect)
sa4_index <- data.frame(
  sa4_code = sa4_shp$SA4_CODE21,
  sa4_idx  = seq_len(nrow(sa4_shp)),
  stringsAsFactors = FALSE
)

# Prepare time variables
panel <- panel |>
  mutate(
    year         = year(week_start),
    week_of_year = isoweek(week_start),
    time_index   = as.numeric(week_start - min(week_start)) / 7,
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  ) |>
  left_join(sa4_index, by = "sa4_code")

listw <- nb2listw(nb, style = "W", zero.policy = TRUE)


# ==============================================================================
# 2. FIT MODELS: ORIGINAL GLM vs BYM2 SPATIAL GLMM
# ==============================================================================

all_comparisons <- list()
all_moran_after <- list()

for (grp in primary_groups) {
  cat("\n", "=" |> strrep(50), "\n")
  cat("Processing:", group_labels[grp], "\n")
  cat("=" |> strrep(50), "\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(sa4_idx)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

  # --- Cross-basis ---
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

  temp_p90 <- quantile(df$tmax_mean, 0.90, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  # --- A. Original GLM (as in script 05) ---
  cat("  Fitting original quasi-Poisson GLM...\n")

  mod_glm <- glm(
    count ~ cb +
      ns(time_index, df = trend_df) +
      sin1 + cos1 + sin2 + cos2 +
      precip_total +
      factor(sa4_code),
    data = df,
    family = quasipoisson(link = "log")
  )

  pred_glm <- crosspred(cb, mod_glm,
                         at = c(temp_p90, temp_p95, temp_p99),
                         cen = temp_ref)

  cat("    Dispersion:", round(summary(mod_glm)$dispersion, 2), "\n")

  # --- B. BYM2 Spatial GLMM via INLA ---
  cat("  Fitting BYM2 spatial GLMM via INLA...\n")

  # Expand cross-basis into design matrix columns
  cb_mat <- as.matrix(cb)
  cb_names <- paste0("cb_", seq_len(ncol(cb_mat)))
  for (j in seq_len(ncol(cb_mat))) {
    df[[cb_names[j]]] <- cb_mat[, j]
  }

  # Trend spline columns
  trend_mat <- ns(df$time_index, df = trend_df)
  trend_names <- paste0("trend_", seq_len(ncol(trend_mat)))
  for (j in seq_len(ncol(trend_mat))) {
    df[[trend_names[j]]] <- trend_mat[, j]
  }

  # INLA formula with BYM2 spatial random effect
  # BYM2 = structured (Besag/ICAR) + unstructured (IID) spatial components
  fixed_terms <- paste(c(cb_names, trend_names,
                         "sin1", "cos1", "sin2", "cos2",
                         "precip_total"),
                       collapse = " + ")

  inla_formula <- as.formula(
    paste0("count ~ ", fixed_terms,
           " + f(sa4_idx, model = 'bym2', graph = '", graph_file, "')")
  )

  mod_inla <- tryCatch(
    inla(inla_formula,
         data = df,
         family = "poisson",
         control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
         control.predictor = list(compute = TRUE),
         verbose = FALSE),
    error = function(e) {
      cat("    INLA failed:", conditionMessage(e), "\n")
      cat("    Trying with simplified control...\n")
      inla(inla_formula,
           data = df,
           family = "poisson",
           control.compute = list(dic = TRUE, waic = TRUE),
           verbose = FALSE)
    }
  )

  cat("    INLA DIC:", round(mod_inla$dic$dic, 1), "\n")
  cat("    INLA WAIC:", round(mod_inla$waic$waic, 1), "\n")

  # Extract cross-basis coefficients from INLA
  inla_fixed <- mod_inla$summary.fixed
  cb_coefs_inla <- inla_fixed[cb_names, "mean"]
  cb_vcov_inla  <- matrix(0, length(cb_names), length(cb_names))
  for (j in seq_along(cb_names)) {
    cb_vcov_inla[j, j] <- inla_fixed[cb_names[j], "sd"]^2
  }

  # Compute cumulative RR at percentiles using INLA coefficients
  # We need to manually compute the cross-basis predictions
  pred_inla_rr <- list()
  for (pctl_name in c("p90", "p95", "p99")) {
    pctl_val <- get(paste0("temp_", pctl_name))
    # Basis values at this temperature (cumulative over all lags)
    bvar <- ns(pctl_val, knots = temp_knots,
               Boundary.knots = range(df$tmax_mean, na.rm = TRUE))
    blag_vals <- ns(0:max_lag, knots = logknots(max_lag, nk = 3),
                    Boundary.knots = c(0, max_lag))
    # Sum over lags for cumulative effect
    blag_sum <- colSums(blag_vals)

    # Reference basis values
    bvar_ref <- ns(temp_ref, knots = temp_knots,
                   Boundary.knots = range(df$tmax_mean, na.rm = TRUE))

    # Tensor product difference
    cb_at  <- as.vector(outer(bvar, blag_sum))
    cb_ref <- as.vector(outer(bvar_ref, blag_sum))
    cb_diff <- cb_at - cb_ref

    log_rr <- sum(cb_diff * cb_coefs_inla)
    var_rr <- as.numeric(t(cb_diff) %*% cb_vcov_inla %*% cb_diff)
    se_rr  <- sqrt(var_rr)

    pred_inla_rr[[pctl_name]] <- c(
      rr    = exp(log_rr),
      rr_lo = exp(log_rr - 1.96 * se_rr),
      rr_hi = exp(log_rr + 1.96 * se_rr)
    )
  }

  # --- C. Residual Moran's I after BYM2 ---
  cat("  Computing residual Moran's I from BYM2...\n")

  # INLA fitted values
  df$fitted_inla <- mod_inla$summary.fitted.values$mean[1:nrow(df)]
  df$resid_inla  <- (df$count - df$fitted_inla) / sqrt(pmax(df$fitted_inla, 1))

  mean_resid_inla <- df |>
    group_by(sa4_code) |>
    summarise(mean_resid = mean(resid_inla, na.rm = TRUE), .groups = "drop")

  sa4_order <- match(sa4_shp$SA4_CODE21, mean_resid_inla$sa4_code)
  resid_vec <- mean_resid_inla$mean_resid[sa4_order]
  valid <- !is.na(resid_vec)

  moran_after <- list(estimate = c(NA, NA), p.value = NA)
  if (sum(valid) >= 10) {
    nb_sub <- subset(nb, valid)
    listw_sub <- nb2listw(nb_sub, style = "W", zero.policy = TRUE)
    moran_after <- moran.test(resid_vec[valid], listw_sub, zero.policy = TRUE)
    cat("    Moran's I (BYM2 residuals):", round(moran_after$estimate[1], 4),
        "  p =", format.pval(moran_after$p.value, digits = 3), "\n")
  }

  # Also compute Moran's I from original GLM for comparison
  df$resid_glm <- as.numeric(residuals(mod_glm, type = "deviance"))
  mean_resid_glm <- df |>
    group_by(sa4_code) |>
    summarise(mean_resid = mean(resid_glm, na.rm = TRUE), .groups = "drop")

  resid_vec_glm <- mean_resid_glm$mean_resid[sa4_order]
  valid_glm <- !is.na(resid_vec_glm)
  moran_before <- list(estimate = c(NA, NA), p.value = NA)
  if (sum(valid_glm) >= 10) {
    nb_sub_g <- subset(nb, valid_glm)
    listw_sub_g <- nb2listw(nb_sub_g, style = "W", zero.policy = TRUE)
    moran_before <- moran.test(resid_vec_glm[valid_glm], listw_sub_g, zero.policy = TRUE)
    cat("    Moran's I (GLM residuals):", round(moran_before$estimate[1], 4),
        "  p =", format.pval(moran_before$p.value, digits = 3), "\n")
  }

  # --- Store results ---
  comp_row <- data.frame(
    group = grp,
    label = group_labels[grp],
    # GLM results
    glm_rr_p90    = round(pred_glm$allRRfit[1], 4),
    glm_rr_p90_lo = round(pred_glm$allRRlow[1], 4),
    glm_rr_p90_hi = round(pred_glm$allRRhigh[1], 4),
    glm_rr_p95    = round(pred_glm$allRRfit[2], 4),
    glm_rr_p95_lo = round(pred_glm$allRRlow[2], 4),
    glm_rr_p95_hi = round(pred_glm$allRRhigh[2], 4),
    glm_rr_p99    = round(pred_glm$allRRfit[3], 4),
    glm_rr_p99_lo = round(pred_glm$allRRlow[3], 4),
    glm_rr_p99_hi = round(pred_glm$allRRhigh[3], 4),
    # BYM2 results
    bym2_rr_p90    = round(pred_inla_rr$p90["rr"], 4),
    bym2_rr_p90_lo = round(pred_inla_rr$p90["rr_lo"], 4),
    bym2_rr_p90_hi = round(pred_inla_rr$p90["rr_hi"], 4),
    bym2_rr_p95    = round(pred_inla_rr$p95["rr"], 4),
    bym2_rr_p95_lo = round(pred_inla_rr$p95["rr_lo"], 4),
    bym2_rr_p95_hi = round(pred_inla_rr$p95["rr_hi"], 4),
    bym2_rr_p99    = round(pred_inla_rr$p99["rr"], 4),
    bym2_rr_p99_lo = round(pred_inla_rr$p99["rr_lo"], 4),
    bym2_rr_p99_hi = round(pred_inla_rr$p99["rr_hi"], 4),
    # Model fit
    inla_dic  = round(mod_inla$dic$dic, 1),
    inla_waic = round(mod_inla$waic$waic, 1),
    glm_dispersion = round(summary(mod_glm)$dispersion, 2),
    n_obs = nrow(df),
    stringsAsFactors = FALSE
  )
  all_comparisons[[grp]] <- comp_row

  all_moran_after[[grp]] <- data.frame(
    group = grp,
    label = group_labels[grp],
    moran_I_glm  = round(moran_before$estimate[1], 4),
    moran_p_glm  = moran_before$p.value,
    moran_I_bym2 = round(moran_after$estimate[1], 4),
    moran_p_bym2 = moran_after$p.value,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

comp_table <- bind_rows(all_comparisons)
write.csv(comp_table,
          file.path(tab_dir, "spatial_car_comparison.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/spatial_car_comparison.csv\n")

moran_table <- bind_rows(all_moran_after)
write.csv(moran_table,
          file.path(tab_dir, "spatial_car_moran_after.csv"),
          row.names = FALSE)
cat("-> Saved: outputs/tables/spatial_car_moran_after.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison plots...\n")

# --- Forest plot: GLM vs BYM2 at 95th percentile ---
forest_data <- comp_table |>
  tidyr::pivot_longer(
    cols = matches("^(glm|bym2)_rr_p95"),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    model = ifelse(grepl("^glm", metric), "GLM (SA4 fixed effects)",
                   "BYM2 (spatial random effects)"),
    stat  = gsub("^(glm|bym2)_rr_p95_?", "", metric),
    stat  = ifelse(stat == "", "rr", stat)
  ) |>
  select(group, label, model, stat, value) |>
  tidyr::pivot_wider(names_from = stat, values_from = value) |>
  mutate(label_model = paste0(label, "\n", model))

p_forest <- ggplot(forest_data, aes(x = rr, y = label_model, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_colour_manual(
    values = c("GLM (SA4 fixed effects)" = "#4575b4",
               "BYM2 (spatial random effects)" = "#d73027"),
    name = "Model"
  ) +
  labs(
    title = "GLM vs BYM2 spatial model: cumulative RR at 95th percentile",
    subtitle = "BYM2 accounts for spatial autocorrelation via Besag-York-Mollié random effects",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "spatial_car_forest.png"),
       p_forest, width = 10, height = 7, dpi = 300)
cat("  -> Saved: spatial_car_forest.png\n")

# --- Moran's I before/after comparison ---
moran_long <- moran_table |>
  tidyr::pivot_longer(
    cols = starts_with("moran_I"),
    names_to = "model",
    values_to = "moran_I"
  ) |>
  mutate(
    model = ifelse(grepl("glm", model), "GLM", "BYM2"),
    model = factor(model, levels = c("GLM", "BYM2"))
  )

p_moran <- ggplot(moran_long, aes(x = model, y = moran_I, fill = model)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0, colour = "grey40") +
  facet_wrap(~ label) +
  scale_fill_manual(values = c(GLM = "#4575b4", BYM2 = "#d73027"), guide = "none") +
  labs(
    title = "Residual spatial autocorrelation: GLM vs BYM2",
    subtitle = "Moran's I of time-averaged residuals (lower = less spatial clustering)",
    x = NULL,
    y = "Moran's I"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(fig_dir, "spatial_car_er_comparison.png"),
       p_moran, width = 10, height = 5, dpi = 300)
cat("  -> Saved: spatial_car_er_comparison.png\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 31 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Comparison of RR estimates at 95th percentile:\n")
for (grp in primary_groups) {
  r <- comp_table[comp_table$group == grp, ]
  cat(sprintf("  %-15s GLM: %.3f (%.3f-%.3f)  BYM2: %.3f (%.3f-%.3f)\n",
              r$label,
              r$glm_rr_p95, r$glm_rr_p95_lo, r$glm_rr_p95_hi,
              r$bym2_rr_p95, r$bym2_rr_p95_lo, r$bym2_rr_p95_hi))
}

cat("\nMoran's I comparison:\n")
print(as.data.frame(moran_table), row.names = FALSE)

cat("\nInterpretation:\n")
cat("  - If BYM2 CIs are wider but point estimates similar: spatial autocorrelation\n")
cat("    inflated precision but not estimates — original findings robust.\n")
cat("  - If BYM2 CIs similar: spatial autocorrelation adequately handled by SA4 FE.\n")
cat("  - BYM2 Moran's I should be closer to 0 than GLM Moran's I.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/spatial_car_comparison.csv\n")
cat("  outputs/tables/spatial_car_moran_after.csv\n")
cat("  outputs/figures/spatial_car_forest.png\n")
cat("  outputs/figures/spatial_car_er_comparison.png\n")
