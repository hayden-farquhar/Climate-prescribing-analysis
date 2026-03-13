# ==============================================================================
# Script 18: Cross-Basis Sensitivity Analysis
# ==============================================================================
# Tests whether primary DLNM results are sensitive to the choice of:
#   A. Temperature spline knot placement (p10/p50/p90 vs p25/p50/p75 vs p33/p66)
#   B. Temperature spline degrees of freedom (3 vs 4 vs 5 df)
#   C. Lag spline specification (2 vs 3 vs 4 log-knots)
#   D. Linear temperature term (no spline)
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/crossbasis_sensitivity.csv
#   outputs/figures/crossbasis_sensitivity_forest.png
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

cat("=" |> strrep(70), "\n")
cat("Script 18: Cross-Basis Sensitivity Analysis\n")
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
    year         = year(week_start),
    sin1 = sin(2 * pi * week_of_year / 52),
    cos1 = cos(2 * pi * week_of_year / 52),
    sin2 = sin(4 * pi * week_of_year / 52),
    cos2 = cos(4 * pi * week_of_year / 52)
  )


# ==============================================================================
# 2. DEFINE SPECIFICATIONS
# ==============================================================================

# Each specification is a list with argvar and arglag for crossbasis()
make_specs <- function(df) {
  tmax <- df$tmax_mean

  specs <- list(
    # --- A. Primary specification (reference) ---
    "Primary (ns, knots p10/p50/p90, 3 lag knots)" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.10, 0.50, 0.90), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),

    # --- B. Alternative knot placements ---
    "Knots p25/p50/p75" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.25, 0.50, 0.75), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),
    "Knots p33/p66" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.33, 0.66), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),
    "Knots p05/p25/p75/p95" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.05, 0.25, 0.75, 0.95), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),

    # --- C. Alternative df ---
    "ns df=3 (no internal knots)" = list(
      argvar = list(fun = "ns", df = 3),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),
    "ns df=5" = list(
      argvar = list(fun = "ns", df = 5),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    ),

    # --- D. Alternative lag specifications ---
    "2 lag knots" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.10, 0.50, 0.90), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 2))
    ),
    "4 lag knots" = list(
      argvar = list(fun = "ns",
                    knots = quantile(tmax, c(0.10, 0.50, 0.90), na.rm = TRUE)),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 4))
    ),

    # --- E. Linear temperature ---
    "Linear temperature" = list(
      argvar = list(fun = "lin"),
      arglag = list(fun = "ns", knots = logknots(max_lag, nk = 3))
    )
  )

  return(specs)
}


# ==============================================================================
# 3. FIT ALL SPECIFICATIONS
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Fitting", length(primary_groups), "groups x 9 specifications = ",
    length(primary_groups) * 9, "models\n")
cat("=" |> strrep(70), "\n")

all_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))
  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)
  temp_p99 <- quantile(df$tmax_mean, 0.99, na.rm = TRUE)

  specs <- make_specs(df)

  for (spec_name in names(specs)) {
    cat("  ", spec_name, "... ")

    spec <- specs[[spec_name]]

    cb <- tryCatch(
      crossbasis(df$tmax_mean, lag = max_lag,
                 argvar = spec$argvar, arglag = spec$arglag),
      error = function(e) NULL
    )
    if (is.null(cb)) { cat("cross-basis failed\n"); next }

    model <- tryCatch(
      glm(count ~ cb +
            ns(time_index, df = trend_df) +
            sin1 + cos1 + sin2 + cos2 +
            precip_total +
            factor(sa4_code),
          data = df,
          family = quasipoisson(link = "log")),
      error = function(e) NULL
    )
    if (is.null(model)) { cat("model failed\n"); next }

    pred_pctl <- tryCatch(
      crosspred(cb, model, at = c(temp_p95, temp_p99), cen = temp_ref),
      error = function(e) NULL
    )
    if (is.null(pred_pctl)) { cat("prediction failed\n"); next }

    # Quasi-AIC (QIC): use deviance + 2*k*dispersion
    qaic <- model$deviance + 2 * length(coef(model)) * summary(model)$dispersion

    result <- data.frame(
      group        = grp,
      label        = group_labels[grp],
      specification = spec_name,
      rr_p95       = round(pred_pctl$allRRfit[1], 4),
      rr_p95_lo    = round(pred_pctl$allRRlow[1], 4),
      rr_p95_hi    = round(pred_pctl$allRRhigh[1], 4),
      rr_p99       = round(pred_pctl$allRRfit[2], 4),
      rr_p99_lo    = round(pred_pctl$allRRlow[2], 4),
      rr_p99_hi    = round(pred_pctl$allRRhigh[2], 4),
      n_params     = length(coef(model)),
      dispersion   = round(summary(model)$dispersion, 2),
      qaic         = round(qaic, 0),
      stringsAsFactors = FALSE
    )

    all_results[[paste0(grp, "_", spec_name)]] <- result
    cat("RR@p95 =", result$rr_p95,
        "(", result$rr_p95_lo, "-", result$rr_p95_hi, ")",
        " QAIC =", result$qaic, "\n")
  }
}


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

sens_table <- bind_rows(all_results)
write.csv(sens_table,
          file.path(tab_dir, "crossbasis_sensitivity.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/crossbasis_sensitivity.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating sensitivity forest plots...\n")

for (grp in primary_groups) {
  df_plot <- sens_table |>
    filter(group == grp) |>
    mutate(
      is_primary = grepl("Primary", specification),
      spec_f = factor(specification, levels = rev(unique(specification)))
    )

  p <- ggplot(df_plot, aes(x = rr_p95, y = spec_f)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5, colour = "steelblue") +
    geom_point(aes(colour = is_primary), size = 2.5) +
    scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "steelblue"),
                        guide = "none") +
    labs(
      title = paste0(group_labels[grp], " — Cross-basis sensitivity"),
      subtitle = "Red = primary specification; blue = alternatives",
      x = "Cumulative RR at 95th percentile (95% CI)",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))

  ggsave(file.path(fig_dir, paste0("crossbasis_sensitivity_", grp, ".png")),
         p, width = 10, height = 6, dpi = 300)
}
cat("  -> Saved: crossbasis_sensitivity_*.png\n")

# Combined panel
plots <- list()
for (grp in primary_groups) {
  df_plot <- sens_table |>
    filter(group == grp) |>
    mutate(
      is_primary = grepl("Primary", specification),
      spec_f = factor(specification, levels = rev(unique(specification)))
    )
  plots[[grp]] <- ggplot(df_plot, aes(x = rr_p95, y = spec_f)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                   height = 0.3, linewidth = 0.5, colour = "steelblue") +
    geom_point(aes(colour = is_primary), size = 2) +
    scale_colour_manual(values = c("TRUE" = "#d73027", "FALSE" = "steelblue"),
                        guide = "none") +
    labs(title = group_labels[grp], x = "Cumulative RR (95% CI)", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

combined <- wrap_plots(plots, ncol = 3)
ggsave(file.path(fig_dir, "crossbasis_sensitivity_panel.png"),
       combined, width = 18, height = 6, dpi = 300)
cat("  -> Saved: crossbasis_sensitivity_panel.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 18 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Cross-basis sensitivity (RR at p95):\n")
print(as.data.frame(
  sens_table |> select(label, specification, rr_p95, rr_p95_lo, rr_p95_hi, qaic)
), row.names = FALSE)

cat("\nOutputs:\n")
cat("  outputs/tables/crossbasis_sensitivity.csv\n")
cat("  outputs/figures/crossbasis_sensitivity_*.png\n")
