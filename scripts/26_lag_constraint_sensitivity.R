# ==============================================================================
# Script 26: Distributed Lag Constraint Sensitivity
# ==============================================================================
# Compares spline-constrained lag structure (primary) against unconstrained
# (indicator) lag models and alternative lag specifications to verify that
# the lag shape assumption does not drive results.
#
# Tests:
#   1. Unconstrained (indicator) lag — one parameter per lag week
#   2. Polynomial DL (degree 3) — classic Almon lag
#   3. Shorter lag (0-8 weeks) with spline
#   4. Longer lag (0-16 weeks) with spline
#   5. Primary model (0-12 weeks, ns lag with 3 log-knots) as reference
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/lag_constraint_sensitivity.csv
#   outputs/figures/lag_constraint_comparison.png
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

primary_groups <- c("cvd_total", "mh_total", "resp_total")
group_labels <- c(
  cvd_total  = "Cardiovascular",
  mh_total   = "Mental Health",
  resp_total = "Respiratory"
)

cat("=" |> strrep(70), "\n")
cat("Script 26: Distributed Lag Constraint Sensitivity\n")
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
# 2. LAG SPECIFICATIONS
# ==============================================================================

lag_specs <- list(
  primary = list(
    name = "Primary (ns lag, 0-12w)",
    max_lag = 12,
    arglag = list(fun = "ns", knots = logknots(12, nk = 3))
  ),
  unconstrained = list(
    name = "Unconstrained (indicator, 0-12w)",
    max_lag = 12,
    arglag = list(fun = "integer")
  ),
  poly3 = list(
    name = "Polynomial DL (degree 3, 0-12w)",
    max_lag = 12,
    arglag = list(fun = "poly", degree = 3)
  ),
  short_lag = list(
    name = "Shorter lag (ns, 0-8w)",
    max_lag = 8,
    arglag = list(fun = "ns", knots = logknots(8, nk = 2))
  ),
  long_lag = list(
    name = "Longer lag (ns, 0-16w)",
    max_lag = 16,
    arglag = list(fun = "ns", knots = logknots(16, nk = 3))
  )
)


# ==============================================================================
# 3. FIT MODELS AND COMPARE
# ==============================================================================

all_results <- list()
all_lag_curves <- list()

for (grp in primary_groups) {
  cat("\n", "=" |> strrep(50), "\n")
  cat(group_labels[grp], "\n")
  cat("=" |> strrep(50), "\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean)) |>
    arrange(sa4_code, week_start)

  temp_knots <- quantile(df$tmax_mean, c(0.10, 0.50, 0.90), na.rm = TRUE)
  n_years  <- n_distinct(df$year)
  trend_df <- max(2, round(n_years * 2))
  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)

  for (spec_name in names(lag_specs)) {
    spec <- lag_specs[[spec_name]]
    cat("  ", spec$name, "...\n")

    cb <- tryCatch(
      crossbasis(
        df$tmax_mean,
        lag = spec$max_lag,
        argvar = list(fun = "ns", knots = temp_knots),
        arglag = spec$arglag
      ),
      error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
    )
    if (is.null(cb)) next

    model <- tryCatch(
      glm(
        count ~ cb +
          ns(time_index, df = trend_df) +
          sin1 + cos1 + sin2 + cos2 +
          precip_total +
          factor(sa4_code),
        data = df,
        family = quasipoisson(link = "log")
      ),
      error = function(e) { cat("    ERROR:", e$message, "\n"); NULL }
    )
    if (is.null(model)) next

    # Cumulative RR at p95
    pred_cum <- crosspred(cb, model, at = temp_p95, cen = temp_ref)

    # Lag-specific RR at p95
    lag_rr <- pred_cum$matRRfit[1, ]
    lag_lo <- pred_cum$matRRlow[1, ]
    lag_hi <- pred_cum$matRRhigh[1, ]

    all_lag_curves[[paste0(grp, "_", spec_name)]] <- data.frame(
      group     = grp,
      label     = group_labels[grp],
      spec      = spec$name,
      spec_key  = spec_name,
      lag       = 0:spec$max_lag,
      rr        = lag_rr,
      rr_lo     = lag_lo,
      rr_hi     = lag_hi,
      stringsAsFactors = FALSE
    )

    # QAIC-like measure
    ll <- -0.5 * deviance(model)
    disp <- summary(model)$dispersion
    k <- length(coef(model))
    qaic <- -2 * ll / disp + 2 * k

    all_results[[paste0(grp, "_", spec_name)]] <- data.frame(
      group       = grp,
      label       = group_labels[grp],
      spec        = spec$name,
      spec_key    = spec_name,
      max_lag     = spec$max_lag,
      n_cb_params = ncol(cb),
      rr_p95      = round(pred_cum$allRRfit[1], 4),
      rr_p95_lo   = round(pred_cum$allRRlow[1], 4),
      rr_p95_hi   = round(pred_cum$allRRhigh[1], 4),
      dispersion  = round(disp, 2),
      qaic        = round(qaic, 1),
      stringsAsFactors = FALSE
    )

    cat("    RR at p95:", round(pred_cum$allRRfit[1], 4),
        "(", round(pred_cum$allRRlow[1], 4), "-",
        round(pred_cum$allRRhigh[1], 4), ")\n")
    cat("    CB params:", ncol(cb), " QAIC:", round(qaic, 1), "\n")
  }
}


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

results_table <- bind_rows(all_results)
write.csv(results_table,
          file.path(tab_dir, "lag_constraint_sensitivity.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/lag_constraint_sensitivity.csv\n")


# ==============================================================================
# 5. VISUALISATION
# ==============================================================================

cat("\nGenerating lag-response comparison plots...\n")

lag_df <- bind_rows(all_lag_curves)

for (grp in primary_groups) {
  df_plot <- lag_df |> filter(group == grp)

  p <- ggplot(df_plot, aes(x = lag, y = rr, colour = spec, fill = spec)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    labs(
      title = paste0(group_labels[grp],
                     " — Lag-response at 95th percentile"),
      subtitle = "Comparing lag constraint specifications",
      x = "Lag (weeks)",
      y = "RR per lag week",
      colour = "Specification",
      fill = "Specification"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 8)) +
    guides(colour = guide_legend(ncol = 2),
           fill = guide_legend(ncol = 2))

  ggsave(file.path(fig_dir, paste0("lag_constraint_", grp, ".png")),
         p, width = 9, height = 6, dpi = 300)
  cat("  -> Saved: lag_constraint_", grp, ".png\n")
}

# Combined cumulative RR comparison
p_forest <- ggplot(results_table,
                    aes(x = rr_p95,
                        y = reorder(paste0(label, ": ", spec), rr_p95))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_p95_lo, xmax = rr_p95_hi),
                 height = 0.3, linewidth = 0.5, colour = "steelblue") +
  geom_point(size = 2, colour = "steelblue") +
  facet_wrap(~ label, scales = "free_y", ncol = 1) +
  labs(
    title = "Cumulative RR at p95 across lag specifications",
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "lag_constraint_comparison.png"),
       p_forest, width = 10, height = 10, dpi = 300)
cat("  -> Saved: lag_constraint_comparison.png\n")


# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 26 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Lag constraint sensitivity (cumulative RR at p95):\n")
print(as.data.frame(
  results_table |>
    select(label, spec, rr_p95, rr_p95_lo, rr_p95_hi, n_cb_params, qaic)
), row.names = FALSE)

cat("\nInterpretation:\n")
cat("  If cumulative RR is consistent across specifications, the lag shape\n")
cat("  assumption does not drive results.\n")
cat("  Unconstrained model provides reference but may be noisy.\n")
cat("  Lower QAIC indicates better fit-complexity trade-off.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/lag_constraint_sensitivity.csv\n")
cat("  outputs/figures/lag_constraint_*.png\n")
cat("  outputs/figures/lag_constraint_comparison.png\n")
