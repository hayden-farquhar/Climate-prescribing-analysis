# ==============================================================================
# Script 29: Newey-West (HAC) Standard Errors
# ==============================================================================
# Computes heteroskedasticity and autocorrelation consistent (HAC) standard
# errors for the primary DLNM cumulative associations using the sandwich
# estimator with Newey-West kernel.
#
# Motivation:
#   Script 17 found lag-1 ACF of 0.29-0.55 in model residuals. Standard
#   quasi-Poisson SEs assume independent observations. If residual
#   autocorrelation is present, CIs may be too narrow. Newey-West SEs
#   account for serial correlation up to a specified lag bandwidth.
#
# Approach:
#   1. Fit primary DLNM models
#   2. Compute Newey-West (HAC) variance-covariance matrix using sandwich pkg
#   3. Re-derive cumulative RR and CI from HAC-adjusted vcov
#   4. Compare standard vs HAC CIs
#
# Inputs:
#   data/processed/panel_weekly_sa4.csv
#
# Outputs:
#   outputs/tables/newey_west_comparison.csv
#   outputs/figures/newey_west_forest.png
#
# Dependencies:
#   install.packages(c("dlnm", "splines", "dplyr", "lubridate",
#                       "sandwich", "lmtest", "ggplot2"))
# ==============================================================================

library(dlnm)
library(splines)
library(dplyr)
library(lubridate)
library(sandwich)
library(lmtest)
library(ggplot2)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
fig_dir     <- file.path(project_dir, "outputs", "figures")
tab_dir     <- file.path(project_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

max_lag <- 12

primary_groups <- c("cvd_total", "mh_total", "resp_total",
                    "resp_relievers", "resp_preventers")
group_labels <- c(
  cvd_total       = "Cardiovascular",
  mh_total        = "Mental Health",
  resp_total      = "Respiratory",
  resp_relievers  = "Resp. Relievers",
  resp_preventers = "Resp. Preventers"
)

cat("=" |> strrep(70), "\n")
cat("Script 29: Newey-West HAC Standard Errors\n")
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
# 2. FIT MODELS AND COMPUTE HAC SEs
# ==============================================================================

all_results <- list()

for (grp in primary_groups) {
  cat("\n--- ", group_labels[grp], " ---\n")

  df <- panel |>
    filter(analysis_group == grp,
           suppressed != "True",
           !is.na(count),
           !is.na(tmax_mean),
           !is.na(precip_total)) |>
    arrange(sa4_code, week_start)

  cat("  Observations:", format(nrow(df), big.mark = ","), "\n")

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

  temp_ref <- median(df$tmax_mean, na.rm = TRUE)
  temp_p95 <- quantile(df$tmax_mean, 0.95, na.rm = TRUE)

  # --- Standard crosspred (model-based SEs) ---
  pred_std <- crosspred(cb, model, at = temp_p95, cen = temp_ref)
  rr_std    <- pred_std$allRRfit[1]
  rr_std_lo <- pred_std$allRRlow[1]
  rr_std_hi <- pred_std$allRRhigh[1]

  cat("  Standard:   RR =", round(rr_std, 4),
      "(", round(rr_std_lo, 4), "-", round(rr_std_hi, 4), ")\n")

  # --- Newey-West HAC vcov ---
  # Bandwidth: use max_lag (12 weeks) as the autocorrelation window
  # This accounts for serial correlation up to 12 weeks
  nw_bandwidth <- max_lag

  # Compute HAC vcov using sandwich package
  # For panel data ordered by (sa4, time), vcovHAC captures within-unit
  # serial correlation
  vcov_hac <- tryCatch(
    vcovHAC(model, order.by = df$time_index, prewhite = FALSE,
            kernel = "Bartlett",
            bw = nw_bandwidth,
            sandwich = TRUE),
    error = function(e) {
      cat("  HAC failed:", e$message, "\n")
      cat("  Trying NeweyWest...\n")
      tryCatch(
        NeweyWest(model, lag = nw_bandwidth, prewhite = FALSE),
        error = function(e2) {
          cat("  NeweyWest also failed:", e2$message, "\n")
          NULL
        }
      )
    }
  )

  if (is.null(vcov_hac)) {
    cat("  Skipping HAC for this group.\n")
    next
  }

  # --- Re-derive cumulative RR with HAC vcov ---
  # The cumulative association at a temperature is a linear combination
  # of cross-basis coefficients. We need the basis values at p95 and
  # the corresponding coefficients + HAC vcov.

  # Get crossbasis coefficient names and positions
  cb_names <- grep("^cb", names(coef(model)), value = TRUE)
  cb_idx   <- which(names(coef(model)) %in% cb_names)
  cb_coefs <- coef(model)[cb_idx]

  # Create the basis values at p95 (centering at median)
  # The cumulative association sums over all lags
  # crosspred stores this in the $allfit attribute
  # We can reconstruct: log(RR) = sum over lags of basis_values * coefs
  # For the cumulative effect, we need the full basis at the target temp
  # summed across lags

  # Use the internal crosspred approach: get the basis matrix at p95
  # and sum across lag dimension
  b_temp <- do.call("onebasis", c(list(x = temp_p95),
                                    attr(cb, "argvar")))
  b_ref  <- do.call("onebasis", c(list(x = temp_ref),
                                    attr(cb, "argvar")))
  b_lag  <- do.call("onebasis", c(list(x = 0:max_lag),
                                    attr(cb, "arglag")))

  # The cumulative prediction vector (length = ncol(cb))
  # is the Kronecker product of (b_temp - b_ref) and colSums(b_lag)
  bvar_diff <- b_temp - b_ref
  blag_sum  <- colSums(b_lag)
  cum_vector <- as.numeric(bvar_diff %x% blag_sum)  # Kronecker product

  # log(RR) cumulative
  log_rr_hac <- sum(cum_vector * cb_coefs)

  # SE from HAC vcov (subset to cb coefficients)
  vcov_cb_hac <- vcov_hac[cb_idx, cb_idx]
  se_hac <- sqrt(as.numeric(t(cum_vector) %*% vcov_cb_hac %*% cum_vector))

  # Standard SE for comparison
  vcov_cb_std <- vcov(model)[cb_idx, cb_idx]
  se_std <- sqrt(as.numeric(t(cum_vector) %*% vcov_cb_std %*% cum_vector))

  # RR and CI
  rr_hac    <- exp(log_rr_hac)
  rr_hac_lo <- exp(log_rr_hac - 1.96 * se_hac)
  rr_hac_hi <- exp(log_rr_hac + 1.96 * se_hac)

  # CI width ratio
  ci_width_std <- log(rr_std_hi) - log(rr_std_lo)
  ci_width_hac <- log(rr_hac_hi) - log(rr_hac_lo)
  ci_ratio <- ci_width_hac / ci_width_std

  cat("  HAC (NW):   RR =", round(rr_hac, 4),
      "(", round(rr_hac_lo, 4), "-", round(rr_hac_hi, 4), ")\n")
  cat("  SE ratio (HAC/std):", round(se_hac / se_std, 3), "\n")
  cat("  CI width ratio:", round(ci_ratio, 3), "\n")

  all_results[[grp]] <- data.frame(
    group          = grp,
    label          = group_labels[grp],
    rr_std         = round(rr_std, 4),
    rr_std_lo      = round(rr_std_lo, 4),
    rr_std_hi      = round(rr_std_hi, 4),
    se_std         = round(se_std, 5),
    rr_hac         = round(rr_hac, 4),
    rr_hac_lo      = round(rr_hac_lo, 4),
    rr_hac_hi      = round(rr_hac_hi, 4),
    se_hac         = round(se_hac, 5),
    se_ratio       = round(se_hac / se_std, 3),
    ci_width_ratio = round(ci_ratio, 3),
    nw_bandwidth   = nw_bandwidth,
    sig_std        = rr_std_lo > 1 | rr_std_hi < 1,
    sig_hac        = rr_hac_lo > 1 | rr_hac_hi < 1,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

nw_table <- bind_rows(all_results)
write.csv(nw_table,
          file.path(tab_dir, "newey_west_comparison.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/newey_west_comparison.csv\n")


# ==============================================================================
# 4. VISUALISATION
# ==============================================================================

cat("\nGenerating comparison forest plot...\n")

# Reshape for plotting
plot_std <- nw_table |>
  transmute(label, se_type = "Model-based", rr = rr_std,
            rr_lo = rr_std_lo, rr_hi = rr_std_hi)
plot_hac <- nw_table |>
  transmute(label, se_type = "Newey-West HAC", rr = rr_hac,
            rr_lo = rr_hac_lo, rr_hi = rr_hac_hi)
plot_df <- bind_rows(plot_std, plot_hac)

p <- ggplot(plot_df, aes(x = rr, y = label, colour = se_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = rr_lo, xmax = rr_hi),
                 height = 0.3, linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Model-based" = "steelblue",
                                  "Newey-West HAC" = "#d73027"),
                      name = "Standard Error") +
  labs(
    title = "Cumulative RR at 95th percentile: model-based vs Newey-West SEs",
    subtitle = paste0("Newey-West bandwidth = ", max_lag, " weeks"),
    x = "Cumulative RR (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "newey_west_forest.png"),
       p, width = 10, height = 5, dpi = 300)
cat("  -> Saved: newey_west_forest.png\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 29 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Newey-West HAC comparison (RR at p95):\n")
print(as.data.frame(
  nw_table |>
    select(label, rr_std, rr_std_lo, rr_std_hi,
           rr_hac, rr_hac_lo, rr_hac_hi, se_ratio, sig_std, sig_hac)
), row.names = FALSE)

cat("\nInterpretation:\n")
cat("  SE ratio > 1: HAC SEs are wider (model-based SEs underestimate uncertainty)\n")
cat("  SE ratio ~ 1: autocorrelation has minimal impact on inference\n")
cat("  If sig_std == TRUE but sig_hac == FALSE: finding may not survive\n")
cat("    correction for serial correlation.\n")
cat("  If both agree: finding is robust to autocorrelation correction.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/newey_west_comparison.csv\n")
cat("  outputs/figures/newey_west_forest.png\n")
