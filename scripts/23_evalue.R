# ==============================================================================
# Script 23: E-Value for Unmeasured Confounding
# ==============================================================================
# Computes E-values for the primary DLNM findings to quantify the strength
# of unmeasured confounding needed to explain away the observed associations.
#
# The E-value is the minimum strength of association (on the RR scale) that
# an unmeasured confounder would need to have with both the exposure and
# outcome to fully explain the observed RR.
#
# Reference: VanderWeele & Ding (2017) Annals of Internal Medicine
#
# Inputs:
#   outputs/tables/dlnm_results_summary.csv
#   outputs/tables/bushfire_did_results.csv
#
# Outputs:
#   outputs/tables/evalues.csv
#
# Dependencies:
#   install.packages(c("dplyr"))
# ==============================================================================

library(dplyr)

project_dir <- here::here()
tab_dir     <- file.path(project_dir, "outputs", "tables")

cat("=" |> strrep(70), "\n")
cat("Script 23: E-Value for Unmeasured Confounding\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. E-VALUE FUNCTION
# ==============================================================================

# E-value formula (VanderWeele & Ding 2017):
# For RR >= 1: E = RR + sqrt(RR * (RR - 1))
# For RR < 1:  compute on 1/RR, which gives the same result
# For the CI limit closest to null: same formula on that limit

compute_evalue <- function(rr, rr_lo = NULL, rr_hi = NULL) {
  # Convert protective associations (RR < 1) to the equivalent RR > 1
  if (rr < 1) {
    rr_use <- 1 / rr
    # For protective: the CI limit closest to null is the upper bound
    ci_use <- if (!is.null(rr_hi)) 1 / rr_hi else NULL
  } else {
    rr_use <- rr
    # For harmful: the CI limit closest to null is the lower bound
    ci_use <- if (!is.null(rr_lo)) rr_lo else NULL
  }

  # E-value for point estimate
  e_point <- rr_use + sqrt(rr_use * (rr_use - 1))

  # E-value for CI limit closest to null
  e_ci <- NA
  if (!is.null(ci_use)) {
    if (ci_use < 1) {
      # CI crosses null — E-value for CI = 1
      e_ci <- 1
    } else {
      e_ci <- ci_use + sqrt(ci_use * (ci_use - 1))
    }
  }

  return(list(e_point = round(e_point, 2), e_ci = round(e_ci, 2)))
}


# ==============================================================================
# 2. COMPUTE E-VALUES FOR DLNM RESULTS
# ==============================================================================

cat("Computing E-values for primary DLNM findings...\n\n")

evalue_results <- list()

# --- Primary DLNM results ---
dlnm_file <- file.path(tab_dir, "dlnm_results_summary.csv")
if (file.exists(dlnm_file)) {
  dlnm <- read.csv(dlnm_file, stringsAsFactors = FALSE)

  for (i in seq_len(nrow(dlnm))) {
    r <- dlnm[i, ]
    ev <- compute_evalue(r$rr_p95, r$rr_p95_lo, r$rr_p95_hi)

    evalue_results[[length(evalue_results) + 1]] <- data.frame(
      source      = "Primary DLNM",
      outcome     = r$label,
      comparison  = "p95 vs median",
      rr          = r$rr_p95,
      rr_lo       = r$rr_p95_lo,
      rr_hi       = r$rr_p95_hi,
      e_point     = ev$e_point,
      e_ci        = ev$e_ci,
      stringsAsFactors = FALSE
    )

    cat("  ", r$label, ": RR =", r$rr_p95,
        " -> E-value =", ev$e_point,
        "(CI:", ev$e_ci, ")\n")
  }
}

# --- Bushfire DiD results ---
bush_file <- file.path(tab_dir, "bushfire_did_results.csv")
if (file.exists(bush_file)) {
  bush <- read.csv(bush_file, stringsAsFactors = FALSE)

  cat("\nBushfire DiD results:\n")
  for (i in seq_len(nrow(bush))) {
    r <- bush[i, ]
    if (!is.na(r$rr_high) && !is.na(r$rr_high_lo)) {
      ev <- compute_evalue(r$rr_high, r$rr_high_lo, r$rr_high_hi)

      evalue_results[[length(evalue_results) + 1]] <- data.frame(
        source      = "Bushfire DiD",
        outcome     = r$label,
        comparison  = "High smoke vs low",
        rr          = r$rr_high,
        rr_lo       = r$rr_high_lo,
        rr_hi       = r$rr_high_hi,
        e_point     = ev$e_point,
        e_ci        = ev$e_ci,
        stringsAsFactors = FALSE
      )

      cat("  ", r$label, ": RR =", r$rr_high,
          " -> E-value =", ev$e_point,
          "(CI:", ev$e_ci, ")\n")
    }
  }
}

# --- Key stratified results ---
# Respiratory at p99 (most extreme)
if (file.exists(dlnm_file)) {
  resp <- dlnm |> filter(group == "resp_total")
  if (nrow(resp) == 1 && !is.na(resp$rr_p99)) {
    ev99 <- compute_evalue(resp$rr_p99, resp$rr_p99_lo, resp$rr_p99_hi)
    evalue_results[[length(evalue_results) + 1]] <- data.frame(
      source      = "Primary DLNM",
      outcome     = "Respiratory",
      comparison  = "p99 vs median",
      rr          = resp$rr_p99,
      rr_lo       = resp$rr_p99_lo,
      rr_hi       = resp$rr_p99_hi,
      e_point     = ev99$e_point,
      e_ci        = ev99$e_ci,
      stringsAsFactors = FALSE
    )
    cat("\n  Respiratory p99: RR =", resp$rr_p99,
        " -> E-value =", ev99$e_point, "\n")
  }
}


# ==============================================================================
# 3. SAVE RESULTS
# ==============================================================================

evalue_table <- bind_rows(evalue_results)
write.csv(evalue_table,
          file.path(tab_dir, "evalues.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/evalues.csv\n")


# ==============================================================================
# 4. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 23 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("E-value summary:\n")
print(as.data.frame(
  evalue_table |>
    select(source, outcome, comparison, rr, e_point, e_ci)
), row.names = FALSE)

cat("\nInterpretation guide:\n")
cat("  E-value = minimum RR that an unmeasured confounder would need with\n")
cat("  BOTH exposure and outcome to explain away the observed association.\n")
cat("  E-value (CI) = same, but to move the CI to include null.\n")
cat("  Higher E-values = more robust to unmeasured confounding.\n")
cat("  E-value < 2: vulnerable to even modest confounders.\n")
cat("  E-value > 3: reasonably robust.\n")

cat("\nOutputs:\n")
cat("  outputs/tables/evalues.csv\n")
