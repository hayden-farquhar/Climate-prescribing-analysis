# ==============================================================================
# Script 22: Multiple Testing Adjustment
# ==============================================================================
# Applies Benjamini-Hochberg (FDR) correction across all primary hypothesis
# tests to account for the number of models fitted.
#
# Collects p-values from:
#   - Primary DLNM models (7 groups, script 05)
#   - Equity interaction tests (3 groups, script 06)
#   - Remoteness interaction tests (3 groups, script 09)
#   - Temporal interaction tests (3 groups, script 13)
#   - Compound interaction tests (3 groups, script 14)
#
# Inputs:
#   outputs/tables/dlnm_results_summary.csv
#   outputs/tables/equity_interaction_tests.csv
#   outputs/tables/remoteness_interaction_tests.csv
#   outputs/tables/temporal_interaction_tests.csv
#   outputs/tables/compound_interaction_tests.csv
#
# Outputs:
#   outputs/tables/multiple_testing_adjusted.csv
#
# Dependencies:
#   install.packages(c("dplyr"))
# ==============================================================================

library(dplyr)
library(dlnm)
library(splines)
library(lubridate)

project_dir <- here::here()
data_dir    <- file.path(project_dir, "data", "processed")
tab_dir     <- file.path(project_dir, "outputs", "tables")

cat("=" |> strrep(70), "\n")
cat("Script 22: Multiple Testing Adjustment\n")
cat("=" |> strrep(70), "\n\n")


# ==============================================================================
# 1. COLLECT P-VALUES FROM PRIMARY DLNM MODELS
# ==============================================================================

# The primary DLNM results table has RR and CI but not explicit p-values.
# Compute p-values from the Wald test: z = log(RR) / SE, where
# SE = (log(RR_hi) - log(RR_lo)) / (2 * 1.96)

cat("Collecting p-values from primary DLNM models...\n")

all_pvals <- list()

dlnm_file <- file.path(tab_dir, "dlnm_results_summary.csv")
if (file.exists(dlnm_file)) {
  dlnm_results <- read.csv(dlnm_file, stringsAsFactors = FALSE)

  for (i in seq_len(nrow(dlnm_results))) {
    r <- dlnm_results[i, ]
    log_rr <- log(r$rr_p95)
    se <- (log(r$rr_p95_hi) - log(r$rr_p95_lo)) / (2 * 1.96)
    z <- log_rr / se
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)

    all_pvals[[length(all_pvals) + 1]] <- data.frame(
      source     = "Primary DLNM (script 05)",
      test       = paste0(r$label, " — RR at p95"),
      estimate   = round(r$rr_p95, 4),
      ci_lo      = round(r$rr_p95_lo, 4),
      ci_hi      = round(r$rr_p95_hi, 4),
      p_value    = p,
      stringsAsFactors = FALSE
    )
  }
  cat("  Primary DLNM:", nrow(dlnm_results), "tests\n")
}


# ==============================================================================
# 2. COLLECT P-VALUES FROM INTERACTION TESTS
# ==============================================================================

cat("Collecting p-values from interaction tests...\n")

interaction_files <- list(
  "Equity interaction (script 06)" = "equity_interaction_tests.csv",
  "Remoteness interaction (script 09)" = "remoteness_interaction_tests.csv",
  "Temporal interaction (script 13)" = "temporal_interaction_tests.csv",
  "Compound interaction (script 14)" = "compound_interaction_tests.csv"
)

for (source_name in names(interaction_files)) {
  fpath <- file.path(tab_dir, interaction_files[[source_name]])
  if (!file.exists(fpath)) {
    cat("  ", source_name, ": file not found, skipping\n")
    next
  }

  int_data <- read.csv(fpath, stringsAsFactors = FALSE)

  for (i in seq_len(nrow(int_data))) {
    r <- int_data[i, ]
    all_pvals[[length(all_pvals) + 1]] <- data.frame(
      source     = source_name,
      test       = paste0(r$label, " — interaction F-test"),
      estimate   = round(r$f_stat, 3),
      ci_lo      = NA,
      ci_hi      = NA,
      p_value    = r$p_value,
      stringsAsFactors = FALSE
    )
  }
  cat("  ", source_name, ":", nrow(int_data), "tests\n")
}


# ==============================================================================
# 3. APPLY BH CORRECTION
# ==============================================================================

pval_table <- bind_rows(all_pvals)
n_tests <- nrow(pval_table)

cat("\nTotal tests collected:", n_tests, "\n")

# Benjamini-Hochberg FDR correction
pval_table$p_bh <- p.adjust(pval_table$p_value, method = "BH")

# Bonferroni for comparison (more conservative)
pval_table$p_bonf <- p.adjust(pval_table$p_value, method = "bonferroni")

# Significance flags
pval_table$sig_nominal  <- pval_table$p_value < 0.05
pval_table$sig_bh       <- pval_table$p_bh < 0.05
pval_table$sig_bonf     <- pval_table$p_bonf < 0.05

cat("\nSignificant at nominal p<0.05:", sum(pval_table$sig_nominal), "/", n_tests, "\n")
cat("Significant after BH correction:", sum(pval_table$sig_bh), "/", n_tests, "\n")
cat("Significant after Bonferroni:", sum(pval_table$sig_bonf), "/", n_tests, "\n")


# ==============================================================================
# 4. SAVE RESULTS
# ==============================================================================

# Round p-values for readability
pval_table$p_value <- signif(pval_table$p_value, 4)
pval_table$p_bh    <- signif(pval_table$p_bh, 4)
pval_table$p_bonf  <- signif(pval_table$p_bonf, 4)

write.csv(pval_table,
          file.path(tab_dir, "multiple_testing_adjusted.csv"),
          row.names = FALSE)
cat("\n-> Saved: outputs/tables/multiple_testing_adjusted.csv\n")


# ==============================================================================
# 5. SUMMARY
# ==============================================================================

cat("\n", "=" |> strrep(70), "\n")
cat("Script 22 complete\n")
cat("=" |> strrep(70), "\n\n")

cat("Multiple testing adjustment results:\n\n")
print(as.data.frame(
  pval_table |>
    select(source, test, p_value, p_bh, sig_bh, p_bonf, sig_bonf)
), row.names = FALSE)

cat("\nTests that survive BH correction:\n")
surviving <- pval_table |> filter(sig_bh)
if (nrow(surviving) > 0) {
  print(as.data.frame(
    surviving |> select(source, test, p_value, p_bh)
  ), row.names = FALSE)
} else {
  cat("  None\n")
}

cat("\nOutputs:\n")
cat("  outputs/tables/multiple_testing_adjusted.csv\n")
