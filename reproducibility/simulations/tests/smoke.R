#!/usr/bin/env Rscript

# Run both fresh authors-method simulations in new temporary directories.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) stop("Run this test with Rscript.")
this_file <- normalizePath(sub("^--file=", "", file_arg))
runner <- normalizePath(file.path(dirname(this_file), "..",
                                  "run_decomposition.R"))
output <- tempfile("segrecp-decomposition-smoke-")
status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(runner), paste0("--output-dir=", shQuote(output))))
if (status != 0L) stop("Fresh decomposition smoke run failed.")
path <- utils::read.csv(file.path(output, "trajectory.csv"))
terminal <- utils::read.csv(file.path(output, "terminal_accuracy.csv"))
checks <- utils::read.csv(file.path(output, "validation_checks.csv"))
comparison <- utils::read.csv(file.path(output, "reference_comparison.csv"))
stopifnot(nrow(path) == 16L, nrow(terminal) == 4L,
          all(checks$passed),
          nrow(comparison) == 1L, comparison$matched_rows == 16L,
          is.finite(comparison$cpca_max_absolute_error_difference),
          !is.na(comparison$checkpoint_rgd_rgn_rank_preserved),
          identical(sort(unique(path$method)), c("RGD", "RGN")),
          identical(sort(unique(path$initialization)), c("CPCA", "HOSVD")),
          identical(sort(unique(path$iteration)), 0:3))
cat("Fresh authors-only decomposition smoke test passed.\n")
cat("Output: ", output, "\n", sep = "")

regression_runner <- normalizePath(file.path(dirname(this_file), "..",
                                              "run_tensor_regression.R"))
regression_output <- tempfile("segrecp-regression-smoke-")
regression_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(shQuote(regression_runner),
    paste0("--output-dir=", shQuote(regression_output)))
)
if (regression_status != 0L) {
  stop("Fresh tensor-regression smoke run failed.")
}
metrics <- utils::read.csv(file.path(regression_output,
                                     "authors_regression_metrics.csv"))
objectives <- utils::read.csv(file.path(regression_output,
                                        "authors_regression_objective.csv"))
stopifnot(
  nrow(metrics) == 2L, nrow(objectives) == 6L,
  identical(sort(metrics$method), c("RGD", "RGN")),
  identical(sort(unique(objectives$iteration)), 0:2),
  all(is.finite(metrics$coefficient_error)),
  all(is.finite(objectives$objective))
)
cat("Fresh authors-only tensor-regression smoke test passed.\n")
cat("Output: ", regression_output, "\n", sep = "")
