#!/usr/bin/env Rscript

# Exercise comparator routing from a detached, data-free SegreCP checkout.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this file with Rscript, not source().", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
comparison_source <- file.path(dirname(script_path), "..", "comparison_methods.R")
package_source <- normalizePath(file.path(dirname(script_path), "..", ".."),
                                mustWork = TRUE)
invisible(parse(file = comparison_source))

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

fixture <- tempfile(pattern = "segrecp-standalone-", tmpdir = tempdir())
dir.create(file.path(fixture, "reproducibility"), recursive = TRUE)
check(file.copy(file.path(package_source, "DESCRIPTION"), fixture),
      "Could not stage package metadata in the temporary checkout")
fixture_script <- file.path(fixture, "reproducibility", "comparison_methods.R")
check(file.copy(comparison_source, fixture_script),
      "Could not stage the comparison script in the temporary checkout")

old_wd <- getwd()
setwd(fixture)
comparison <- new.env(parent = globalenv())
source(fixture_script, local = comparison)
check(identical(comparison$comparison_repository_root(),
                normalizePath(fixture, mustWork = TRUE)),
      "The comparator must locate a standalone SegreCP source checkout")
catalog <- comparison$comparison_method_catalog()
check(identical(catalog$repository_worker[[1L]],
                "reproducibility/comparison_methods.R") &&
        identical(catalog$required_checkout[[1L]], "SegreCP source checkout") &&
        all(catalog$required_checkout[-1L] == "complete manuscript repository"),
      "The method catalog must distinguish standalone CP-ALS from external workers")

missing_worker <- tryCatch(
  comparison$run_tensor_regression_comparisons(execute = FALSE),
  error = conditionMessage
)
check(is.character(missing_worker) && length(missing_worker) == 1L &&
        grepl("requires the complete manuscript repository", missing_worker,
              fixed = TRUE) &&
        grepl("rgn_benchmark/run_tensor_on_tensor_baseline_grid.R",
              missing_worker, fixed = TRUE),
      "A detached checkout must report the missing manuscript worker clearly")
setwd(old_wd)

if (requireNamespace("SegreCP", quietly = TRUE) &&
    requireNamespace("rTensor", quietly = TRUE)) {
  set.seed(11L)
  observed <- array(stats::rnorm(3L * 4L * 2L), c(3L, 4L, 2L))
  initial <- SegreCP::cp_initialize(observed, rank = 1L,
                                    method = "hosvd")$comps
  fit <- comparison$fit_comparison_cp_als(observed, initial, iterations = 2L)
  check(identical(fit$method, "CP-ALS") &&
          nrow(fit$trajectory) == 3L &&
          all(is.finite(fit$trajectory$objective)),
        "Standalone CP-ALS must complete two finite synthetic sweeps")
  cat("PASS: detached checkout routing and synthetic CP-ALS fit.\n")
} else {
  cat("PASS: detached checkout routing. SKIP synthetic CP-ALS fit: ",
      "install SegreCP and rTensor first.\n", sep = "")
}
