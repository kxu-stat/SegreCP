#!/usr/bin/env Rscript

# Deterministically rebuild the rebuttal-stage tensor-regression tables from
# the paired replicate-level records shipped with SegreCP.

args <- commandArgs(trailingOnly = TRUE)
output_position <- match("--output", args)
if (is.na(output_position)) {
  output_dir <- file.path(getwd(), "SegreCP_rebuttal_table_rebuild")
  consumed <- integer()
} else {
  if (output_position == length(args)) {
    stop("Missing value after --output.", call. = FALSE)
  }
  output_dir <- path.expand(args[[output_position + 1L]])
  consumed <- c(output_position, output_position + 1L)
}
if (length(setdiff(seq_along(args), consumed))) {
  stop("Unknown command-line argument(s): ",
       paste(args[setdiff(seq_along(args), consumed)], collapse = ", "),
       call. = FALSE)
}

if (!requireNamespace("SegreCP", quietly = TRUE)) {
  stop("Install SegreCP before running this validator.", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install the suggested package 'digest' for provenance hashes.",
       call. = FALSE)
}
data_dir <- system.file(
  "extdata", "rebuttal_tensor_regression", package = "SegreCP"
)
if (!nzchar(data_dir)) {
  stop("The installed package does not contain the rebuttal data bundle.",
       call. = FALSE)
}
if (dir.exists(output_dir) && length(list.files(
  output_dir, all.files = TRUE, no.. = TRUE
))) {
  stop("Refusing to overwrite nonempty output directory: ", output_dir,
       call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_bundle <- function(file) {
  utils::read.csv(file.path(data_dir, file), stringsAsFactors = FALSE,
                  check.names = FALSE)
}

grid_raw <- read_bundle("final_seven_method_grid_raw.csv")
checkpoint_raw <- read_bundle("checkpoint_metrics_raw.csv")
balanced_raw <- read_bundle("final_seven_method_replication_balanced.csv")
timing_raw <- read_bundle("final_seven_method_serial_timing_raw.csv")
stored_checks <- read_bundle("final_validation_checks.csv")
checkpoint_manifest <- read_bundle("checkpoint_source_manifest.csv")
bundle_manifest <- read_bundle("package_subset_manifest.csv")

mean_se <- function(data, metric, groups) {
  formula <- stats::as.formula(paste(metric, "~", paste(groups, collapse = "+")))
  means <- stats::aggregate(formula, data, mean)
  ses <- stats::aggregate(formula, data, function(x) stats::sd(x) / sqrt(length(x)))
  names(means)[[ncol(means)]] <- "mean"
  names(ses)[[ncol(ses)]] <- "se"
  merge(means, ses, by = groups, sort = FALSE)
}

panel_a <- mean_se(
  checkpoint_raw, "coefficient_error", c("method", "iteration")
)
panel_b <- mean_se(
  grid_raw[grid_raw$coherence == 0.75, ],
  "clean_test_rmse", c("display_method", "noise_sd")
)
panel_c_coefficient <- mean_se(
  balanced_raw, "coefficient_error", "display_method"
)
panel_c_clean <- mean_se(
  balanced_raw, "clean_test_rmse", "display_method"
)
panel_c_noisy <- mean_se(
  balanced_raw, "test_rmse", "display_method"
)
panel_c_time <- mean_se(
  timing_raw, "reported_fit_sec", "display_method"
)
panel_c <- Reduce(function(x, y) merge(
  x, y, by = "display_method", sort = FALSE
), list(
  transform(panel_c_coefficient,
            coefficient_error = mean, coefficient_error_se = se)[
              c("display_method", "coefficient_error", "coefficient_error_se")
            ],
  transform(panel_c_clean,
            clean_test_rmse = mean, clean_test_rmse_se = se)[
              c("display_method", "clean_test_rmse", "clean_test_rmse_se")
            ],
  transform(panel_c_noisy,
            noisy_test_rmse = mean, noisy_test_rmse_se = se)[
              c("display_method", "noisy_test_rmse", "noisy_test_rmse_se")
            ],
  transform(panel_c_time,
            fit_sec = mean, fit_sec_se = se)[
              c("display_method", "fit_sec", "fit_sec_se")
            ]
))

methods <- c("CP-ALS", "RRR", "Tucker-RGD", "HOPLS", "TT-IHT", "RGD", "RGN")
expected_grid <- expand.grid(
  coherence = c(0, 0.2, 0.5, 0.75),
  noise_sd = c(0, 0.5, 1)
)
pair_count <- table(interaction(
  grid_raw$coherence, grid_raw$noise_sd, grid_raw$replication, drop = TRUE
))
seed_count <- stats::aggregate(
  seed ~ coherence + noise_sd + replication, grid_raw,
  function(x) length(unique(x))
)
finite_columns <- c(
  "coefficient_error", "train_rmse", "test_rmse", "clean_test_rmse",
  "relative_prediction_error", "elapsed_s"
)

time_mean <- setNames(panel_c$fit_sec, panel_c$display_method)
runtime_tiers_hold <-
  max(time_mean[c("RGD", "RGN")]) < time_mean[["CP-ALS"]] &&
  time_mean[["CP-ALS"]] < time_mean[["HOPLS"]] &&
  time_mean[["HOPLS"]] < time_mean[["TT-IHT"]] &&
  time_mean[["TT-IHT"]] < min(time_mean[c("RRR", "Tucker-RGD")])

derived_checkpoint <- checkpoint_manifest[
  checkpoint_manifest$role == "derived_checkpoint_archive",
]
checkpoint_hash_matches <- nrow(derived_checkpoint) == 1L &&
  identical(
    digest::digest(
      file.path(data_dir, "checkpoint_metrics_raw.csv"),
      algo = "sha256", file = TRUE, serialize = FALSE
    ),
    derived_checkpoint$sha256[[1L]]
  )
bundle_paths <- file.path(data_dir, bundle_manifest$file)
bundle_hashes_match <- all(file.exists(bundle_paths)) && identical(
  unname(vapply(
    bundle_paths, digest::digest, character(1L), algo = "sha256",
    file = TRUE, serialize = FALSE
  )),
  bundle_manifest$sha256
)

checks <- data.frame(
  check = c(
    "stored validation checks pass",
    "exact 1680-row paired grid",
    "exact seven-method set",
    "exact coherence-noise grid",
    "seven methods in every cell-replication pair",
    "one seed in every cell-replication pair",
    "all reported metrics finite",
    "complete 360-row checkpoint archive",
    "checkpoint archive SHA-256 matches its provenance manifest",
    "all installed rebuttal-bundle SHA-256 hashes match",
    "runtime tier ordering preserved"
  ),
  passed = c(
    all(stored_checks$passed),
    nrow(grid_raw) == 1680L,
    setequal(unique(grid_raw$display_method), methods),
    setequal(
      interaction(grid_raw$coherence, grid_raw$noise_sd, drop = TRUE),
      interaction(expected_grid$coherence, expected_grid$noise_sd, drop = TRUE)
    ),
    all(pair_count == 7L),
    all(seed_count$seed == 1L),
    all(vapply(grid_raw[finite_columns], function(x) all(is.finite(x)),
               logical(1L))),
    nrow(checkpoint_raw) == 360L &&
      all(table(checkpoint_raw$method, checkpoint_raw$iteration) == 20L),
    checkpoint_hash_matches,
    bundle_hashes_match,
    runtime_tiers_hold
  ),
  stringsAsFactors = FALSE
)
if (!all(checks$passed)) {
  print(checks, row.names = FALSE)
  stop("Rebuttal table validation failed.", call. = FALSE)
}

utils::write.csv(panel_a, file.path(output_dir, "panel_a_checkpoints.csv"),
                 row.names = FALSE)
utils::write.csv(panel_b, file.path(output_dir, "panel_b_high_coherence.csv"),
                 row.names = FALSE)
utils::write.csv(panel_c, file.path(output_dir, "panel_c_accuracy_timing.csv"),
                 row.names = FALSE)
utils::write.csv(checks, file.path(output_dir, "validation_checks.csv"),
                 row.names = FALSE)
writeLines(c(
  "SegreCP rebuttal-stage tensor-regression table rebuild",
  sprintf("Passed: %d/%d checks", sum(checks$passed), nrow(checks)),
  capture.output(print(checks, row.names = FALSE)),
  "",
  "Wall-clock seconds are reference-machine measurements. Fresh runs should",
  "validate the separated runtime tiers, not exact seconds or within-tier order."
), file.path(output_dir, "VALIDATION_REPORT.txt"))

cat("\nRecomputed Panel A:\n")
print(panel_a, row.names = FALSE, digits = 6L)
cat("\nRecomputed Panel B:\n")
print(panel_b, row.names = FALSE, digits = 6L)
cat("\nRecomputed Panel C:\n")
print(panel_c, row.names = FALSE, digits = 6L)
cat("\nValidation:\n")
print(checks, row.names = FALSE)
cat("\nSaved deterministic table rebuild to:\n",
    normalizePath(output_dir), "\n")
