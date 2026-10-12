#!/usr/bin/env Rscript

# Fresh authors-method fits for the paper's tensor-on-tensor simulation.
# This source-only script generates data locally and calls SegreCP's public API.
# The installed package does not contain competing-method implementations.

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript reproducibility/simulations/run_tensor_regression.R",
  "  --output-dir=/absolute/private/path [--full]",
  "  [--cell-id=1] [--replication=1] [--sweeps=2]",
  "",
  "Default: one paired design cell, one replication, two RGD/RGN sweeps.",
  "--full: all 12 cells x 20 replications x 10 sweeps (480 method fits).",
  "Smoke-only selectors cannot be combined with --full. A new output",
  "directory outside the package checkout is required; files are never",
  "overwritten. These are fresh SegreCP fits, not an exact replay of the",
  "historical seven-method runner or its timing measurements.",
  sep = "\n"
)
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
known <- args == "--full" |
  vapply(args, function(z) any(startsWith(z, c(
    "--output-dir=", "--cell-id=", "--replication=", "--sweeps="
  ))), logical(1L))
if (!length(args) || any(!known)) stop(usage, call. = FALSE)
argument <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Repeated argument: ", prefix, call. = FALSE)
  if (!length(hit)) return(default)
  substring(hit, nchar(prefix) + 1L)
}
whole_number <- function(x, label) {
  if (!grepl("^[1-9][0-9]*$", x)) {
    stop(label, " must be a positive integer.", call. = FALSE)
  }
  out <- suppressWarnings(as.integer(x))
  if (is.na(out)) stop(label, " is too large.", call. = FALSE)
  out
}
full <- "--full" %in% args
if (full && any(vapply(args, function(z) any(startsWith(
    z, c("--cell-id=", "--replication=", "--sweeps=")
  )), logical(1L)))) {
  stop("--full uses the fixed paper grid; omit smoke-only selectors.",
       call. = FALSE)
}
output_arg <- argument("--output-dir=")
if (is.null(output_arg) || !startsWith(output_arg, "/")) {
  stop("Pass an absolute --output-dir outside the SegreCP source tree.",
       call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                              mustWork = TRUE)

# Resolve an existing ancestor so a symlink cannot hide a package-local path.
canonical_target <- function(path) {
  pieces <- strsplit(path, "/", fixed = TRUE)[[1L]]
  stack <- character()
  for (piece in pieces) {
    if (!nzchar(piece) || piece == ".") next
    if (piece == "..") {
      if (length(stack)) stack <- head(stack, -1L)
    } else stack <- c(stack, piece)
  }
  path <- paste0("/", paste(stack, collapse = "/"))
  ancestor <- path
  suffix <- character()
  while (!file.exists(ancestor)) {
    suffix <- c(basename(ancestor), suffix)
    ancestor <- dirname(ancestor)
  }
  for (piece in suffix) {
    ancestor <- file.path(normalizePath(ancestor, mustWork = TRUE), piece)
  }
  ancestor
}
output_dir <- canonical_target(output_arg)
if (identical(output_dir, "/") || identical(output_dir, package_root) ||
    startsWith(output_dir, paste0(package_root, "/")) ||
    startsWith(package_root, paste0(output_dir, "/")) ||
    file.exists(output_dir)) {
  stop("Choose a new output directory outside the package checkout; ",
       "existing files are never overwritten.", call. = FALSE)
}
if (!requireNamespace("SegreCP", quietly = TRUE)) {
  stop("Install SegreCP before running this simulation.", call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install digest to record generated-data SHA-256 fingerprints.",
       call. = FALSE)
}

settings_path <- system.file(
  "extdata", "rebuttal_tensor_regression", "final_settings.csv",
  package = "SegreCP", mustWork = TRUE
)
settings <- utils::read.csv(settings_path, stringsAsFactors = FALSE)
if (!identical(names(settings), c("parameter", "value")) ||
    anyDuplicated(settings$parameter)) {
  stop("The installed simulation settings file is malformed.", call. = FALSE)
}
setting <- function(name) {
  z <- settings$value[match(name, settings$parameter)]
  if (length(z) != 1L || is.na(z)) {
    stop("Missing simulation setting: ", name, call. = FALSE)
  }
  z
}
parse_dims <- function(x) as.integer(strsplit(x, "x", fixed = TRUE)[[1L]])
output_dims <- parse_dims(setting("output_dimensions"))
input_dims <- parse_dims(setting("input_dimensions"))
dims_legacy <- c(output_dims, input_dims)
rank <- as.integer(setting("rank"))
kappa <- as.numeric(setting("kappa"))
n_train <- as.integer(setting("n_train"))
n_test <- as.integer(setting("n_test"))
coherence <- as.numeric(strsplit(setting("coherence"), ",", fixed = TRUE)[[1L]])
noise_sd <- as.numeric(strsplit(setting("noise_sd"), ",", fixed = TRUE)[[1L]])
replications <- as.integer(setting("replications"))
paper_sweeps <- as.integer(setting("iterations"))
master_seed <- as.integer(setting("master_seed"))
if (!identical(dims_legacy, rep(10L, 5L)) || rank != 3L || kappa != 30 ||
    n_train != 1000L || n_test != 1000L ||
    !identical(coherence, c(0, 0.2, 0.5, 0.75)) ||
    !identical(noise_sd, c(0, 0.5, 1)) || replications != 20L ||
    paper_sweeps != 10L || master_seed != 20260829L) {
  stop("The installed settings differ from this pinned simulation design.",
       call. = FALSE)
}

cell_grid <- expand.grid(coherence = coherence, noise_sd = noise_sd)
cell_grid$cell_id <- seq_len(nrow(cell_grid))
cell_grid$coherence_id <- match(cell_grid$coherence, coherence)
cell_grid$noise_id <- match(cell_grid$noise_sd, noise_sd)
if (full) {
  jobs <- merge(cell_grid, data.frame(replication = seq_len(replications)))
  jobs <- jobs[order(jobs$cell_id, jobs$replication), ]
  sweeps <- paper_sweeps
} else {
  cell_id <- whole_number(argument("--cell-id=", "1"), "--cell-id")
  replication <- whole_number(argument("--replication=", "1"), "--replication")
  sweeps <- whole_number(argument("--sweeps=", "2"), "--sweeps")
  if (cell_id > nrow(cell_grid) || replication > replications ||
      sweeps > paper_sweeps) {
    stop("Smoke selectors exceed the paper's 12 cells, 20 replications, ",
         "or 10 sweeps.", call. = FALSE)
  }
  jobs <- cell_grid[cell_id, , drop = FALSE]
  jobs$replication <- replication
}
rownames(jobs) <- NULL

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

# This is the small data generator from the benchmark, independent of its
# solver source. The five factors are ordered output-first during generation.
coherent_vectors <- function(p, r, rho) {
  gram <- stats::toeplitz(rho^(0:(r - 1L)))
  columns <- rbind(diag(1, r), matrix(0, p - r, r)) %*% chol(gram)
  columns <- apply(columns, 2L, function(z) z / sqrt(sum(z^2)))
  rotation <- qr.Q(qr(matrix(rnorm(p * p), nrow = p, ncol = p)))
  rotation %*% columns
}
make_truth <- function(rho) {
  factors <- lapply(dims_legacy, coherent_vectors, r = rank, rho = rho)
  lapply(seq_len(rank), function(i) {
    list(lambda = kappa^(0.5 - (i - 1L) / (rank - 1L)),
         U = lapply(factors, function(z) z[, i]))
  })
}
truth_array <- function(components) {
  value <- Reduce(`+`, lapply(components, function(z) {
    z$lambda * Reduce(kronecker, rev(z$U))
  }))
  array(value, dim = dims_legacy)
}

metric_rows <- vector("list", nrow(jobs) * 2L)
objective_rows <- vector("list", nrow(jobs) * 2L)
fingerprint_rows <- vector("list", nrow(jobs))
counter <- 0L
for (job_id in seq_len(nrow(jobs))) {
  row <- jobs[job_id, ]
  seed <- master_seed + row$coherence_id * 10000L + row$replication
  set.seed(seed)
  truth <- make_truth(row$coherence)
  truth_old_order <- truth_array(truth)
  # SegreCP documents input modes first; transpose the benchmark's
  # output-first tensor rather than changing the generated data.
  truth_public <- aperm(truth_old_order, c(3L, 4L, 5L, 1L, 2L))
  truth_norm <- sqrt(sum(truth_public^2))

  p_in <- prod(input_dims)
  q_out <- prod(output_dims)
  X_train <- matrix(rnorm(n_train * p_in, sd = 1 / sqrt(n_train)), n_train)
  X_test <- matrix(rnorm(n_test * p_in, sd = 1 / sqrt(n_train)), n_test)
  truth_matrix <- matrix(truth_old_order, nrow = q_out)
  Y_train_clean <- X_train %*% t(truth_matrix)
  Y_test_clean <- X_test %*% t(truth_matrix)
  Z_train <- matrix(rnorm(n_train * q_out), nrow = n_train)
  Z_test <- matrix(rnorm(n_test * q_out), nrow = n_test)
  Y_train <- Y_train_clean + row$noise_sd * Z_train / sqrt(n_train)
  Y_test <- Y_test_clean + row$noise_sd * Z_test / sqrt(n_train)

  fingerprint_rows[[job_id]] <- data.frame(
    cell_id = row$cell_id, replication = row$replication, seed = seed,
    coherence = row$coherence, noise_sd = row$noise_sd,
    data_sha256 = digest::digest(list(
      truth = truth, X_train = X_train, Y_train = Y_train,
      X_test = X_test, Y_test = Y_test
    ), algo = "sha256", serialize = TRUE),
    base_data_sha256 = digest::digest(list(
      truth = truth, X_train = X_train, Y_train_clean = Y_train_clean,
      X_test = X_test, Y_test_clean = Y_test_clean,
      Z_train = Z_train, Z_test = Z_test
    ), algo = "sha256", serialize = TRUE),
    stringsAsFactors = FALSE
  )

  # Historical CPCA groups the two output modes. In the public input-first
  # order these are modes 4 and 5. The same CPCA object starts RGD and RGN.
  moment <- array(as.numeric(crossprod(X_train, Y_train)),
                  dim = c(input_dims, output_dims))
  init_elapsed <- system.time({
    initial <- SegreCP::cp_initialize(
      moment, rank = rank, method = "cpca", split = c(4L, 5L)
    )
  })[["elapsed"]]

  for (method in c("RGD", "RGN")) {
    control <- SegreCP::segre_control(
      profile = "manuscript", max_iter = sweeps, tolerance = 0,
      step_size = 0.5,
      step_rule = if (method == "RGD") "exact_quadratic" else "fixed",
      component_order = "decreasing_weight", ridge = 0,
      normalize_regression = "none"
    )
    fit_elapsed <- system.time({
      fit <- SegreCP::segre_regress(
        X_train, Y_train, rank = rank,
        input_dims = input_dims, output_dims = output_dims,
        method = tolower(method), init = initial$comps,
        control = control
      )
    })[["elapsed"]]
    estimate <- stats::coef(fit)
    pred_train <- matrix(as.numeric(stats::predict(fit, X_train)), nrow = n_train)
    pred_test <- matrix(as.numeric(stats::predict(fit, X_test)), nrow = n_test)
    if (!all(is.finite(estimate)) || !all(is.finite(pred_test))) {
      stop("Nonfinite result in cell ", row$cell_id, ", replication ",
           row$replication, ", ", method, call. = FALSE)
    }
    counter <- counter + 1L
    metric_rows[[counter]] <- data.frame(
      mode = if (full) "full" else "smoke",
      cell_id = row$cell_id, coherence_id = row$coherence_id,
      noise_id = row$noise_id, replication = row$replication, seed = seed,
      coherence = row$coherence, noise_sd = row$noise_sd,
      method = method, rank = rank, sweeps = fit$iterations,
      coefficient_error = sqrt(sum((estimate - truth_public)^2)) / truth_norm,
      train_rmse = sqrt(mean((Y_train - pred_train)^2)),
      test_rmse = sqrt(mean((Y_test - pred_test)^2)),
      clean_test_rmse = sqrt(mean((Y_test_clean - pred_test)^2)),
      relative_prediction_error =
        sqrt(sum((Y_test_clean - pred_test)^2)) /
        max(sqrt(sum(Y_test_clean^2)), .Machine$double.eps),
      initialization_elapsed_s = unname(init_elapsed),
      fit_elapsed_s = unname(fit_elapsed),
      end_to_end_elapsed_s = unname(init_elapsed + fit_elapsed),
      stringsAsFactors = FALSE
    )
    trajectory <- SegreCP::convergence_data(fit)
    objective_rows[[counter]] <- data.frame(
      mode = if (full) "full" else "smoke",
      cell_id = row$cell_id, replication = row$replication, seed = seed,
      coherence = row$coherence, noise_sd = row$noise_sd,
      method = method, iteration = trajectory$iteration,
      objective = trajectory$objective,
      objective_ratio = trajectory$objective_ratio,
      stringsAsFactors = FALSE
    )
  }
  message(sprintf("Completed cell %d/%d, replication %d (%d/%d jobs).",
                  row$cell_id, nrow(cell_grid), row$replication,
                  job_id, nrow(jobs)))
}

metrics <- do.call(rbind, metric_rows)
objectives <- do.call(rbind, objective_rows)
fingerprints <- do.call(rbind, fingerprint_rows)
if (file.exists(output_dir)) {
  stop("Output directory appeared while fitting; refusing to overwrite it.",
       call. = FALSE)
}
if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Could not create output directory: ", output_dir, call. = FALSE)
}
utils::write.csv(metrics, file.path(output_dir, "authors_regression_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(objectives,
                 file.path(output_dir, "authors_regression_objective.csv"),
                 row.names = FALSE)
utils::write.csv(fingerprints,
                 file.path(output_dir, "generated_data_fingerprints.csv"),
                 row.names = FALSE)
utils::write.csv(data.frame(
  parameter = c("mode", "package_version", "rank", "input_dimensions",
                "output_dimensions", "n_train", "n_test", "replications_run",
                "sweeps", "master_seed", "initializer", "cpca_public_split",
                "RGD_step_rule", "component_order", "RGN_ridge",
                "normalization"),
  value = c(if (full) "full" else "smoke",
            as.character(utils::packageVersion("SegreCP")), rank,
            paste(input_dims, collapse = "x"),
            paste(output_dims, collapse = "x"), n_train, n_test,
            if (full) replications else 1L, sweeps, master_seed,
            "shared CPCA", "4,5", "exact_quadratic",
            "decreasing_weight", "0", "none")
), file.path(output_dir, "run_settings.csv"), row.names = FALSE)
utils::capture.output(utils::sessionInfo(),
                      file = file.path(output_dir, "sessionInfo.txt"))
message("Wrote fresh authors-method results to ", output_dir)
