#!/usr/bin/env Rscript

# Refit RGD and RGN on the processed 2018–2019 trade tensor.
# This source-only runner reads four local inputs and never downloads data.

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]),
                             mustWork = TRUE)
package_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                              mustWork = TRUE)
usage <- paste(
  "Usage: Rscript fit_authors_methods.R [--data-dir=PATH] [--output-dir=PATH]",
  "[--folds=1,2,3,4,5] [--repetitions=1] [--check-inputs] [--smoke]",
  "\n",
  "Refit the rank-three RGD/RGN trade analysis on the four prepared inputs.",
  "The data directory is outside the package source tree; --check-inputs",
  "checks the inputs without fitting. --smoke runs one",
  "fold, one repetition, and one iteration; it is not a paper result.",
  sep = "\n"
)
if (identical(args, "--help")) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
known <- c("--check-inputs", "--smoke")
values <- c("--data-dir=", "--output-dir=", "--folds=", "--repetitions=")
if (any(!vapply(args, function(arg) {
  arg %in% known || any(startsWith(arg, values))
}, logical(1L))) || anyDuplicated(args)) {
  stop(usage, call. = FALSE)
}
arg_value <- function(prefix, default = "") {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Repeated argument: ", prefix, call. = FALSE)
  if (!length(hit)) return(default)
  substring(hit, nchar(prefix) + 1L)
}
check_only <- "--check-inputs" %in% args
smoke <- "--smoke" %in% args
if (check_only && smoke) {
  stop("Use either --check-inputs or --smoke.", call. = FALSE)
}
if (smoke && any(vapply(args, function(arg) {
  any(startsWith(arg, c("--folds=", "--repetitions=")))
}, logical(1L)))) {
  stop("Smoke mode always runs fold 1 once for one iteration.", call. = FALSE)
}
data_dir <- arg_value("--data-dir=", Sys.getenv("SEGRECP_REAL_DATA_DIR"))
if (!nzchar(data_dir)) {
  stop("Set SEGRECP_REAL_DATA_DIR or pass --data-dir=PATH.", call. = FALSE)
}
data_dir <- normalizePath(data_dir, mustWork = TRUE)
output_dir <- arg_value("--output-dir=")
if (!check_only && !nzchar(output_dir)) {
  stop("Pass --output-dir=PATH for a fit run.", call. = FALSE)
}
if (!check_only) {
  # Resolve an existing ancestor so relative paths and symlinks cannot hide an
  # output directory inside the package source tree.
  target <- path.expand(output_dir)
  if (!startsWith(target, "/")) target <- file.path(getwd(), target)
  stack <- character()
  for (piece in strsplit(target, "/", fixed = TRUE)[[1L]]) {
    if (!nzchar(piece) || piece == ".") next
    if (piece == "..") stack <- head(stack, -1L) else stack <- c(stack, piece)
  }
  target <- paste0("/", paste(stack, collapse = "/"))
  ancestor <- target
  suffix <- character()
  while (!file.exists(ancestor)) {
    suffix <- c(basename(ancestor), suffix)
    ancestor <- dirname(ancestor)
  }
  output_dir <- if (length(suffix)) {
    do.call(file.path, as.list(c(normalizePath(ancestor, mustWork = TRUE), suffix)))
  } else normalizePath(ancestor, mustWork = TRUE)
  if (identical(output_dir, package_root) ||
      startsWith(output_dir, paste0(package_root, "/"))) {
    stop("Choose an output directory outside the SegreCP package source tree.",
         call. = FALSE)
  }
}

analysis_dir <- file.path(data_dir, "analysis_official")
inputs <- c(
  tensor = file.path(analysis_dir, "processed", "tensors_real.RData"),
  design = file.path(analysis_dir, "processed", "preprocessing_sensitivity_inputs",
                     "v3_scope_sqrt_raw_center_X_raw.csv"),
  folds = file.path(analysis_dir, "outputs", "tables",
                    "official_rank3_fold_assignments.csv"),
  registry = file.path(analysis_dir, "outputs", "tables",
                       "preprocessing_sensitivity_registry.csv")
)
missing <- names(inputs)[!file.exists(inputs)]
if (length(missing)) {
  stop("Missing processed input(s): ", paste(missing, collapse = ", "),
       call. = FALSE)
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required to check the prepared inputs.", call. = FALSE)
}

tensor_env <- new.env(parent = baseenv())
load(inputs[["tensor"]], envir = tensor_env)
required_objects <- c("Y", "exporters", "meta")
if (!all(vapply(required_objects, exists, logical(1L),
                envir = tensor_env, inherits = FALSE))) {
  stop("The tensor file lacks Y, exporters, or meta.", call. = FALSE)
}
Y <- tensor_env$Y
exporters <- tensor_env$exporters
meta <- tensor_env$meta
if (!is.numeric(Y) || !identical(as.integer(dim(Y)), c(80L, 25L, 21L, 8L)) ||
    any(!is.finite(Y))) {
  stop("Y must be a finite 80 x 25 x 21 x 8 numeric array.", call. = FALSE)
}
if (!is.data.frame(exporters) || !"exp_iso3" %in% names(exporters)) {
  stop("The exporter table has no exp_iso3 key.", call. = FALSE)
}
exporter_key <- as.character(exporters$exp_iso3)
if (length(exporter_key) != 80L || anyNA(exporter_key) ||
    any(!nzchar(exporter_key)) || anyDuplicated(exporter_key)) {
  stop("Exporter keys must be 80 distinct nonempty ISO3 codes.", call. = FALSE)
}
if (!is.list(meta) || !identical(meta$data_status, "official_real") ||
    !identical(as.integer(meta$selection_year), 2018L) ||
    !identical(as.integer(meta$covariate_year), 2018L) ||
    !identical(as.integer(meta$outcome_year), 2019L)) {
  stop("Tensor metadata do not describe the 2018-to-2019 official analysis.",
       call. = FALSE)
}

design <- read.csv(inputs[["design"]], check.names = FALSE)
features <- c("tariff", "ntb_proxy", "log_gdp_partner", "log_distance",
              "rta_scope_count", "world_demand")
if (!identical(names(design), c("exp_iso3", features)) ||
    nrow(design) != 80L || anyNA(design$exp_iso3) ||
    anyDuplicated(design$exp_iso3) ||
    !setequal(design$exp_iso3, exporter_key) ||
    !all(vapply(design[features], is.numeric, logical(1L)))) {
  stop("The trade design must contain one row per exporter and six numeric features ",
       "in the registered order.", call. = FALSE)
}
design <- design[match(exporter_key, design$exp_iso3), , drop = FALSE]
X <- as.matrix(design[features])
storage.mode(X) <- "double"
if (any(!is.finite(X))) stop("The trade design contains nonfinite values.", call. = FALSE)

folds <- read.csv(inputs[["folds"]], check.names = FALSE)
if (!all(c("exporter", "fold") %in% names(folds)) ||
    nrow(folds) != 80L || anyNA(folds$exporter) ||
    anyDuplicated(folds$exporter) ||
    !setequal(folds$exporter, exporter_key) ||
    !is.numeric(folds$fold) || anyNA(folds$fold) ||
    any(folds$fold != as.integer(folds$fold)) ||
    !identical(as.integer(table(factor(folds$fold, levels = 1:5))),
               rep(16L, 5L))) {
  stop("The fold file must assign each exporter once to five balanced folds.",
       call. = FALSE)
}
folds <- folds[match(exporter_key, folds$exporter), , drop = FALSE]
fold_id <- as.integer(folds$fold)

registry <- read.csv(inputs[["registry"]], check.names = FALSE)
required_registry <- c("variant_id", "rta_encoding", "aggregation_weight",
                       "tariff_transform", "x_scaling", "response_fit_scale",
                       "feature_count")
if (!all(required_registry %in% names(registry))) {
  stop("The preprocessing registry lacks required specification columns.",
       call. = FALSE)
}
spec <- registry[registry$variant_id == "v3_scope_sqrt_raw_center", , drop = FALSE]
expected <- c(rta_encoding = "scope_0_1_2", aggregation_weight = "sqrt_value",
              tariff_transform = "raw", x_scaling = "mean_sd",
              response_fit_scale = "center_only")
if (nrow(spec) != 1L || !identical(as.integer(spec$feature_count), 6L) ||
    !all(vapply(names(expected), function(name) {
      identical(as.character(spec[[name]]), expected[[name]])
    }, logical(1L)))) {
  stop("The registered trade preprocessing settings disagree with this runner.",
       call. = FALSE)
}

input_record <- data.frame(
  input = names(inputs),
  relative_path = substring(inputs, nchar(data_dir) + 2L),
  bytes = unname(file.info(inputs)$size),
  row.names = NULL
)
input_digest <- vapply(inputs, digest::digest, character(1L),
                       algo = "sha256", file = TRUE)
input_record$sha256 <- unname(input_digest)
reference_path <- file.path(dirname(script_path), "reference_derived_sha256.csv")
if (!file.exists(reference_path)) {
  stop("Missing reference for the four prepared inputs: ", reference_path,
       call. = FALSE)
}
reference <- read.csv(reference_path, check.names = FALSE)
if (!identical(names(reference), c("local_path", "bytes", "sha256")) ||
    nrow(reference) != 4L || anyDuplicated(reference$local_path) ||
    !setequal(reference$local_path, input_record$relative_path)) {
  stop("The prepared-input reference has an unexpected schema or path set.",
       call. = FALSE)
}
reference <- reference[match(input_record$relative_path, reference$local_path), ]
if (any(reference$bytes != input_record$bytes) ||
    any(tolower(reference$sha256) != unname(input_digest))) {
  stop("Prepared inputs differ from the recorded trade analysis.", call. = FALSE)
}
message("Prepared inputs match the recorded analysis: 80 exporters, five folds, ",
        "six features, and an 80 x 25 x 21 x 8 response.")
if (check_only) {
  print(input_record[c("input", "bytes")], row.names = FALSE)
  quit(save = "no", status = 0L)
}
prefix <- if (smoke) "authors_methods_smoke" else "authors_methods"
output_files <- stats::setNames(
  file.path(output_dir, paste0(prefix, "_",
                               c("per_fold.csv", "summary.csv", "objective.csv",
                                 "checkpoints.csv", "settings.rds", "inputs.csv",
                                 "session.txt"))),
  c("per_fold", "summary", "objective", "checkpoints", "settings", "inputs",
    "session")
)
existing_outputs <- output_files[file.exists(output_files)]
if (length(existing_outputs)) {
  stop("Result file(s) already exist: ",
       paste(basename(existing_outputs), collapse = ", "),
       ". Choose a new --output-dir to preserve the previous run.",
       call. = FALSE)
}
if (!requireNamespace("SegreCP", quietly = TRUE) ||
    utils::packageVersion("SegreCP") < numeric_version("0.1.5")) {
  stop("Install SegreCP 0.1.5 or later before fitting the authors' methods.",
       call. = FALSE)
}

# Standardize features and center responses within each training fold; the
# solver applies no additional regression normalization.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
        sample.kind = "Rejection")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
iterations <- if (smoke) 1L else 30L
checkpoints <- if (smoke) 1L else c(3L, 5L, 10L, 30L)
fold_text <- arg_value("--folds=", "1,2,3,4,5")
repetition_text <- arg_value("--repetitions=", "1")
if (!smoke && (!grepl("^[1-5](,[1-5])*$", fold_text) ||
               !grepl("^[1-9][0-9]*$", repetition_text))) {
  stop("Folds must be comma-separated integers 1 to 5; repetitions positive.",
       call. = FALSE)
}
folds_to_run <- if (smoke) 1L else as.integer(strsplit(fold_text, ",",
                                                    fixed = TRUE)[[1L]])
repetitions <- if (smoke) 1L else suppressWarnings(as.integer(repetition_text))
if (!length(folds_to_run) || anyNA(folds_to_run) ||
    any(!folds_to_run %in% 1:5) || anyDuplicated(folds_to_run) ||
    length(repetitions) != 1L || is.na(repetitions) ||
    repetitions < 1L || repetitions > 20L) {
  stop("Folds must be distinct integers 1 to 5; repetitions must be 1 to 20.",
       call. = FALSE)
}
control <- SegreCP::segre_control(
  profile = "practical", max_iter = iterations, tolerance = 0,
  step_size = 0.5, step_rule = "backtracking",
  component_order = "decreasing_weight", ridge = 0,
  normalize_regression = "none"
)

rows <- list()
paths <- list()
checkpoint_rows <- list()
fold_settings <- list()
position <- 0L
for (fold in folds_to_run) {
  test_id <- which(fold_id == fold)
  train_id <- which(fold_id != fold)
  x_mean <- colMeans(X[train_id, , drop = FALSE])
  x_sd <- pmax(apply(X[train_id, , drop = FALSE], 2L, sd), 1e-12)
  Xtr <- sweep(sweep(X[train_id, , drop = FALSE], 2L, x_mean, "-"),
               2L, x_sd, "/")
  Xte <- sweep(sweep(X[test_id, , drop = FALSE], 2L, x_mean, "-"),
               2L, x_sd, "/")
  y_mean <- apply(Y[train_id, , , , drop = FALSE], 4L, mean)
  Ytr <- sweep(Y[train_id, , , , drop = FALSE], 4L, y_mean, "-")
  Yte <- sweep(Y[test_id, , , , drop = FALSE], 4L, y_mean, "-")
  Ytr_matrix <- matrix(as.numeric(Ytr), nrow = length(train_id))
  Yte_matrix <- matrix(as.numeric(Yte), nrow = length(test_id))
  # CPCA uses X'Y/n; the solver retains the original training-data loss scale.
  initialization_seconds <- system.time({
    moment <- array(crossprod(Xtr / sqrt(nrow(Xtr)),
                              Ytr_matrix / sqrt(nrow(Xtr))),
                    dim = c(6L, 25L, 21L, 8L))
    initial <- SegreCP::cp_initialize(moment, rank = 3L, method = "cpca",
                                     split = c(1L, 2L), seed = 2549L)
  })[["elapsed"]]
  fold_settings[[as.character(fold)]] <- list(
    train_id = train_id, test_id = test_id, x_mean = x_mean, x_sd = x_sd,
    y_mean = y_mean, initialization = initial,
    initialization_seconds = as.numeric(initialization_seconds)
  )

  for (rep_id in seq_len(repetitions)) {
    # Rotate method order across repeats to reduce systematic timing bias.
    methods <- if ((fold + rep_id) %% 2L) c("RGD", "RGN") else c("RGN", "RGD")
    for (method in methods) {
    elapsed <- system.time({
      fit <- SegreCP::segre_regress(
        X = Xtr, Y = Ytr, rank = 3L, input_dims = 6L,
        output_dims = c(25L, 21L, 8L), method = tolower(method),
        init = initial$comps, design_variance = 1,
        control = control, keep_history = TRUE
      )
    })[["elapsed"]]
    if (fit$iterations != iterations) {
      stop(method, " did not complete the requested iterations in fold ", fold,
           call. = FALSE)
    }
    fitted_train <- matrix(as.numeric(predict(fit, Xtr)), nrow = length(train_id))
    fitted_test <- matrix(as.numeric(predict(fit, Xte)), nrow = length(test_id))
    if (any(!is.finite(fitted_train)) || any(!is.finite(fitted_test))) {
      stop(method, " produced nonfinite predictions in fold ", fold,
           call. = FALSE)
    }
    train_error <- Ytr_matrix - fitted_train
    test_error <- Yte_matrix - fitted_test
    position <- position + 1L
    rows[[position]] <- data.frame(
      fold = fold, repetition = rep_id, method = method, rank = 3L,
      iterations = fit$iterations,
      n_train = length(train_id), n_test = length(test_id),
      train_rmse = sqrt(mean(train_error^2)),
      test_rmse = sqrt(mean(test_error^2)),
      train_mae = mean(abs(train_error)), test_mae = mean(abs(test_error)),
      train_r2 = 1 - sum(train_error^2) / sum(Ytr_matrix^2),
      test_r2 = 1 - sum(test_error^2) / sum(Yte_matrix^2),
      fit_seconds = as.numeric(elapsed), package_version =
        as.character(utils::packageVersion("SegreCP"))
    )
    path <- SegreCP::convergence_data(fit)
    paths[[position]] <- data.frame(
      fold = fold, repetition = rep_id, method = method,
      path[, c("iteration", "objective", "objective_ratio")],
      check.names = FALSE
    )
    checkpoint_rows[[position]] <- do.call(rbind, lapply(checkpoints, function(k) {
      coefficient <- matrix(SegreCP::reconstruct_cp(
        fit$component_history[[k + 1L]]), nrow = ncol(Xtr))
      train_error <- Ytr_matrix - Xtr %*% coefficient
      test_error <- Yte_matrix - Xte %*% coefficient
      data.frame(
        fold = fold, repetition = rep_id, method = method, iteration = k,
        train_rmse = sqrt(mean(train_error^2)), test_rmse = sqrt(mean(test_error^2)),
        train_mae = mean(abs(train_error)), test_mae = mean(abs(test_error)),
        train_r2 = 1 - sum(train_error^2) / sum(Ytr_matrix^2),
        test_r2 = 1 - sum(test_error^2) / sum(Yte_matrix^2)
      )
    }))
    message(sprintf("fold %d repeat %d %s: test RMSE %.6f (%d iterations)",
                    fold, rep_id, method, rows[[position]]$test_rmse,
                    fit$iterations))
    }
  }
}

per_fold <- do.call(rbind, rows)
objective <- do.call(rbind, paths)
accuracy_columns <- c("train_rmse", "test_rmse", "train_mae", "test_mae",
                      "train_r2", "test_r2")
repeat_groups <- split(per_fold, interaction(per_fold$fold, per_fold$method),
                       drop = TRUE)
if (any(vapply(repeat_groups, function(x) {
  any(vapply(x[accuracy_columns], function(z) diff(range(z)) > 1e-10,
             logical(1L)))
}, logical(1L)))) {
  stop("Repeated fits gave different accuracy on the same fold.", call. = FALSE)
}
fold_accuracy <- per_fold[per_fold$repetition == 1L,
                          c("fold", "method", accuracy_columns)]
summary <- do.call(rbind, lapply(split(fold_accuracy, fold_accuracy$method), function(x) {
  fit_times <- per_fold$fit_seconds[per_fold$method == x$method[[1L]]]
  data.frame(
    method = x$method[[1L]], rank = 3L, folds = nrow(x),
    repetitions_per_fold = repetitions,
    mean_test_rmse = mean(x$test_rmse),
    se_test_rmse = if (nrow(x) > 1L) sd(x$test_rmse) / sqrt(nrow(x)) else NA_real_,
    mean_test_r2 = mean(x$test_r2),
    se_test_r2 = if (nrow(x) > 1L) sd(x$test_r2) / sqrt(nrow(x)) else NA_real_,
    mean_test_mae = mean(x$test_mae),
    mean_fit_seconds = mean(fit_times)
  )
}))
summary <- summary[order(summary$mean_test_rmse), , drop = FALSE]
summary$rmse_rank_among_authors_methods <- seq_len(nrow(summary))
rownames(summary) <- NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(per_fold, output_files[["per_fold"]], row.names = FALSE)
write.csv(summary, output_files[["summary"]], row.names = FALSE)
write.csv(objective, output_files[["objective"]], row.names = FALSE)
write.csv(do.call(rbind, checkpoint_rows), output_files[["checkpoints"]],
          row.names = FALSE)
saveRDS(list(
  control = control, checkpoints = checkpoints, seed = 2549L, RNGkind = RNGkind(),
  rank = 3L, input_dims = 6L, output_dims = c(25L, 21L, 8L),
  initializer = "CPCA", initializer_normalization = "sqrt_n",
  initializer_split = c(1L, 2L), design_variance = 1,
  metric_reference = "training-fold response-channel mean for train and test R2",
  timing_scope = "solver including history; excludes shared CPCA, preprocessing and metrics",
  package_version = as.character(utils::packageVersion("SegreCP")),
  installed_package = find.package("SegreCP"),
  runner_sha256 = digest::digest(script_path, algo = "sha256", file = TRUE),
  folds = fold_settings
), output_files[["settings"]])
write.csv(input_record, output_files[["inputs"]], row.names = FALSE)
writeLines(capture.output(sessionInfo()), output_files[["session"]])
message("Saved authors-method results under ", normalizePath(output_dir),
        if (smoke) " (smoke mode; not a paper result)." else ".")
