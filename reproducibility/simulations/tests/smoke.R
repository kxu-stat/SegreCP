#!/usr/bin/env Rscript
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
directory <- normalizePath(file.path(dirname(script), ".."))
source(file.path(directory, "helpers.R"))

# Compare contracted diagnostic predictions with independent dense products.
set.seed(104L)
metric_checks <- 0L
cases <- list(list(input = c(3L, 2L), output = c(2L, 3L)),
              list(input = c(2L, 1L, 3L), output = 4L),
              list(input = 3L, output = c(1L, 2L, 3L)),
              list(input = c(2L, 3L), output = integer()))
for (case in cases) for (rank in c(1L, 3L)) for (sigma in c(0, .2)) {
  dims <- c(case$input, case$output)
  components <- lapply(seq_len(rank), function(i) list(lambda = i / rank,
    U = lapply(dims, function(p) { u <- rnorm(p); u / sqrt(sum(u^2)) })))
  coefficient <- simulation_reconstruct(components)
  B <- matrix(as.numeric(coefficient), prod(case$input), prod(case$output))
  X_train <- matrix(rnorm(9L * nrow(B)), 9L)
  X_test <- matrix(rnorm(7L * nrow(B)), 7L)
  train_prediction <- X_train %*% B
  test_prediction <- X_test %*% B
  data <- list(X_train = X_train, X_test = X_test, truth = coefficient,
    input_dims = case$input, output_dims = case$output,
    Y_train = train_prediction + sigma * matrix(rnorm(length(train_prediction)), 9L),
    Y_test = test_prediction + sigma * matrix(rnorm(length(test_prediction)), 7L),
    Y_test_clean = test_prediction)
  actual <- unlist(simulation_metrics(components, data, "regression"))
  train <- train_prediction - data$Y_train
  test <- test_prediction - data$Y_test
  expected <- c(relative_recovery_error = 0, objective = .5 * sum(train^2),
    train_RMSE = sqrt(mean(train^2)), test_RMSE = sqrt(mean(test^2)),
    clean_test_RMSE = 0, relative_clean_prediction_error = 0)
  stopifnot(identical(names(actual), names(expected)),
            all(abs(actual - expected) <= 1e-12 * pmax(1, abs(expected))))
  metric_checks <- metric_checks + 1L
}
cat("PASS: ", metric_checks, " contracted-metric cases match dense predictions at 1e-12 tolerance.\n", sep = "")

for (task in c("decomposition", "regression")) {
  output <- tempfile(paste0("segrecp-", task, "-signal4-smoke-"))
  runner <- file.path(directory, if (task == "decomposition") "run_decomposition.R" else "run_tensor_regression.R")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(runner), paste0("--output-dir=", shQuote(output))))
  stopifnot(status == 0L)
  settings <- readRDS(file.path(output, "settings.rds"))
  paths <- utils::read.csv(file.path(output, "trajectories.csv"))
  fits <- utils::read.csv(file.path(output, "fits.csv"))
  inputs <- utils::read.csv(file.path(output, "inputs.csv"))
  stopifnot(settings$study == "signal4_formal20", settings$signal_multiplier == 4,
    settings$ridge == 0, settings$master_seed == 741103000L,
    settings$component_order == "decreasing_weight", settings$rank == 3L,
    settings$coherence == .5, settings$kappa == 30,
    nrow(paths) == 4L * (settings$iterations + 1L), nrow(fits) == 4L,
    all(fits$ridge == 0), all(fits$actual_iterations == settings$iterations),
    all(fits$engine == settings$engine), all(fits$stopping_reason == "max_iter"),
    identical(sort(unique(paths$noise_sd)), c(0L, 1L)),
    identical(sort(unique(paths$method)), c("RGD", "RGN")),
    all(is.finite(paths$relative_recovery_error)), all(is.finite(paths$objective)),
    all(nchar(inputs$sha256) == 64L))
  for (path in list.files(file.path(output, "fits"), full.names = TRUE)) {
    fit <- readRDS(path)
    stopifnot(length(fit$component_history) == settings$iterations + 1L,
              fit$control$ridge == 0)
  }
  if (task == "regression") stopifnot(all(is.finite(paths$train_RMSE)),
    all(is.finite(paths$test_RMSE)), all(is.finite(paths$clean_test_RMSE)))
  cat("PASS: ", task, " signal4 paired smoke; results: ", output, "\n", sep = "")
}
