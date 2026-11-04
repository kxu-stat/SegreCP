#!/usr/bin/env Rscript
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
directory <- normalizePath(file.path(dirname(script), ".."))
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
