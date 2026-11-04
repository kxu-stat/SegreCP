#!/usr/bin/env Rscript

# Optional regression comparisons; solver adapters stay in the manuscript repo.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(),
  value = TRUE)[[1L]]))
source(file.path(dirname(script), "helpers.R"))
source(file.path(dirname(script), "..", "comparison_methods.R"))
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  cat("Usage: Rscript run_comparisons.R --output-dir=PATH [--full]\n",
      "  [--methods=CP-ALS,RRR,Tucker-RGD,TT-IHT,HOPLS]\n",
      "Default: CP-ALS and supplementary RRR, first dimension/replication,\n",
      "two iterations at both noise levels. --full uses all dimensions,\n",
      "20 replications and ten iterations. Requires the manuscript workspace.\n",
      sep = "")
  quit(status = 0L)
}
option <- function(prefix, default = NULL) {
  values <- args[startsWith(args, prefix)]
  if (length(values) > 1L) stop("Repeated option: ", prefix)
  if (length(values)) substring(values, nchar(prefix) + 1L) else default
}
if (any(!grepl("^--(full$|output-dir=|methods=)", args))) stop("Unknown option")
output <- option("--output-dir=")
if (is.null(output) || !nzchar(output) || file.exists(output))
  stop("Supply a new --output-dir")
methods <- strsplit(option("--methods=", "CP-ALS,RRR"), ",", fixed = TRUE)[[1L]]
stopifnot(length(methods) > 0L, !anyDuplicated(methods),
  all(methods %in% c("CP-ALS", "RRR", "Tucker-RGD", "TT-IHT", "HOPLS")))
root <- .comparison_manuscript_root()
if (!requireNamespace("SegreCP", quietly = TRUE) ||
    utils::packageVersion("SegreCP") < "0.1.5") stop("Install SegreCP >= 0.1.5.")
settings <- simulation_settings("regression", full = "--full" %in% args)
dir.create(output, recursive = TRUE)
settings$methods <- methods
settings$engine <- "experiments/current_package/regression_adapters.R"
settings$adapter <- file.path(root, "experiments/current_package/regression_adapters.R")
settings$adapter_sha256 <- digest::digest(settings$adapter, file = TRUE, algo = "sha256")
saveRDS(settings, file.path(output, "settings.rds"))
metrics <- objectives <- controls <- list()
for (p_index in settings$p_indices) for (replication in settings$replications) {
  P <- settings$P[[p_index]]
  seed <- settings$master_seed + settings$task_offset + 1000L * p_index + replication
  base <- simulation_regression_base(P, seed)
  for (sigma in settings$noise_sd) {
    data <- simulation_regression_case(base, sigma)
    initial <- simulation_regression_initialization(data, seed)
    for (method in methods) {
      key <- paste(P, replication, sigma, method, sep = "_")
      message(key)
      fit <- fit_comparison_tensor_regression(method, data$X_train, data$Y_train,
        c(P, 10L, 10L), c(10L, 10L), settings$rank, initial,
        iterations = settings$iterations, seed = seed)
      coefficient <- matrix(fit$coefficient, ncol(data$X_train))
      train <- data$X_train %*% coefficient
      test <- data$X_test %*% coefficient
      metrics[[key]] <- data.frame(P = P, replication = replication, seed = seed,
        noise_sd = sigma, method = method, requested_iterations = settings$iterations,
        actual_iterations = fit$actual_iterations,
        latent_blocks = if (method == "HOPLS") fit$metadata$latent_blocks else NA_integer_,
        relative_recovery_error = sqrt(sum((fit$coefficient - data$truth)^2) / sum(data$truth^2)),
        train_RMSE = sqrt(mean((train - matrix(data$Y_train, nrow(train)))^2)),
        test_RMSE = sqrt(mean((test - matrix(data$Y_test, nrow(test)))^2)),
        ridge = fit$metadata$regularization, stop_reason = fit$stop_reason)
      objectives[[key]] <- data.frame(P = P, replication = replication,
        noise_sd = sigma, method = method, fit$objective)
      controls[[key]] <- fit$metadata
      stopifnot(fit$finite, fit$metadata$regularization == 0)
      write.csv(do.call(rbind, metrics), file.path(output, "metrics.csv"), row.names = FALSE)
      write.csv(do.call(rbind, objectives), file.path(output, "objectives.csv"), row.names = FALSE)
      saveRDS(controls, file.path(output, "controls.rds"))
    }
  }
}
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
