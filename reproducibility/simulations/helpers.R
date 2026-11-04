# Fixed signal4_formal20 design. Sourcing this file does not run an experiment.
simulation_settings <- function(task, full = FALSE, p_index = 1L,
                                replication = 1L, iterations = NULL) {
  if (!requireNamespace("SegreCP", quietly = TRUE) ||
      utils::packageVersion("SegreCP") < "0.1.5") stop("Install SegreCP >= 0.1.5.")
  task <- match.arg(task, c("decomposition", "regression"))
  P <- if (task == "decomposition") c(30L, 60L, 90L, 120L) else c(40L, 60L, 100L)
  iterations <- if (full) 10L else if (!is.null(iterations)) iterations else
    if (task == "decomposition") 3L else 2L
  stopifnot(p_index %in% seq_along(P), replication %in% 1:20,
            iterations %in% 1:10)
  list(study = "signal4_formal20", task = task,
       mode = if (full) "full" else "smoke", P = P,
       p_indices = if (full) seq_along(P) else p_index,
       replications = if (full) 1:20 else replication,
       iterations = as.integer(iterations), rank = 3L, kappa = 30,
       coherence = .5, signal_multiplier = 4, noise_sd = c(0, 1),
       master_seed = 741103000L,
       task_offset = if (task == "decomposition") 0L else 100000L,
       cpca_split = if (task == "decomposition") c(1L, 2L) else c(4L, 5L),
       n_train = 1000L, n_test = 1000L, output_dims = c(10L, 10L),
       ridge = 0, component_order = "decreasing_weight",
       normalize_regression = "none", methods = c("RGD", "RGN"),
       engine = if (task == "decomposition") "SegreCP::segre_decompose" else
         "SegreCP::segre_regress",
       package_version = as.character(utils::packageVersion("SegreCP")),
       package_path = normalizePath(find.package("SegreCP")))
}

simulation_hash <- function(x) {
  digest::digest(x, algo = "sha256", serialize = TRUE, serializeVersion = 3L)
}

simulation_reconstruct <- function(components) {
  dims <- vapply(components[[1L]]$U, length, integer(1L))
  result <- array(0, dims)
  for (component in components) {
    product <- component$U[[1L]]
    for (mode in seq.int(2L, length(dims)))
      product <- outer(product, component$U[[mode]])
    result <- result + component$lambda * product
  }
  result
}

simulation_decomposition_base <- function(P, seed) {
  set.seed(seed)
  dims <- c(P, 30L, 30L, 30L)
  factors <- lapply(dims, function(p) {
    gaussian <- matrix(stats::rnorm(p * 3L), p, 3L)
    qr.Q(qr(gaussian)) %*% chol(toeplitz(.5^(0:2)))
  })
  truth <- lapply(seq_len(3L), function(i) {
    list(lambda = 1 * sqrt(3L) * 30^(.5 - (i - 1L) / (3L - 1L)) *
           max(dims)^(length(dims) / 4), U = lapply(factors, function(A) A[, i]))
  })
  signal <- SegreCP::reconstruct_cp(truth)
  noise <- array(stats::rnorm(prod(dims), sd = 1), dims)
  list(truth = truth, signal = signal, noise = noise)
}

simulation_decomposition_case <- function(base, sigma) {
  truth <- lapply(base$truth, function(component) {
    component$lambda <- 4 * component$lambda
    component
  })
  list(observed = 4 * base$signal + sigma * base$noise,
       truth = simulation_reconstruct(truth), components = truth)
}

simulation_regression_base <- function(P, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  input_dims <- c(P, 10L, 10L); output_dims <- c(10L, 10L)
  generation_dims <- c(output_dims, input_dims)
  factors <- lapply(generation_dims, function(p) {
    gram <- stats::toeplitz(.5^(0:2))
    columns <- rbind(diag(1, 3L), matrix(0, p - 3L, 3L)) %*% chol(gram)
    columns <- apply(columns, 2L, function(z) z / sqrt(sum(z^2)))
    rotation <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
    rotation %*% columns
  })
  internal <- lapply(seq_len(3L), function(i) {
    list(lambda = 30^(.5 - (i - 1L) / (3L - 1L)),
         U = lapply(factors, function(A) A[, i]))
  })
  value <- numeric(prod(generation_dims))
  for (component in internal)
    value <- value + component$lambda *
      Reduce(kronecker, rev(lapply(component$U, as.numeric)))
  truth <- aperm(array(value, generation_dims), c(3L, 4L, 5L, 1L, 2L))
  components <- lapply(internal, function(component) {
    component$U <- component$U[c(3L, 4L, 5L, 1L, 2L)]
    component
  })
  p_in <- prod(input_dims); q_out <- prod(output_dims)
  X_train <- matrix(stats::rnorm(1000L * p_in, sd = 1 / sqrt(1000L)), 1000L, p_in)
  X_test <- matrix(stats::rnorm(1000L * p_in, sd = 1 / sqrt(1000L)), 1000L, p_in)
  coefficient <- matrix(as.numeric(truth), p_in, q_out)
  Y_train_clean <- unname(X_train %*% coefficient)
  Y_test_clean <- unname(X_test %*% coefficient)
  Z_train <- matrix(stats::rnorm(1000L * q_out), 1000L, q_out)
  Z_test <- matrix(stats::rnorm(1000L * q_out), 1000L, q_out)
  Y_train <- Y_train_clean + 1 * Z_train / sqrt(1000L)
  Y_test <- Y_test_clean + 1 * Z_test / sqrt(1000L)
  response <- function(Y) array(as.numeric(Y), c(nrow(Y), output_dims))
  data <- list(X_train = X_train, X_test = X_test,
    Y_train = response(Y_train), Y_test = response(Y_test),
    Y_train_matrix = Y_train, Y_test_matrix = Y_test,
    Y_train_clean = response(Y_train_clean), Y_test_clean = response(Y_test_clean),
    Y_test_clean_matrix = Y_test_clean, truth = truth, components = components,
    input_dims = input_dims, output_dims = output_dims)
  # Reuse the realized residual arrays before multiplying the signal.
  list(data = data, unit_noise = list(train = data$Y_train - data$Y_train_clean,
                                     test = data$Y_test - data$Y_test_clean))
}

simulation_regression_case <- function(base, sigma) {
  data <- base$data
  data$truth <- 4 * base$data$truth
  data$components <- lapply(base$data$components, function(component) {
    component$lambda <- 4 * component$lambda
    component
  })
  data$Y_train_clean <- 4 * base$data$Y_train_clean
  data$Y_test_clean <- 4 * base$data$Y_test_clean
  data$Y_train <- data$Y_train_clean + sigma * base$unit_noise$train
  data$Y_test <- data$Y_test_clean + sigma * base$unit_noise$test
  data$Y_train_matrix <- matrix(data$Y_train, nrow(data$X_train))
  data$Y_test_matrix <- matrix(data$Y_test, nrow(data$X_test))
  data$Y_test_clean_matrix <- matrix(data$Y_test_clean, nrow(data$X_test))
  data
}

simulation_regression_initialization <- function(data, seed) {
  moment <- crossprod(data$X_train * 1, data$Y_train_matrix * 1) / 1
  moment <- array(as.numeric(moment), c(data$input_dims, data$output_dims))
  initial <- SegreCP::cp_initialize(moment, 3L, method = "cpca",
                                  split = c(4L, 5L), seed = seed)
  structure(list(components = initial$comps,
    initialization = initial[names(initial) != "comps"], moment = moment,
    metadata = list(method = "CPCA", split = c(4L, 5L), design_variance = 1,
      normalization = "none", scale = 1, mode_order = "input_then_output")),
    class = "regression_cpca_initialization")
}

simulation_metrics <- function(components, data, task) {
  coefficient <- simulation_reconstruct(components)
  result <- data.frame(relative_recovery_error =
    sqrt(sum((coefficient - data$truth)^2) / sum(data$truth^2)))
  if (task == "decomposition") {
    result$objective <- .5 * sum((coefficient - data$observed)^2)
  } else {
    input_modes <- seq_along(data$input_dims)
    product <- function(U) if (length(U)) Reduce(kronecker, rev(U)) else 1
    input <- do.call(cbind, lapply(components, function(z) product(z$U[input_modes])))
    output <- do.call(cbind, lapply(components, function(z) product(z$U[-input_modes])))
    weights <- vapply(components, `[[`, numeric(1L), "lambda")
    predict <- function(X) tcrossprod(sweep(X %*% input, 2L, weights, "*"), output)
    train <- predict(data$X_train) - matrix(data$Y_train, nrow(data$X_train))
    test_prediction <- predict(data$X_test)
    test <- test_prediction - matrix(data$Y_test, nrow(data$X_test))
    clean <- test_prediction - matrix(data$Y_test_clean, nrow(data$X_test))
    result$objective <- .5 * sum(train^2)
    result$train_RMSE <- sqrt(mean(train^2))
    result$test_RMSE <- sqrt(mean(test^2))
    result$clean_test_RMSE <- sqrt(mean(clean^2))
    result$relative_clean_prediction_error <- sqrt(sum(clean^2) / sum(data$Y_test_clean^2))
  }
  result
}

run_signal4_simulation <- function(task, script) {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste("Usage: Rscript", basename(script),
    "--output-dir=/new/path [--full] [--p-index=1] [--replication=1] [--iterations=N]",
    "\nDefault: first P, one paired replication, noise 0 and 1; 3 decomposition or 2 regression iterations.",
    "\n--full: all P, 20 replications, 10 iterations. Selectors apply only to smoke runs.")
  if ("--help" %in% args) { cat(usage, "\n"); return(invisible(NULL)) }
  value <- function(name, default = NULL) {
    matches <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(matches) > 1L) stop("Repeated option: ", name)
    if (length(matches)) sub(paste0("^--", name, "="), "", matches) else default
  }
  full <- "--full" %in% args
  known <- grepl("^--(output-dir|p-index|replication|iterations)=", args) | args == "--full"
  if (!length(args) || !all(known) ||
      (full && any(grepl("^--(p-index|replication|iterations)=", args)))) stop(usage)
  settings <- simulation_settings(task, full, as.integer(value("p-index", "1")),
    as.integer(value("replication", "1")),
    if (is.null(value("iterations"))) NULL else as.integer(value("iterations")))
  output <- value("output-dir")
  if (is.null(output) || !startsWith(output, "/") || file.exists(output))
    stop("Choose a new absolute --output-dir outside the package checkout.")
  package_root <- normalizePath(file.path(dirname(script), "..", ".."))
  parent <- normalizePath(dirname(output), mustWork = TRUE)
  output <- file.path(parent, basename(output))
  if (startsWith(output, paste0(package_root, "/"))) stop("Keep outputs outside the package checkout.")
  dir.create(output)
  for (folder in c("cases", "fits")) dir.create(file.path(output, folder))
  saveRDS(settings, file.path(output, "settings.rds"), version = 3)
  utils::write.csv(data.frame(setting = names(settings),
    value = vapply(settings, paste, character(1L), collapse = ";")),
    file.path(output, "settings.csv"), row.names = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1")
  paths <- fits <- inputs <- list()
  for (p_index in settings$p_indices) for (replication in settings$replications) {
    P <- settings$P[[p_index]]
    seed <- settings$master_seed + settings$task_offset + 1000L * p_index + replication
    base <- if (task == "decomposition") simulation_decomposition_base(P, seed) else
      simulation_regression_base(P, seed)
    for (sigma in settings$noise_sd) {
      key <- sprintf("P%03d_rep%02d_c4_sigma%d", P, replication, sigma)
      data <- if (task == "decomposition") simulation_decomposition_case(base, sigma) else
        simulation_regression_case(base, sigma)
      initialization <- if (task == "decomposition")
        SegreCP::cp_initialize(data$observed, 3L, method = "cpca", split = c(1L, 2L)) else
        simulation_regression_initialization(data, seed)
      initial <- if (task == "decomposition") initialization$comps else initialization$components
      fields <- if (task == "decomposition") c("observed", "truth", "components") else
        c("X_train", "X_test", "Y_train", "Y_test", "Y_train_clean", "Y_test_clean", "truth", "components")
      hashes <- vapply(data[fields], simulation_hash, character(1L))
      input <- data.frame(key = key, P = P, replication = replication, seed = seed,
        noise_sd = sigma, object = names(hashes), sha256 = unname(hashes),
        initialization_sha256 = simulation_hash(initialization),
        initial_components_sha256 = simulation_hash(initial))
      inputs[[key]] <- input
      saveRDS(list(seed = seed, noise_sd = sigma, P = P, input_hashes = hashes,
        truth_components = data$components, initialization = initialization),
        file.path(output, "cases", paste0(key, ".rds")), version = 3)
      for (method in settings$methods) {
        control <- SegreCP::segre_control(profile = "manuscript",
          max_iter = settings$iterations, tolerance = 0,
          step_size = if (task == "decomposition" && method == "RGN") 1 else .5,
          step_rule = if (task == "regression" && method == "RGD") "exact_quadratic" else "fixed",
          component_order = "decreasing_weight", ridge = 0, normalize_regression = "none")
        fit <- if (task == "decomposition")
          SegreCP::segre_decompose(data$observed, 3L, tolower(method), init = initial,
            control = control, keep_history = TRUE) else
          SegreCP::segre_regress(data$X_train, data$Y_train, 3L, data$input_dims,
            data$output_dims, method = tolower(method), init = initial, seed = seed,
            design_variance = 1, control = control, keep_history = TRUE)
        fit_key <- paste(key, method, sep = "__")
        path <- do.call(rbind, lapply(fit$component_history, simulation_metrics, data = data, task = task))
        paths[[fit_key]] <- data.frame(task = task, P = P, replication = replication,
          seed = seed, signal_multiplier = 4, noise_sd = sigma, method = method,
          iteration = seq_len(nrow(path)) - 1L, path)
        fits[[fit_key]] <- data.frame(key = fit_key, method = method,
          actual_iterations = fit$iterations, ridge = fit$control$ridge,
          stopping_reason = fit$stopping_reason, engine = settings$engine)
        saveRDS(list(component_history = fit$component_history, control = fit$control,
          accepted_steps = fit$accepted_steps, objective = fit$objective,
          iterations = fit$iterations, stopping_reason = fit$stopping_reason),
          file.path(output, "fits", paste0(fit_key, ".rds")), version = 3)
      }
      message("Completed ", task, " ", key)
    }
    utils::write.csv(do.call(rbind, paths), file.path(output, "trajectories.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, fits), file.path(output, "fits.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, inputs), file.path(output, "inputs.csv"), row.names = FALSE)
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(output, "sessionInfo.txt"))
  message("Saved fresh authors-method results to ", output)
  invisible(output)
}
