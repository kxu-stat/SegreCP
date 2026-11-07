# Simulations of Section 5 and Appendix C.
#
# Run from the repository root:
#   Rscript reproducibility/simulations/run_simulations.R <study> <out_dir> [replications]
# study:
#   regression           tensor-on-tensor regression, 20 iterations
#                        (Figures 1, 2 and 5, Table 6)
#   regression_long      CP-ALS, RGD and RGN for 100 iterations on the same
#                        data sets (Table 4, Figure 5)
#   decomposition        CP decomposition, 20 iterations (Figures 3 and 4)
#   decomposition_long   CP decomposition, 100 iterations, on the same data
#                        sets (coherence-0.75 column of Figures 3 and 4,
#                        Tables 5 and 7, Figure 5)
#   initialization       CPCA, T-HOSVD and random starts in the settings with
#                        noise (Tables 8 and 9)
#   timing_regression, timing_decomposition
#                        running time per iteration (Tables 3-5, Figure 5);
#                        run with replications 1:5 on one core (SEGRECP_CORES=1,
#                        OPENBLAS_NUM_THREADS=1) with nothing else running.
# replications defaults to 1:20. Cases run in parallel on SEGRECP_CORES cores
# (default 1).
#
# Each (coherence, snr, replication) case is saved as one .rds file in
# out_dir/cases_<study>/ (existing files are skipped, so a run can be resumed).
# At the end all cases are collected into
#   out_dir/raw_<study>.rds       every case, with all retained raw data;
#   out_dir/results_<study>.csv   one row per case, method and iteration.
#
# Design: rank 3, coherence {0, .25, .5, .75}, signal-to-noise ratio of the
# weakest component s in {Inf, 30, 10, 3}, shared CPCA start for CP-ALS, RGD
# and RGN, ridge 0. The seeds are set in decomposition_data() and
# regression_data() below.
#
# Raw data kept for every case: settings and seed, true and initial CP
# components, per-iteration error and objective of every method, the CP
# component history of RGD / RGN (and of CP-ALS for decomposition), the final
# components of the other CP methods, and wall-clock seconds of every fit. Full
# tensors are not kept; they can be rebuilt from the seed.

args <- commandArgs(trailingOnly = TRUE)
study <- args[[1]]
out_dir <- args[[2]]
replications <- if (length(args) >= 3) eval(parse(text = args[[3]])) else 1:20
case_dir <- file.path(out_dir, paste0("cases_", study))
dir.create(case_dir, recursive = TRUE, showWarnings = FALSE)

library(SegreCP)
code <- file.path("reproducibility", "code")
source(file.path(code, "regression_adapters.R"))
source(file.path(code, "simulations.R"))
source(file.path(code, "decomposition_comparators.R"))
source(file.path(code, "tucker_rgn.R"))

ITERATIONS <- if (study %in% c("regression_long", "decomposition_long")) 100 else 20
COHERENCE <- c(0, 0.25, 0.5, 0.75)          # also used in the seeds
# regression_long runs only the three methods that share the CPCA start
REGRESSION_METHODS <- if (study == "regression_long") c("CP-ALS", "RGD", "RGN") else
  c("CP-ALS", "RGD", "RGN", "RRR", "Tucker-RGD", "TT-IHT")
SNR <- if (study == "initialization") c(30, 10, 3) else c(Inf, 30, 10, 3)
RANK <- 3
KAPPA <- 30

relative_error <- function(estimate, truth) {
  sqrt(sum((estimate - truth)^2) / sum(truth^2))
}

# Evaluate expr and return its value with the elapsed seconds.
timed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- expr
  list(value = value, seconds = proc.time()[["elapsed"]] - start)
}

path_rows <- function(method, error, objective = NA) {
  data.frame(method = method, iteration = seq_along(error) - 1,
             error = error, objective = objective)
}

khatri_rao <- function(matrices) {      # first matrix varies slowest
  out <- matrices[[1]]
  for (B in matrices[-1]) {
    out <- do.call(cbind, lapply(seq_len(ncol(out)), function(j) {
      kronecker(out[, j], B[, j])
    }))
  }
  out
}

# CP-ALS: the update of rTensor::cp (same unfolding, Khatri-Rao product and
# normal equations), started from the shared CPCA start. The convergence check
# of rTensor::cp (a full reconstruction after every mode) is left out.
cp_als_decomposition <- function(Y, init, iterations) {
  tnsr <- rTensor::as.tensor(Y)
  d <- tnsr@num_modes
  factors <- lapply(seq_len(d), function(m) sapply(init, function(z) z$U[[m]]))
  weights <- sapply(init, `[[`, "lambda")
  unfolded <- lapply(seq_len(d), function(m) rTensor::rs_unfold(tnsr, m = m)@data)
  as_components <- function() lapply(seq_along(weights), function(i) {
    list(lambda = weights[i], U = lapply(factors, function(A) A[, i]))
  })
  history <- list(as_components())
  for (t in seq_len(iterations)) {
    for (m in seq_len(d)) {
      gram <- rTensor::hadamard_list(lapply(factors[-m], crossprod))
      updated <- unfolded[[m]] %*%
        rTensor::khatri_rao_list(factors[-m], reverse = TRUE) %*% solve(gram)
      norms <- sqrt(colSums(updated^2))
      factors[[m]] <- sweep(updated, 2, norms, "/")
      weights <- norms
    }
    history[[t + 1]] <- as_components()
  }
  history
}

# ------------------------------------------------------- decomposition (d = 4)
decomposition_data <- function(coherence, snr, replication) {
  dims <- c(60, 30, 30, 30)
  seed <- 20261104 + 1000 * 4 + 100 * match(coherence, COHERENCE) + replication
  set.seed(seed)
  factors <- .coherent_factors(dims, RANK, coherence)
  lambda <- sqrt(RANK) * dims[1] *
    KAPPA^(0.5 - (seq_len(RANK) - 1) / (RANK - 1))
  truth_components <- lapply(seq_len(RANK), function(i) {
    list(lambda = lambda[i], U = lapply(factors, function(A) A[, i]))
  })
  truth <- SegreCP::reconstruct_cp(truth_components)
  noise <- array(rnorm(prod(dims)), dims)
  df <- 1 + sum(dims - 1)
  sigma <- if (is.infinite(snr)) 0 else min(lambda) / (snr * sqrt(df))
  list(dims = dims, seed = seed, sigma = sigma, truth = truth,
       truth_components = truth_components, Y = truth + sigma * noise)
}

decomposition_case <- function(coherence, snr, replication) {
  data <- decomposition_data(coherence, snr, replication)
  Y <- data$Y
  init <- SegreCP::cp_initialize(Y, RANK, method = "cpca")$comps
  path_of <- function(method, history) {
    est <- lapply(history, SegreCP::reconstruct_cp)
    path_rows(method, sapply(est, relative_error, truth = data$truth),
              sapply(est, function(A) 0.5 * sum((Y - A)^2)))
  }
  fits <- list(
    "CP-ALS" = timed(cp_als_decomposition(Y, init, ITERATIONS)),
    "RGD" = timed(SegreCP::segre_decompose(Y, RANK, "rgd", init = init,
      control = segre_control(max_iter = ITERATIONS),
      keep_history = TRUE)$component_history),
    "RGN" = timed(SegreCP::segre_decompose(Y, RANK, "rgn", init = init,
      control = segre_control(max_iter = ITERATIONS),
      keep_history = TRUE)$component_history))
  paths <- lapply(names(fits), function(m) path_of(m, fits[[m]]$value))
  components <- lapply(fits, `[[`, "value")
  seconds <- sapply(fits, `[[`, "seconds")
  for (method in c("Orth-ALS", "ICO", "NLS-GN", "GN-CG", "CP-ALS+ELS",
                   "Hybrid AMDM")) {
    fit <- timed(tryCatch(fit_decomposition_comparator(method, Y, init,
      ITERATIONS, truth = data$truth_components, diagnostics = TRUE),
      error = function(e) { message(method, ": ", conditionMessage(e)); NULL }))
    if (is.null(fit$value)) next
    paths[[length(paths) + 1]] <- path_rows(method,
      fit$value$objective$relative_frobenius_error, fit$value$objective$objective)
    components[[method]] <- fit$value$components       # final components only
    seconds[method] <- fit$seconds
  }
  list(settings = data.frame(task = "decomposition_d4", coherence = coherence,
         snr = snr, replication = replication, seed = data$seed,
         sigma = data$sigma, dims = paste(data$dims, collapse = "x"),
         rank = RANK, kappa = KAPPA, iterations = ITERATIONS),
       paths = do.call(rbind, paths),
       truth_components = data$truth_components, init_components = init,
       components = components, seconds = seconds)
}

# ---------------------------------------------------- tensor-on-tensor regression
regression_data <- function(coherence, snr, replication) {
  input_dims <- c(10, 10, 10)
  output_dims <- c(10, 10)
  n <- 1000
  seed <- 20261204 + 100 * match(coherence, COHERENCE) + replication
  df <- 1 + sum(c(input_dims, output_dims) - 1)
  sigma <- if (is.infinite(snr)) 0 else KAPPA^(-0.5) / (snr * sqrt(df / n))
  data <- .tensor_regression_data(input_dims, output_dims, RANK, KAPPA,
    coherence, sigma, n, n, seed)
  c(data, list(input_dims = input_dims, output_dims = output_dims, n = n,
               seed = seed, sigma = sigma))
}

regression_case <- function(coherence, snr, replication) {
  data <- regression_data(coherence, snr, replication)
  start <- regression_start(data$X_train, data$Y_train, data$input_dims,
                            data$output_dims, RANK, split = c(4, 5))
  # RGD step 1 = sigma_X^(-2): the design entries have variance 1/n.
  fit <- function(method) {
    fit_regression_method(method, data$X_train, data$Y_train, data$input_dims,
      data$output_dims, rank = RANK, start = start, iterations = ITERATIONS,
      seed = data$seed)
  }
  paths <- list(); components <- list(); seconds <- numeric(0)
  for (method in REGRESSION_METHODS) {
    f <- timed(tryCatch(fit(method), error = function(e) {
      message(method, ": ", conditionMessage(e)); NULL }))
    if (is.null(f$value)) next
    history <- f$value$coefficient_history
    if (is.null(history)) {
      history <- lapply(f$value$component_history, SegreCP::reconstruct_cp)
    }
    objective <- f$value$objective
    if (length(objective) != length(history)) objective <- NA     # RRR
    paths[[length(paths) + 1]] <- path_rows(method,
      sapply(history, relative_error, truth = data$truth), objective)
    components[[method]] <- if (method == "CP-ALS") f$value$components else
      f$value$component_history                        # NULL for Tucker / TT
    seconds[method] <- f$seconds
  }
  if (study == "regression") {
    tr <- timed(tucker_rgn_tot(data$X_train, matrix(data$Y_train, data$n),
                               data$input_dims, data$output_dims, rep(RANK, 5),
                               ITERATIONS))
    paths[[length(paths) + 1]] <- path_rows("Tucker-RGN",
      sapply(tr$value$coefficient_history, relative_error, truth = data$truth))
    seconds["Tucker-RGN"] <- tr$seconds
    hopls <- timed(tryCatch(fit("HOPLS"), error = function(e) NULL))
    if (!is.null(hopls$value)) {
      paths[[length(paths) + 1]] <- data.frame(method = "HOPLS", iteration = NA,
        error = relative_error(hopls$value$coefficient, data$truth),
        objective = NA)
      seconds["HOPLS"] <- hopls$seconds
    }
  }
  list(settings = data.frame(task = "regression", coherence = coherence,
         snr = snr, replication = replication, seed = data$seed,
         sigma = data$sigma, n = data$n,
         dims = paste(c(data$input_dims, data$output_dims), collapse = "x"),
         rank = RANK, kappa = KAPPA, iterations = ITERATIONS),
       paths = do.call(rbind, paths),
       truth_components = data$components, init_components = start,
       components = components, seconds = seconds)
}

# ------------------------------------------------------------------------ timing
# Studies timing_decomposition and timing_regression: the fits above with the
# error diagnostics switched off. Run them one case at a time on one thread
# (SEGRECP_CORES=1, OPENBLAS_NUM_THREADS=1). Each row of `paths` is the elapsed
# seconds of one fit of ITERATIONS iterations; the shared CPCA start is timed
# as its own row ("CPCA start") and is not included in the other rows.
decomposition_timing_case <- function(coherence, snr, replication) {
  data <- decomposition_data(coherence, snr, replication)
  Y <- data$Y
  cpca <- timed(SegreCP::cp_initialize(Y, RANK, method = "cpca")$comps)
  init <- cpca$value
  seconds <- c("CPCA start" = cpca$seconds)
  seconds["CP-ALS"] <- timed(cp_als_decomposition(Y, init, ITERATIONS))$seconds
  for (method in c("RGD", "RGN")) {
    seconds[method] <- timed(SegreCP::segre_decompose(Y, RANK, tolower(method),
      init = init, control = segre_control(max_iter = ITERATIONS)))$seconds
  }
  for (method in c("Orth-ALS", "ICO", "NLS-GN", "GN-CG", "CP-ALS+ELS",
                   "Hybrid AMDM")) {
    seconds[method] <- timed(fit_decomposition_comparator(method, Y, init,
      ITERATIONS, truth = NULL, diagnostics = FALSE))$seconds
  }
  list(settings = data.frame(task = "decomposition_d4", coherence = coherence,
         snr = snr, replication = replication, seed = data$seed,
         iterations = ITERATIONS),
       paths = data.frame(method = names(seconds), seconds = seconds))
}

regression_timing_case <- function(coherence, snr, replication) {
  data <- regression_data(coherence, snr, replication)
  cpca <- timed(regression_start(data$X_train, data$Y_train, data$input_dims,
                                 data$output_dims, RANK, split = c(4, 5)))
  seconds <- c("CPCA start" = cpca$seconds)
  for (method in c("CP-ALS", "RGD", "RGN", "RRR", "Tucker-RGD", "TT-IHT",
                   "HOPLS")) {
    prepared <- prepare_regression_fit(method, data$X_train, data$Y_train,
      data$input_dims, data$output_dims, rank = RANK, start = cpca$value,
      iterations = ITERATIONS, seed = data$seed, record_path = FALSE)
    seconds[method] <- timed(run_regression_fit(prepared))$seconds
  }
  seconds["Tucker-RGN"] <- timed(tucker_rgn_tot(data$X_train,
    matrix(data$Y_train, data$n), data$input_dims, data$output_dims,
    rep(RANK, 5), ITERATIONS))$seconds
  list(settings = data.frame(task = "regression", coherence = coherence,
         snr = snr, replication = replication, seed = data$seed,
         iterations = ITERATIONS),
       paths = data.frame(method = names(seconds), seconds = seconds))
}

# ------------------------------------------------ initialization sensitivity
# Study initialization: CP-ALS, RGD and RGN started from CPCA, from T-HOSVD and
# from random unit factors (weights fitted to the data, seed = data seed), on
# the decomposition and regression data sets with noise (same seeds as the
# studies above), 20 iterations. Error = relative Frobenius error.
STARTS <- c("CPCA", "T-HOSVD", "Random")
initialization_case <- function(coherence, snr, replication) {
  paths <- list()
  add_path <- function(problem, start, method, history) {
    error <- sapply(history, function(B) relative_error(B, truth))
    paths[[length(paths) + 1]] <<- data.frame(problem = problem,
      start = start, path_rows(method, error))
  }
  # decomposition
  data <- decomposition_data(coherence, snr, replication)
  truth <- data$truth
  for (start in STARTS) {
    init <- SegreCP::cp_initialize(data$Y, RANK, seed = data$seed,
      method = switch(start, "CPCA" = "cpca", "T-HOSVD" = "hosvd",
                      "Random" = "random"))$comps
    add_path("decomposition", start, "CP-ALS", lapply(
      cp_als_decomposition(data$Y, init, ITERATIONS), SegreCP::reconstruct_cp))
    for (method in c("RGD", "RGN")) {
      fit <- SegreCP::segre_decompose(data$Y, RANK, tolower(method), init = init,
        control = segre_control(max_iter = ITERATIONS), keep_history = TRUE)
      add_path("decomposition", start, method,
               lapply(fit$component_history, SegreCP::reconstruct_cp))
    }
  }
  # regression
  data <- regression_data(coherence, snr, replication)
  truth <- data$truth
  for (start in STARTS) {
    init <- regression_start(data$X_train, data$Y_train, data$input_dims,
      data$output_dims, RANK, split = c(4, 5), seed = data$seed,
      method = switch(start, "CPCA" = "cpca", "T-HOSVD" = "hosvd",
                      "Random" = "random"))
    for (method in c("CP-ALS", "RGD", "RGN")) {
      fit <- fit_regression_method(method, data$X_train, data$Y_train,
        data$input_dims, data$output_dims, rank = RANK, start = init,
        iterations = ITERATIONS)
      history <- fit$coefficient_history
      if (is.null(history)) {
        history <- lapply(fit$component_history, SegreCP::reconstruct_cp)
      }
      add_path("regression", start, method, history)
    }
  }
  list(settings = data.frame(task = "initialization", coherence = coherence,
         snr = snr, replication = replication, iterations = ITERATIONS),
       paths = do.call(rbind, paths))
}

# --------------------------------------------------------------------------- run
cases <- expand.grid(replication = replications, snr = SNR,
                     coherence = COHERENCE)
run_case <- switch(study, decomposition = decomposition_case,
                   decomposition_long = decomposition_case,
                   regression = regression_case,
                   regression_long = regression_case,
                   timing_decomposition = decomposition_timing_case,
                   timing_regression = regression_timing_case,
                   initialization = initialization_case)
case_file <- function(cs) file.path(case_dir, sprintf("coh%.2f_snr%s_rep%02d.rds",
                                    cs$coherence, cs$snr, cs$replication))
invisible(parallel::mclapply(seq_len(nrow(cases)), function(k) {
  file <- case_file(cases[k, ])
  if (file.exists(file)) return(NULL)
  res <- tryCatch(run_case(cases$coherence[k], cases$snr[k],
                           cases$replication[k]),
    error = function(e) { message("case ", k, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) saveRDS(res, file)
  NULL
}, mc.cores = as.integer(Sys.getenv("SEGRECP_CORES", "1")), mc.preschedule = FALSE))

# Collect every finished case of this study (all replications run so far).
files <- sort(list.files(case_dir, pattern = "[.]rds$", full.names = TRUE))
raw <- lapply(files, readRDS)
saveRDS(raw, file.path(out_dir, paste0("raw_", study, ".rds")))
results <- do.call(rbind, lapply(raw, function(r) {
  cbind(r$settings[c("task", "coherence", "snr", "replication")], r$paths)
}))
write.csv(results, file.path(out_dir, paste0("results_", study, ".csv")),
          row.names = FALSE)
message(length(files), " cases collected for ", study)
