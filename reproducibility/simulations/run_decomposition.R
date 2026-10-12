#!/usr/bin/env Rscript

# Fresh RGD/RGN rerun of the rebuttal-stage CP decomposition design.
# This script generates tensors anew; it does not read stored trajectories.

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript run_decomposition.R --output-dir=PATH [--full]",
  "",
  "Without --full: one replication and three sweeps (a quick smoke run).",
  "With --full: the rebuttal design of 20 replications and 10 sweeps.",
  "The output directory must not already exist. SegreCP must be installed.",
  sep = "\n"
)
if ("--help" %in% args) {
  cat(usage, "\n")
  quit(status = 0L)
}
value_args <- grep("^--output-dir=", args, value = TRUE)
if (length(value_args) != 1L ||
    length(args) != length(value_args) + as.integer("--full" %in% args) ||
    anyDuplicated(args)) {
  stop(usage, call. = FALSE)
}
output_dir <- sub("^--output-dir=", "", value_args)
if (!nzchar(output_dir) || file.exists(output_dir)) {
  stop("Provide a new, non-existing --output-dir path.", call. = FALSE)
}
if (!requireNamespace("SegreCP", quietly = TRUE)) {
  stop("Install SegreCP before running this script.", call. = FALSE)
}

full_run <- "--full" %in% args
replications <- if (full_run) 20L else 1L
sweeps <- if (full_run) 10L else 3L
dims <- c(60L, 30L, 30L)
rank <- 3L
kappa <- 30
coherence <- 0.5
noise_sd <- 1
signal_multiplier <- 3
master_seed <- 20260829L
cpca_split <- 1L
rgd_step <- 0.5
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# A random orthonormal embedding with the prescribed Toeplitz component Gram.
generate_coherent_vectors <- function(p, r, rho) {
  gaussian <- matrix(stats::rnorm(p * r), nrow = p, ncol = r)
  qr.Q(qr(gaussian)) %*% chol(toeplitz(rho^(0:(r - 1L))))
}

generate_truth <- function() {
  factors <- lapply(dims, generate_coherent_vectors,
                    r = rank, rho = coherence)
  lapply(seq_len(rank), function(i) {
    exponent <- 0.5 - (i - 1L) / (rank - 1L)
    list(
      lambda = signal_multiplier * sqrt(rank) * kappa^exponent *
        max(dims)^(length(dims) / 4),
      U = lapply(factors, function(A) A[, i])
    )
  })
}

# Order the HOSVD start by signed core weights for the first sweep.
hosvd_first_order <- function(observed) {
  factors <- lapply(seq_along(dims), function(mode) {
    permutation <- c(mode, setdiff(seq_along(dims), mode))
    unfolding <- matrix(aperm(observed, permutation), nrow = dims[[mode]])
    svd(unfolding, nu = rank, nv = 0L)$u[, seq_len(rank), drop = FALSE]
  })
  signed_weights <- vapply(seq_len(rank), function(i) {
    direction <- Reduce(
      kronecker,
      rev(lapply(factors, function(A) A[, i]))
    )
    sum(as.numeric(observed) * as.numeric(direction))
  }, numeric(1L))
  order(signed_weights, decreasing = TRUE)
}

fit_path <- function(observed, truth, initial, method, initialization) {
  fit_sweeps <- function(start, n, ordering) {
    control <- SegreCP::segre_control(
      profile = "manuscript", max_iter = n, tolerance = 0,
      step_size = if (method == "RGD") rgd_step else 1,
      step_rule = "fixed", component_order = ordering
    )
    SegreCP::segre_decompose(
      observed, rank, method = tolower(method), init = start,
      control = control, keep_history = TRUE
    )
  }
  if (initialization == "CPCA") {
    fit <- fit_sweeps(initial, sweeps, "decreasing_weight")
    return(SegreCP::convergence_data(fit, truth))
  }
  first_order <- hosvd_first_order(observed)
  first <- fit_sweeps(initial[first_order], 1L, "fixed")
  path <- SegreCP::convergence_data(first, truth)
  if (sweeps > 1L) {
    later <- fit_sweeps(first$components, sweeps - 1L,
                        "decreasing_weight")
    later_path <- SegreCP::convergence_data(later, truth)
    path <- rbind(path, later_path[-1L, , drop = FALSE])
    path$iteration <- 0:sweeps
  }
  path
}

trajectory <- vector("list", replications * 4L)
terminal <- vector("list", replications * 4L)
row_id <- 0L
for (replication in seq_len(replications)) {
  seed <- master_seed + replication
  set.seed(seed)
  truth <- generate_truth()
  truth_vector <- as.numeric(SegreCP::reconstruct_cp(truth))
  observed_vector <- truth_vector +
    stats::rnorm(prod(dims), sd = noise_sd)
  observed <- array(observed_vector, dim = dims)
  initializers <- list(
    CPCA = SegreCP::cp_initialize(
      observed, rank, method = "cpca", split = cpca_split
    )$comps,
    HOSVD = SegreCP::cp_initialize(
      observed, rank, method = "hosvd"
    )$comps
  )
  for (initialization in names(initializers)) {
    for (method in c("RGD", "RGN")) {
      row_id <- row_id + 1L
      path <- fit_path(observed, truth, initializers[[initialization]],
                       method, initialization)
      trajectory[[row_id]] <- data.frame(
        replication = replication, seed = seed,
        initialization = initialization, method = method,
        iteration = path$iteration,
        relative_frobenius_error = path$relative_frobenius_error,
        objective = path$objective
      )
      terminal[[row_id]] <- data.frame(
        replication = replication, seed = seed,
        initialization = initialization, method = method,
        final_relative_frobenius_error =
          tail(path$relative_frobenius_error, 1L),
        final_objective = tail(path$objective, 1L)
      )
    }
  }
  message(sprintf("Completed replication %d/%d", replication, replications))
}
trajectory <- do.call(rbind, trajectory)
terminal <- do.call(rbind, terminal)
expected_rows <- replications * 2L * 2L * (sweeps + 1L)
checks <- data.frame(
  check = c("expected trajectory rows", "finite trajectories",
            "paired RGD/RGN fits", "nonincreasing observed loss",
            "unique replication seeds"),
  passed = c(
    nrow(trajectory) == expected_rows,
    all(is.finite(trajectory$relative_frobenius_error)) &&
      all(is.finite(trajectory$objective)),
    all(vapply(split(terminal$method,
                     list(terminal$replication, terminal$initialization),
                     drop = TRUE),
               function(z) identical(sort(z), c("RGD", "RGN")), logical(1L))),
    all(vapply(split(trajectory, list(trajectory$replication,
                                     trajectory$initialization,
                                     trajectory$method), drop = TRUE),
               function(z) all(diff(z$objective[order(z$iteration)]) <=
                                1e-8 * max(1, abs(z$objective[[1L]]))),
               logical(1L))),
    !anyDuplicated(unique(trajectory[c("replication", "seed")])$seed)
  )
)

# Compare fresh trajectories with the rebuttal-stage records.
compare_reference <- function(directory) {
  reference_file <- file.path(directory, "trajectory_raw.csv")
  if (!file.exists(reference_file)) return(NULL)
  reference <- utils::read.csv(reference_file, stringsAsFactors = FALSE)
  reference <- reference[
    reference$method %in% c("RGD", "RGN") &
      reference$replication <= replications &
      reference$iteration <= sweeps, , drop = FALSE
  ]
  keys <- c("replication", "seed", "initialization", "method", "iteration")
  joined <- merge(
    trajectory[c(keys, "relative_frobenius_error")],
    reference[c(keys, "relative_frobenius_error")],
    by = keys, suffixes = c("_fresh", "_reference")
  )
  if (nrow(joined) != nrow(trajectory)) {
    stop("The rebuttal trajectory has incomplete keys.", call. = FALSE)
  }
  joined$error_difference <- abs(
    joined$relative_frobenius_error_fresh -
      joined$relative_frobenius_error_reference
  )
  means <- aggregate(
    relative_frobenius_error ~ initialization + method + iteration,
    trajectory, mean
  )
  reference_means <- aggregate(
    relative_frobenius_error ~ initialization + method + iteration,
    reference, mean
  )
  means <- merge(means, reference_means,
                 by = c("initialization", "method", "iteration"),
                 suffixes = c("_fresh", "_reference"))
  means$absolute_difference <- abs(
    means$relative_frobenius_error_fresh -
      means$relative_frobenius_error_reference
  )
  checkpoints <- means$iteration %in% c(3L, 5L, 10L)
  checkpoint_difference <- if (any(checkpoints)) {
    max(means$absolute_difference[checkpoints])
  } else NA_real_
  checkpoint_rank_preserved <- NA
  if (any(checkpoints)) {
    groups <- split(means[checkpoints, , drop = FALSE],
                    list(means$initialization[checkpoints],
                         means$iteration[checkpoints]), drop = TRUE)
    checkpoint_rank_preserved <- all(vapply(groups, function(z) {
      identical(
        z$method[order(z$relative_frobenius_error_fresh, z$method)],
        z$method[order(z$relative_frobenius_error_reference, z$method)]
      )
    }, logical(1L)))
  }
  data.frame(
    scope = if (full_run) "full_design" else "smoke_subset",
    matched_rows = nrow(joined),
    cpca_max_absolute_error_difference = max(joined$error_difference[
      joined$initialization == "CPCA"
    ]),
    hosvd_initial_max_absolute_error_difference = max(
      joined$error_difference[joined$initialization == "HOSVD" &
                                joined$iteration == 0L]
    ),
    mean_path_max_absolute_error_difference = max(means$absolute_difference),
    checkpoint_max_absolute_error_difference = checkpoint_difference,
    checkpoint_rgd_rgn_rank_preserved = checkpoint_rank_preserved
  )
}
reference_comparison <- compare_reference(system.file(
  "extdata", "rebuttal_init", "decomp_ref", package = "SegreCP"
))

dir.create(output_dir, recursive = TRUE)
utils::write.csv(trajectory, file.path(output_dir, "trajectory.csv"),
                 row.names = FALSE)
utils::write.csv(terminal, file.path(output_dir, "terminal_accuracy.csv"),
                 row.names = FALSE)
utils::write.csv(checks, file.path(output_dir, "validation_checks.csv"),
                 row.names = FALSE)
if (!is.null(reference_comparison)) {
  utils::write.csv(reference_comparison,
                   file.path(output_dir, "reference_comparison.csv"),
                   row.names = FALSE)
}
settings <- data.frame(
  setting = c("mode", "dimensions", "rank", "condition_number",
              "coherence", "noise_sd", "signal_multiplier", "replications",
              "sweeps", "master_seed", "rng_kind", "cpca_split", "rgd_step",
              "package_version"),
  value = c(if (full_run) "full" else "smoke", paste(dims, collapse = "x"),
            rank, kappa, coherence, noise_sd, signal_multiplier,
            replications, sweeps, master_seed, paste(RNGkind(), collapse = ";"),
            cpca_split, rgd_step, as.character(utils::packageVersion("SegreCP")))
)
utils::write.csv(settings, file.path(output_dir, "settings.csv"),
                 row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "session_info.txt"))
message("Wrote fresh author-method results to: ", normalizePath(output_dir))
if (!all(checks$passed)) {
  stop("At least one numerical validation check failed; inspect the output.",
       call. = FALSE)
}
