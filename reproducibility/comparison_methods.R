#!/usr/bin/env Rscript

# Standalone comparison-method entry point for manuscript experiments.
#
# IMPORTANT SCOPE BOUNDARY
# ------------------------
# This file lives in the repository-only reproducibility/ directory. That
# directory is excluded from package builds, so this script is not installed,
# exported, or supported as part of the SegreCP R-package API. It contains
# experiment-only comparison code. The authors' RGD, RGN, CPCA, and HOSVD
# routines belong in the package and in inst/scripts/ours_methods.R; they must
# not be copied into this file.
#
# The archived rebuttal studies used several task-specific comparators.  Rather
# than silently presenting them as package estimators, the catalog below records
# the exact repository worker for each experiment.  The CP-ALS path that was
# formerly mixed into the package's initialization runner is implemented here so
# that the public package contains no competing optimizer.

.comparison_script_path <- local({
  # `source()` stores the sourced file in an `ofile` field on one of the active
  # evaluation frames. Searching all frames is more reliable than assuming a
  # fixed frame number and makes this entry point independent of `getwd()`.
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)
  }), use.names = FALSE)
  source_files <- source_files[!is.na(source_files) & nzchar(source_files)]
  if (!length(source_files)) {
    NA_character_
  } else {
    tryCatch(
      normalizePath(utils::tail(source_files, 1L), mustWork = TRUE),
      error = function(e) NA_character_
    )
  }
})

.comparison_ancestors <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  ancestors <- character()
  repeat {
    ancestors <- c(ancestors, path)
    parent <- dirname(path)
    if (identical(parent, path)) break
    path <- parent
  }
  ancestors
}

comparison_repository_root <- function() {
  starting_directories <- c(
    if (!is.na(.comparison_script_path)) dirname(.comparison_script_path)
    else character(),
    getwd()
  )
  candidates <- unique(unlist(
    lapply(starting_directories, .comparison_ancestors),
    use.names = FALSE
  ))
  marker <- file.path(
    candidates,
    "software", "SegreCP", "DESCRIPTION"
  )
  hit <- which(file.exists(marker))
  if (!length(hit)) {
    stop(
      "Could not locate the manuscript repository root. Source this file from ",
      "its repository checkout or set the working directory to that root.",
      call. = FALSE
    )
  }
  normalizePath(candidates[[hit[[1L]]]], mustWork = TRUE)
}

#' Catalog the experiment-only comparison implementations
#'
#' The table is a routing record, not a statement that all third-party methods
#' have identical APIs or stopping rules.  Each benchmark worker fixes its own
#' rank, initialization, iteration budget, and timing boundary as documented in
#' the corresponding rebuttal-stage experiment.
comparison_method_catalog <- function() {
  data.frame(
    task = c(
      "CP decomposition", "CP decomposition", "CP decomposition",
      "tensor regression", "tensor regression", "tensor regression",
      "tensor regression", "tensor regression"
    ),
    method = c(
      "CP-ALS", "Orth-ALS", "ICO", "CP-ALS", "CP-RRR", "Tucker-RGD",
      "TT-IHT", "HOPLS"
    ),
    repository_worker = c(
      "software/SegreCP/reproducibility/comparison_methods.R",
      "rgn_benchmark_decomp/run_worker_decomp.R",
      "rgn_benchmark_decomp/run_worker_decomp.R",
      rep("rgn_benchmark/run_tensor_on_tensor_baseline_grid.R", 5L)
    ),
    package_api = FALSE,
    stringsAsFactors = FALSE
  )
}

.comparison_components_from_factors <- function(factors, weights) {
  lapply(seq_along(weights), function(component) {
    list(
      lambda = as.numeric(weights[[component]]),
      U = lapply(factors, function(A) as.numeric(A[, component]))
    )
  })
}

#' Run the experiment-specific CP-ALS comparison path
#'
#' This is the CP-ALS update used for the paired initialization-sensitivity
#' experiment.  The caller supplies a realized common initialization, so CP-ALS
#' can be compared with RGD/RGN without charging initialization time to any
#' optimizer.  Every row records the observed least-squares objective; when
#' simulated truth is supplied, whole-tensor relative reconstruction error is
#' recorded as an oracle diagnostic.
#'
#' @param observed Numeric array or `rTensor::Tensor`.
#' @param initialization Common CP component list, normally produced by CPCA or
#'   HOSVD in SegreCP.
#' @param iterations Positive number of complete ALS sweeps.
#' @param truth Optional simulated CP truth for error diagnostics.
#' @return A list containing final components and a sweep-aligned trajectory.
fit_comparison_cp_als <- function(
    observed,
    initialization,
    iterations = 10L,
    truth = NULL) {
  if (!requireNamespace("rTensor", quietly = TRUE)) {
    stop("The comparison CP-ALS path requires the 'rTensor' package.",
         call. = FALSE)
  }
  if (!requireNamespace("SegreCP", quietly = TRUE)) {
    stop("Install SegreCP to use its common component representation.",
         call. = FALSE)
  }
  observed <- SegreCP::as_segre_array(observed)
  if (length(iterations) != 1L || !is.numeric(iterations) ||
      !is.finite(iterations) || iterations != floor(iterations) ||
      iterations < 1L) {
    stop("`iterations` must be one positive integer.", call. = FALSE)
  }
  iterations <- as.integer(iterations)
  rank <- length(initialization)
  dims <- dim(observed)
  if (!rank || any(vapply(initialization, function(z) {
    !is.list(z) || length(z$U) != length(dims)
  }, logical(1L)))) {
    stop("`initialization` must be a nonempty SegreCP component list.",
         call. = FALSE)
  }

  factors <- lapply(seq_along(dims), function(mode) {
    do.call(cbind, lapply(initialization, function(z) z$U[[mode]]))
  })
  weights <- vapply(initialization, `[[`, numeric(1L), "lambda")
  unfolded <- lapply(seq_along(dims), function(mode) {
    rTensor::k_unfold(rTensor::as.tensor(observed), m = mode)@data
  })

  trajectory <- data.frame(
    iteration = 0:iterations,
    objective = NA_real_,
    relative_frobenius_error = NA_real_
  )
  record_state <- function(iteration) {
    components <- .comparison_components_from_factors(factors, weights)
    fitted <- SegreCP::reconstruct_cp(components)
    trajectory$objective[[iteration + 1L]] <<-
      0.5 * sum((fitted - observed)^2)
    if (!is.null(truth)) {
      trajectory$relative_frobenius_error[[iteration + 1L]] <<-
        SegreCP::cp_relative_error(components, truth)
    }
    components
  }

  components <- record_state(0L)
  for (iteration in seq_len(iterations)) {
    for (mode in seq_along(dims)) {
      gram <- Reduce("*", lapply(factors[-mode], crossprod))
      khatri_rao <- rTensor::khatri_rao_list(
        factors[-mode], reverse = TRUE
      )
      updated <- unfolded[[mode]] %*% khatri_rao %*% solve(gram)
      norms <- sqrt(colSums(updated * updated))
      if (any(!is.finite(norms)) || any(norms <= .Machine$double.eps)) {
        stop(
          "CP-ALS encountered a singular component; the comparison path ",
          "cannot continue from this initialization.",
          call. = FALSE
        )
      }
      factors[[mode]] <- sweep(updated, 2L, norms, "/")
      weights <- norms
    }
    components <- record_state(iteration)
  }

  list(
    method = "CP-ALS",
    components = components,
    trajectory = trajectory
  )
}

#' Launch the rebuttal-stage tensor-regression comparator grid
#'
#' Only comparison methods are requested.  RGD and RGN are intentionally absent
#' from the `--methods` argument; their experiment-facing entry point is the
#' package script.  Set `execute = FALSE` to inspect the exact command without
#' starting an expensive benchmark.
run_tensor_regression_comparisons <- function(
    extra_args = character(),
    execute = FALSE) {
  root <- comparison_repository_root()
  runner <- file.path(
    root, "rgn_benchmark", "run_tensor_on_tensor_baseline_grid.R"
  )
  methods <- paste(
    c("CP-ALS", "CP-RRR", "Tucker-RGD", "TT-IHT", "HOPLS"),
    collapse = ","
  )
  command <- c(
    runner,
    "--methods", methods,
    as.character(extra_args)
  )
  if (!isTRUE(execute)) {
    return(list(executable = file.path(R.home("bin"), "Rscript"),
                arguments = command))
  }
  status <- system2(file.path(R.home("bin"), "Rscript"), command)
  if (!identical(status, 0L)) {
    stop("The tensor-regression comparison worker failed with status ", status,
         ".", call. = FALSE)
  }
  invisible(status)
}
