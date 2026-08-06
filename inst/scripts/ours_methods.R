#!/usr/bin/env Rscript

# Experiment-facing entry point for the methods proposed in the manuscript.
#
# This script intentionally exposes only three pieces of the SegreCP workflow:
#   1. CPCA, HOSVD, random, or user-supplied initialization;
#   2. Riemannian gradient descent (RGD); and
#   3. Riemannian Gauss--Newton (RGN).
#
# The numerical algorithms themselves live in the modular R/C++ package
# sources.  Keeping this thin entry point separate from the implementation has
# two advantages: experiments share one documented interface, while package
# internals remain testable and maintainable.  No competing estimator is
# implemented or loaded here.

.require_segrecp <- function() {
  if (!requireNamespace("SegreCP", quietly = TRUE)) {
    stop(
      "Install SegreCP before sourcing inst/scripts/ours_methods.R.",
      call. = FALSE
    )
  }
}

.ours_method <- function(method) {
  method <- toupper(as.character(method))
  if (length(method) != 1L || !method %in% c("RGD", "RGN")) {
    stop("`method` must be either 'RGD' or 'RGN'.", call. = FALSE)
  }
  tolower(method)
}

.ours_initialization <- function(initialization) {
  # Lists and functions are deliberately passed through: the public package API
  # supports both user-supplied component lists and custom initializer functions.
  if (!is.character(initialization)) return(initialization)
  initialization <- tolower(initialization)
  if (length(initialization) != 1L ||
      !initialization %in% c("cpca", "hosvd", "random")) {
    stop(
      "Character `initialization` must be 'cpca', 'hosvd', or 'random'.",
      call. = FALSE
    )
  }
  initialization
}

#' Initialize a CP model for an authors' method
#'
#' @param x Numeric base-R array or `rTensor::Tensor`.
#' @param rank Positive target CP rank.
#' @param initialization `"cpca"`, `"hosvd"`, `"random"`, a component list,
#'   or a custom initializer function.
#' @param cpca_split Optional one-based mode subset for CPCA.
#' @param seed Optional local seed used only by random initialization.
#' @param ... Arguments forwarded to a custom initializer function.
#' @return The complete initialization record returned by
#'   `SegreCP::cp_initialize()`.
initialize_ours <- function(
    x,
    rank,
    initialization = "cpca",
    cpca_split = NULL,
    seed = NULL,
    ...) {
  .require_segrecp()
  SegreCP::cp_initialize(
    x = x,
    rank = rank,
    method = .ours_initialization(initialization),
    split = cpca_split,
    seed = seed,
    ...
  )
}

#' Fit the authors' RGD or RGN CP decomposition
#'
#' @param x Observed tensor as a base-R array or `rTensor::Tensor`.
#' @param rank Positive target CP rank.
#' @param method `"RGD"` or `"RGN"`.
#' @param initialization Any initialization supported by `initialize_ours()`.
#' @param cpca_split Optional CPCA mode subset.
#' @param init_args Named arguments for a custom initializer.
#' @param seed Optional local random-initialization seed.
#' @param control A `SegreCP::segre_control()` object.  When `NULL`, the exact
#'   manuscript profile is used with `max_iter` complete sweeps.
#' @param max_iter Number of sweeps used only when `control` is `NULL`.
#' @param keep_history Retain every completed-sweep iterate for convergence
#'   diagnostics.
#' @return A `segre_cp_decomposition` fit.
fit_ours_decomposition <- function(
    x,
    rank,
    method = "RGN",
    initialization = "cpca",
    cpca_split = NULL,
    init_args = list(),
    seed = NULL,
    control = NULL,
    max_iter = 10L,
    keep_history = TRUE) {
  .require_segrecp()
  if (is.null(control)) {
    control <- SegreCP::segre_control(
      profile = "manuscript",
      max_iter = max_iter,
      tolerance = 0
    )
  }
  SegreCP::segre_decompose(
    x = x,
    rank = rank,
    method = .ours_method(method),
    init = .ours_initialization(initialization),
    init_split = cpca_split,
    init_args = init_args,
    seed = seed,
    control = control,
    keep_history = keep_history
  )
}

#' Fit the authors' RGD or RGN tensor regression
#'
#' This wrapper covers both scalar-on-tensor and tensor-on-tensor responses.
#' Public coefficient modes always follow the manuscript convention: input
#' modes first, followed by output modes.
#'
#' @param X Tensor covariates in any format accepted by `SegreCP::segre_regress`.
#' @param Y Scalar or tensor responses in any supported format.
#' @param rank Positive target CP rank.
#' @param input_dims,output_dims Public coefficient-mode dimensions.
#' @param method `"RGD"` or `"RGN"`.
#' @param initialization Any initialization supported by `initialize_ours()`.
#' @param cpca_split Optional CPCA mode subset in public mode order.
#' @param init_args Named arguments for a custom initializer.
#' @param seed Optional local random-initialization seed.
#' @param design_variance Positive entry variance used in the spectral moment.
#' @param control A `SegreCP::segre_control()` object.  When `NULL`, the exact
#'   manuscript profile is used with `max_iter` complete sweeps.
#' @param max_iter Number of sweeps used only when `control` is `NULL`.
#' @param keep_history Retain every completed-sweep iterate.
#' @return A `segre_cp_regression` fit.
fit_ours_regression <- function(
    X,
    Y,
    rank,
    input_dims = NULL,
    output_dims = NULL,
    method = "RGN",
    initialization = "cpca",
    cpca_split = NULL,
    init_args = list(),
    seed = NULL,
    design_variance = 1,
    control = NULL,
    max_iter = 10L,
    keep_history = TRUE) {
  .require_segrecp()
  if (is.null(control)) {
    control <- SegreCP::segre_control(
      profile = "manuscript",
      max_iter = max_iter,
      tolerance = 0
    )
  }
  SegreCP::segre_regress(
    X = X,
    Y = Y,
    rank = rank,
    input_dims = input_dims,
    output_dims = output_dims,
    method = .ours_method(method),
    init = .ours_initialization(initialization),
    init_split = cpca_split,
    init_args = init_args,
    seed = seed,
    design_variance = design_variance,
    control = control,
    keep_history = keep_history
  )
}

#' Return an auditable optimization trajectory
#'
#' @param fit Fitted decomposition or regression object from this script.
#' @param truth Optional simulated truth in the representation accepted by
#'   `SegreCP::convergence_data()`.
#' @return A data frame aligned to initialization and completed sweeps.
ours_convergence <- function(fit, truth = NULL) {
  .require_segrecp()
  SegreCP::convergence_data(fit, truth = truth)
}

