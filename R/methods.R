.format_response <- function(value, output_dims) {
  n <- nrow(value)
  if (!length(output_dims)) return(as.numeric(value[, 1L]))
  if (length(output_dims) == 1L) return(unname(value))
  array(as.numeric(value), dim = c(n, output_dims))
}

#' Print SegreCP controls
#'
#' @param x A `segre_control` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print segre_control
#' @export
print.segre_control <- function(x, ...) {
  cat(sprintf("SegreCP control (%s profile)\n", x$profile))
  cat(sprintf("  scheme: %s; component order: %s\n",
              x$scheme, x$component_order))
  cat(sprintf("  max sweeps: %d; tolerance: %g; RGN ridge: %g\n",
              x$max_iter, x$tolerance, x$ridge))
  invisible(x)
}

#' Print a SegreCP fit
#'
#' @param x A fitted object from [segre_decompose()] or [segre_regress()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print segre_cp_fit
#' @export
print.segre_cp_fit <- function(x, ...) {
  label <- if (x$task == "decomposition") "CP decomposition" else
    "CP tensor regression"
  cat(sprintf("SegreCP %s fit\n", label))
  cat(sprintf("  method: %s; profile: %s; rank: %d\n",
              toupper(x$method), x$profile, x$rank))
  cat(sprintf("  dimensions: %s\n", paste(x$dims, collapse = " x ")))
  cat(sprintf("  sweeps: %d; stop: %s\n", x$iterations,
              x$stopping_reason))
  cat(sprintf("  objective: %.8g -> %.8g\n",
              x$objective[[1L]], utils::tail(x$objective, 1L)))
  invisible(x)
}

#' Summarize a SegreCP fit
#'
#' @param object A fitted object from [segre_decompose()] or [segre_regress()].
#' @param ... Unused.
#' @return A `summary.segre_cp_fit` list.
#' @export
summary.segre_cp_fit <- function(object, ...) {
  condition <- if (object$task == "regression" && object$method == "rgn") {
    z <- object$normal_condition[is.finite(object$normal_condition)]
    if (length(z)) c(min = min(z), median = stats::median(z), max = max(z)) else
      c(min = NA_real_, median = NA_real_, max = NA_real_)
  } else NULL
  structure(list(
    task = object$task,
    method = object$method,
    profile = object$profile,
    rank = object$rank,
    dims = object$dims,
    iterations = object$iterations,
    converged = object$converged,
    stopping_reason = object$stopping_reason,
    initial_objective = object$objective[[1L]],
    final_objective = utils::tail(object$objective, 1L),
    objective_ratio = utils::tail(object$objective, 1L) /
      max(object$objective[[1L]], .Machine$double.eps),
    normal_condition = condition
  ), class = "summary.segre_cp_fit")
}

#' Print a SegreCP fit summary
#'
#' @param x An object returned by [summary.segre_cp_fit()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print summary.segre_cp_fit
#' @export
print.summary.segre_cp_fit <- function(x, ...) {
  cat(sprintf("%s / %s (%s profile)\n",
              if (x$task == "decomposition") "CP decomposition" else
                "CP tensor regression",
              toupper(x$method), x$profile))
  cat(sprintf("  rank %d; dimensions %s\n", x$rank,
              paste(x$dims, collapse = " x ")))
  cat(sprintf("  objective %.8g -> %.8g after %d sweeps\n",
              x$initial_objective, x$final_objective, x$iterations))
  cat(sprintf("  stop: %s\n", x$stopping_reason))
  if (!is.null(x$normal_condition)) {
    cat(sprintf("  restricted-normal condition: min %.4g, median %.4g, max %.4g\n",
                x$normal_condition[[1L]], x$normal_condition[[2L]],
                x$normal_condition[[3L]]))
  }
  invisible(x)
}

#' Extract a fitted CP coefficient tensor
#'
#' @param object A `segre_cp_fit` object.
#' @param format A base R `"array"` or `"rTensor"` Tensor.
#' @param ... Unused.
#' @return The dense fitted tensor.
#' @export
coef.segre_cp_fit <- function(object, format = c("array", "rTensor"), ...) {
  reconstruct_cp(object$components, match.arg(format))
}

#' @export
as.array.segre_cp_fit <- function(x, ...) {
  reconstruct_cp(x$components, "array")
}

#' Convert to an rTensor Tensor
#'
#' @param x A base R tensor array or a fitted `segre_cp_fit` object.
#' @param ... Unused.
#' @return An `rTensor::Tensor`.
#' @export
as_rtensor <- function(x, ...) UseMethod("as_rtensor")

#' @export
as_rtensor.default <- function(x, ...) {
  rTensor::as.tensor(as_segre_array(x))
}

#' @export
as_rtensor.segre_cp_fit <- function(x, ...) {
  reconstruct_cp(x$components, "rTensor")
}

#' Predict from a SegreCP regression fit
#'
#' @param object A `segre_cp_regression` fit.
#' @param newdata New tensor covariates in any format accepted by
#'   [segre_regress()].
#' @param ... Unused.
#' @return A vector for scalar response, a matrix for a one-mode tensor response,
#'   or an array with dimensions `c(n, output_dims)` for multiple output modes.
#' @export
predict.segre_cp_regression <- function(object, newdata, ...) {
  design <- .as_design_matrix(newdata, object$input_dims)
  internal <- .public_to_internal_components(
    object$components, length(object$input_dims), length(object$output_dims)
  )
  split <- .split_internal_dims(c(object$output_dims, object$input_dims),
                                length(object$output_dims))
  value <- Reduce(`+`, lapply(internal, function(comp) {
    .apply_rank1_operator(comp, design$X, split)$prediction
  }))
  .format_response(value, object$output_dims)
}

#' @export
fitted.segre_cp_regression <- function(object, ...) {
  .format_response(object$fitted_values, object$output_dims)
}

#' @export
residuals.segre_cp_regression <- function(object, ...) {
  .format_response(object$residual_values, object$output_dims)
}

#' Extract an empirical convergence trajectory
#'
#' Returns the objective value at initialization and after every completed
#' sweep. If `truth` is supplied, the fit must have been run with
#' `keep_history = TRUE`; the function then also computes the relative
#' Frobenius error of the reconstructed CP tensor without forming a dense
#' tensor.
#'
#' The optional reconstruction error is an empirical whole-tensor metric. It
#' is not the permutation-matched component metric used in some theoretical
#' recursions in the manuscript. The CP-Gram identity can lose relative
#' accuracy through cancellation extremely close to exact recovery; use a
#' direct dense difference or a stable matched-component diagnostic when the
#' numerical floor itself is the estimand.
#'
#' @param object A fitted object from [segre_decompose()] or [segre_regress()].
#' @param truth Optional ground-truth CP component list or fitted object.
#' @return A data frame with one row for iteration zero and each completed
#'   sweep.
#' @examples
#' dims <- c(3, 4, 2)
#' factors <- lapply(dims, function(p) matrix(c(1, rep(0, p - 1)), p, 1))
#' truth <- cp_components(2, factors)
#' x <- reconstruct_cp(truth)
#' fit <- segre_decompose(
#'   x, 1, method = "rgn", init = "hosvd",
#'   control = segre_control(max_iter = 2), keep_history = TRUE
#' )
#' convergence_data(fit, truth)
#' @export
convergence_data <- function(object, truth = NULL) {
  if (!inherits(object, "segre_cp_fit")) {
    stop("`object` must be a SegreCP fitted object.", call. = FALSE)
  }
  objective <- as.numeric(object$objective)
  n_path <- length(objective)
  if (!n_path || any(!is.finite(objective))) {
    stop("The fitted objective trajectory is empty or nonfinite.",
         call. = FALSE)
  }
  initial_scale <- max(abs(objective[[1L]]), .Machine$double.eps)
  relative_change <- c(
    NA_real_,
    if (n_path > 1L) {
      (objective[-n_path] - objective[-1L]) /
        pmax(abs(objective[-n_path]), .Machine$double.eps)
    } else numeric()
  )
  out <- data.frame(
    iteration = seq.int(0L, n_path - 1L),
    objective = objective,
    objective_ratio = objective / initial_scale,
    relative_objective_decrease = relative_change,
    task = rep(object$task, n_path),
    method = rep(object$method, n_path),
    profile = rep(object$profile, n_path),
    stringsAsFactors = FALSE
  )
  if (!is.null(truth)) {
    history <- object$component_history
    if (is.null(history)) {
      stop("Refit with `keep_history = TRUE` to compute reconstruction error.",
           call. = FALSE)
    }
    if (length(history) != n_path) {
      stop("The component and objective histories have inconsistent lengths.",
           call. = FALSE)
    }
    out$relative_frobenius_error <- vapply(
      history, cp_relative_error, numeric(1L), truth = truth
    )
  }
  out
}
