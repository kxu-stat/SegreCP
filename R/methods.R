# Print, summary, extraction and prediction methods.

# n x q matrix of responses -> vector, matrix or array, as the input.
.format_response <- function(value, output_dims) {
  n <- nrow(value)
  if (!length(output_dims)) return(as.numeric(value[, 1L]))
  if (length(output_dims) == 1L) return(unname(value))
  array(as.numeric(value), dim = c(n, output_dims))
}

#' Print SegreCP options
#'
#' @param x A `segre_control` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print segre_control
#' @export
print.segre_control <- function(x, ...) {
  stopping <- if (x$tolerance > 0) {
    sprintf("stop when the relative decrease is at most %g", x$tolerance)
  } else "no early stopping"
  order <- if (x$component_order == "fixed") "fixed order" else
    "decreasing weight"
  step <- if (is.null(x$step_size)) {
    "0.5 (decomposition), 1 / design_variance (regression)"
  } else if (length(x$step_size) == 1L) {
    format(x$step_size)
  } else sprintf("schedule of length %d", length(x$step_size))
  rule <- if (!is.null(x$step_rule)) x$step_rule else
    if (x$profile == "practical") {
      "fixed (decomposition), backtracking (regression)"
    } else "fixed"
  cat(sprintf("SegreCP options (%s profile)\n", x$profile))
  cat(sprintf("  iterations: at most %d, %s\n", x$max_iter, stopping))
  cat(sprintf("  component updates: one at a time, %s\n", order))
  cat(sprintf("  RGD step: %s\n", step))
  cat(sprintf("  RGD step rule: %s\n", rule))
  cat(sprintf("  RGN ridge: %g\n", x$ridge))
  cat(sprintf("  regression normalization: %s\n", x$normalize_regression))
  invisible(x)
}

.fit_title <- function(x) {
  task <- if (x$task == "decomposition") "CP decomposition" else
    "CP tensor regression"
  sprintf("%s by %s, rank %d (%s profile)", task, toupper(x$method), x$rank,
          x$profile)
}

.fit_dimensions <- function(x) {
  if (x$task == "decomposition") {
    return(sprintf("dimensions: %s", paste(x$dims, collapse = " x ")))
  }
  input <- paste(x$input_dims, collapse = " x ")
  if (!length(x$output_dims)) {
    return(sprintf("coefficient dimensions: %s", input))
  }
  sprintf("coefficient dimensions: %s (input) x %s (output)", input,
          paste(x$output_dims, collapse = " x "))
}

.fit_stop <- function(x) {
  if (x$stopping_reason == "relative_objective_change") {
    "converged: relative decrease of the objective at most the tolerance"
  } else "max_iter reached"
}

.fit_weights <- function(x) {
  vapply(x$components, `[[`, numeric(1L), "lambda")
}

#' Print a SegreCP fit
#'
#' @param x A fit from [segre_decompose()] or [segre_regress()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print segre_cp_fit
#' @export
print.segre_cp_fit <- function(x, ...) {
  cat(.fit_title(x), "\n", sep = "")
  cat("  ", .fit_dimensions(x), "\n", sep = "")
  cat(sprintf("  iterations: %d (%s)\n", x$iterations, .fit_stop(x)))
  cat(sprintf("  objective: %s -> %s\n", format(x$objective[[1L]], digits = 6),
              format(utils::tail(x$objective, 1L), digits = 6)))
  cat(sprintf("  weights: %s\n",
              paste(format(.fit_weights(x), digits = 4), collapse = "  ")))
  invisible(x)
}

#' Summarize a SegreCP fit
#'
#' @param object A fit from [segre_decompose()] or [segre_regress()].
#' @param x An object returned by `summary()`.
#' @param ... Unused.
#' @return `summary()` returns a list of class `summary.segre_cp_fit` with the
#'   task, method, profile, rank, dimensions, component weights, number of
#'   iterations, stopping reason, the objective at the start and at the end and
#'   their ratio, and for regression RGN the range of the estimated condition
#'   numbers of the normal equations. Its `print()` method returns `x`
#'   invisibly.
#' @examples
#' set.seed(1)
#' x <- array(rnorm(60), c(5, 4, 3))
#' summary(segre_decompose(x, rank = 1, control = segre_control(max_iter = 5)))
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
    input_dims = object$input_dims,
    output_dims = object$output_dims,
    weights = .fit_weights(object),
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

#' @rdname summary.segre_cp_fit
#' @method print summary.segre_cp_fit
#' @export
print.summary.segre_cp_fit <- function(x, ...) {
  cat(.fit_title(x), "\n", sep = "")
  cat("  ", .fit_dimensions(x), "\n", sep = "")
  cat(sprintf("  iterations: %d (%s)\n", x$iterations, .fit_stop(x)))
  cat(sprintf("  objective: %s -> %s (ratio %s)\n",
              format(x$initial_objective, digits = 6),
              format(x$final_objective, digits = 6),
              format(x$objective_ratio, digits = 4)))
  cat(sprintf("  weights: %s\n",
              paste(format(x$weights, digits = 6), collapse = "  ")))
  if (!is.null(x$normal_condition)) {
    cat(sprintf(paste("  condition number of the RGN normal equations:",
                      "min %.4g, median %.4g, max %.4g\n"),
                x$normal_condition[[1L]], x$normal_condition[[2L]],
                x$normal_condition[[3L]]))
  }
  invisible(x)
}

#' Estimated tensor of a SegreCP fit
#'
#' `coef()` and `as.array()` return the estimated tensor as a dense array: the
#' low-rank tensor for [segre_decompose()] and the coefficient tensor, with
#' dimensions `c(input_dims, output_dims)`, for [segre_regress()]. The CP
#' components themselves are in `object$components`.
#'
#' @param object,x A fit from [segre_decompose()] or [segre_regress()].
#' @param format `"array"` or `"rTensor"`.
#' @param ... Unused.
#' @return A dense array, or an `rTensor::Tensor` for `format = "rTensor"`.
#' @export
coef.segre_cp_fit <- function(object, format = c("array", "rTensor"), ...) {
  reconstruct_cp(object$components, match.arg(format))
}

#' @rdname coef.segre_cp_fit
#' @export
as.array.segre_cp_fit <- function(x, ...) {
  reconstruct_cp(x$components, "array")
}

#' Convert to an rTensor Tensor
#'
#' @param x A numeric array, or a fit from [segre_decompose()] or
#'   [segre_regress()] (converted to its estimated tensor).
#' @param ... Unused.
#' @return An `rTensor::Tensor`.
#' @seealso [as_segre_array()] for the conversion in the other direction.
#' @examples
#' as_rtensor(array(1:8, c(2, 2, 2)))
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

#' Predictions, fitted values and residuals of a SegreCP regression
#'
#' @param object A fit from [segre_regress()].
#' @param newdata Design for the new observations, in any format accepted by
#'   [segre_regress()].
#' @param ... Unused.
#' @return A vector for a scalar response, a matrix for a response with one
#'   mode, and otherwise an array with dimensions `c(n, output_dims)`.
#'   `fitted()` and `residuals()` refer to the data used for the fit, on their
#'   original scale.
#' @examples
#' set.seed(2)
#' X <- matrix(rnorm(60 * 6), 60)
#' y <- as.numeric(X %*% rep(c(1, -1), 3)) + rnorm(60, sd = 0.1)
#' fit <- segre_regress(X, y, rank = 1, input_dims = c(2, 3),
#'                      control = segre_control(max_iter = 5))
#' predict(fit, X[1:2, ])
#' head(residuals(fit))
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

#' @rdname predict.segre_cp_regression
#' @export
fitted.segre_cp_regression <- function(object, ...) {
  .format_response(object$fitted_values, object$output_dims)
}

#' @rdname predict.segre_cp_regression
#' @export
residuals.segre_cp_regression <- function(object, ...) {
  .format_response(object$residual_values, object$output_dims)
}

#' Objective and estimation error at each iteration
#'
#' Returns the objective at the start and after each iteration of a fit. If
#' `truth` is given, the fit must have been run with `keep_history = TRUE`, and
#' the relative error of the estimated tensor after each iteration is added
#' (see [cp_relative_error()]). This error is computed from inner products of
#' the factors, so values below about `1e-7` are not accurate. For such small
#' errors, compare the dense tensors from [reconstruct_cp()].
#'
#' @param object A fit from [segre_decompose()] or [segre_regress()].
#' @param truth Optional true CP components (see [cp_components()]).
#' @return A data frame with one row per iteration (0 is the start): the
#'   objective, its ratio to the starting value, its relative decrease from the
#'   previous iteration, and if `truth` is given the relative Frobenius error
#'   `relative_frobenius_error`.
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
