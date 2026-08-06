#' Control RGD and RGN fitting
#'
#' The `"manuscript"` profile follows the displayed algorithms: sequential
#' Gauss--Seidel residual refresh, fixed component order, fixed RGD step size,
#' exact structured T-HOSVD retraction, a fixed iteration cap, and no RGN ridge.
#' The `"practical"` profile enables weight-based component ordering, an
#' observable loss-change stopping rule, regression RGD backtracking, and a
#' small explicitly recorded RGN ridge.
#'
#' @param profile `"manuscript"` or `"practical"`.
#' @param max_iter Maximum number of complete component sweeps.
#' @param tolerance Relative loss-change tolerance. Zero disables early stopping.
#' @param step_size RGD step size or schedule. `NULL` selects 0.5.
#' @param step_rule `"fixed"`, `"diminishing"`, `"backtracking"`, or
#'   `"exact_quadratic"`. `NULL` uses a profile/task-specific default.
#' @param component_order `"fixed"` or `"decreasing_weight"`. `NULL` uses the
#'   profile default.
#' @param ridge Nonnegative RGN ridge. A positive value defines stabilized,
#'   regularized RGN rather than the manuscript's unregularized normal equation.
#' @param normalize_regression If `"sqrt_n"`, divide both the design and response
#'   by `sqrt(n)` before fitting. Use `"none"` for data already expressed in the
#'   normalized-operator convention.
#' @param trace Print one line per completed sweep.
#' @param armijo,backtrack_factor,max_backtracks Backtracking controls.
#' @return A `segre_control` object.
#' @export
segre_control <- function(
    profile = c("manuscript", "practical"),
    max_iter = 100L,
    tolerance = NULL,
    step_size = NULL,
    step_rule = NULL,
    component_order = NULL,
    ridge = NULL,
    normalize_regression = c("sqrt_n", "none"),
    trace = FALSE,
    armijo = 1e-4,
    backtrack_factor = 0.5,
    max_backtracks = 20L) {
  profile <- match.arg(profile)
  normalize_regression <- match.arg(normalize_regression)
  max_iter <- .as_whole_number(max_iter, "max_iter")
  if (is.null(tolerance)) tolerance <- if (profile == "manuscript") 0 else 1e-6
  if (is.null(step_size)) step_size <- 0.5
  if (is.null(component_order)) {
    component_order <- if (profile == "manuscript") "fixed" else
      "decreasing_weight"
  }
  component_order <- match.arg(component_order,
                               c("fixed", "decreasing_weight"))
  if (is.null(ridge)) ridge <- if (profile == "manuscript") 0 else 1e-8
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one nonnegative finite number.", call. = FALSE)
  }
  if (!is.numeric(step_size) || !length(step_size) || any(!is.finite(step_size)) ||
      any(step_size <= 0)) {
    stop("`step_size` must contain positive finite values.", call. = FALSE)
  }
  if (!is.null(step_rule)) {
    step_rule <- match.arg(step_rule,
      c("fixed", "diminishing", "backtracking", "exact_quadratic"))
  }
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("`ridge` must be one nonnegative finite number.", call. = FALSE)
  }
  if (profile == "manuscript" && ridge != 0) {
    warning("A positive ridge changes manuscript RGN into regularized RGN.",
            call. = FALSE)
  }
  if (length(armijo) != 1L || !is.finite(armijo) || armijo <= 0 || armijo >= 1) {
    stop("`armijo` must lie strictly between zero and one.", call. = FALSE)
  }
  if (length(backtrack_factor) != 1L || !is.finite(backtrack_factor) ||
      backtrack_factor <= 0 || backtrack_factor >= 1) {
    stop("`backtrack_factor` must lie strictly between zero and one.",
         call. = FALSE)
  }
  max_backtracks <- .as_whole_number(max_backtracks, "max_backtracks")
  structure(list(
    profile = profile,
    scheme = "sequential",
    max_iter = max_iter,
    tolerance = as.numeric(tolerance),
    step_size = as.numeric(step_size),
    step_rule = step_rule,
    component_order = component_order,
    ridge = as.numeric(ridge),
    normalize_regression = normalize_regression,
    trace = isTRUE(trace),
    armijo = armijo,
    backtrack_factor = backtrack_factor,
    max_backtracks = max_backtracks
  ), class = "segre_control")
}

.resolve_control <- function(control, task = c("decomposition", "regression")) {
  task <- match.arg(task)
  if (is.null(control)) control <- segre_control()
  if (!inherits(control, "segre_control")) {
    stop("`control` must be created by `segre_control()`.", call. = FALSE)
  }
  out <- control
  if (is.null(out$step_rule)) {
    out$step_rule <- if (task == "regression" && out$profile == "practical") {
      "backtracking"
    } else "fixed"
  }
  if (length(out$step_size) > 1L && length(out$step_size) < out$max_iter) {
    stop("A step-size schedule must have at least `max_iter` entries.",
         call. = FALSE)
  }
  out
}

.component_update_order <- function(components, control) {
  if (control$component_order == "fixed") return(seq_along(components))
  order(vapply(components, `[[`, numeric(1L), "lambda"), decreasing = TRUE)
}

.step_at <- function(control, iteration) {
  step <- if (length(control$step_size) == 1L) control$step_size else
    control$step_size[[iteration]]
  if (control$step_rule == "diminishing") step <- step / sqrt(iteration)
  step
}
