#' Options for RGD and RGN
#'
#' Sets the number of iterations, the stopping rule, the RGD step size and the
#' order of the component updates for [segre_decompose()] and
#' [segre_regress()]. One iteration updates each of the CP components once,
#' in turn.
#'
#' The two profiles set defaults for the options left `NULL`:
#' * `"manuscript"`: the algorithms analyzed in the paper. Components are
#'   updated in the order 1, ..., r, RGD uses a fixed step, and `max_iter`
#'   iterations are run.
#' * `"practical"`: in each iteration the components are updated in
#'   decreasing order of their weights, the fit stops when the relative
#'   decrease of the objective is at most `tolerance = 1e-6`, and regression
#'   RGD uses backtracking.
#'
#' @param profile `"manuscript"` or `"practical"`.
#' @param max_iter Maximum number of iterations.
#' @param tolerance The fit stops when the objective decreases by a relative
#'   amount of at most `tolerance` in an iteration. Zero (the default of the
#'   manuscript profile) runs all `max_iter` iterations.
#' @param step_size RGD step size: one number, or one per iteration. `NULL`
#'   gives 0.5 for decomposition and `1 / design_variance` (that is,
#'   \eqn{\sigma_X^{-2}}) for regression, the fixed steps covered by the
#'   convergence results of the paper. RGN does not use a step size.
#' @param step_rule How RGD uses `step_size`: `"fixed"`; `"diminishing"`
#'   (`step_size / sqrt(t)` in iteration `t`); `"backtracking"` (start from
#'   `step_size` and shrink it until the Armijo condition holds); or, for
#'   regression only, `"exact_quadratic"` (the step that minimizes the
#'   objective along the projected gradient before the retraction, restricted
#'   to `[1e-10, 100]`). `NULL` gives `"backtracking"` for regression with the
#'   practical profile and `"fixed"` otherwise.
#' @param component_order `"fixed"` (1, ..., r) or `"decreasing_weight"`
#'   (largest current weight first). `NULL` uses the profile default.
#' @param ridge Nonnegative ridge added to the RGN normal equation in
#'   regression. The default 0 gives RGN as in the paper.
#' @param normalize_regression `"sqrt_n"` divides the design and the response
#'   by `sqrt(n)`, so that the objective is an average over observations.
#'   `"none"` uses the data as given, for data already divided by `sqrt(n)`.
#' @param trace If `TRUE`, report the objective after every iteration.
#' @param armijo,backtrack_factor,max_backtracks Backtracking: the search
#'   stops at a step `alpha` that decreases the objective by at least
#'   `armijo * alpha` times the squared norm of the projected gradient, and
#'   otherwise multiplies `alpha` by `backtrack_factor`, for at most
#'   `max_backtracks` trials. The update uses the trial step with the smallest
#'   objective, and leaves the component unchanged if no trial decreases it.
#' @return An object of class `segre_control`.
#' @examples
#' segre_control()
#' segre_control(profile = "practical", max_iter = 200)
#' segre_control(max_iter = 50, step_size = 0.3)
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
  if (is.null(component_order)) {
    component_order <- if (profile == "manuscript") "fixed" else
      "decreasing_weight"
  }
  component_order <- match.arg(component_order,
                               c("fixed", "decreasing_weight"))
  if (is.null(ridge)) ridge <- 0
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be one nonnegative finite number.", call. = FALSE)
  }
  if (!is.null(step_size) && (!is.numeric(step_size) || !length(step_size) ||
      any(!is.finite(step_size)) || any(step_size <= 0))) {
    stop("`step_size` must contain positive finite values.", call. = FALSE)
  }
  if (!is.null(step_rule)) {
    step_rule <- match.arg(step_rule,
      c("fixed", "diminishing", "backtracking", "exact_quadratic"))
  }
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("`ridge` must be one nonnegative finite number.", call. = FALSE)
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
    step_size = if (is.null(step_size)) NULL else as.numeric(step_size),
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

.resolve_control <- function(control, task = c("decomposition", "regression"),
                             design_variance = 1) {
  task <- match.arg(task)
  if (is.null(control)) control <- segre_control()
  if (!inherits(control, "segre_control")) {
    stop("`control` must be created by `segre_control()`.", call. = FALSE)
  }
  out <- control
  if (is.null(out$step_size)) {
    out$step_size <- if (task == "regression") 1 / design_variance else 0.5
  }
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
