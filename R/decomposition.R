#' CP decomposition by RGD or RGN
#'
#' Estimates a tensor of CP rank `rank` from a noisy observation
#' \deqn{Y = \sum_{i=1}^r \lambda_i\, u_{i1} \otimes \cdots \otimes u_{id} + E}{Y = sum_i lambda_i u_i1 o ... o u_id + E}
#' by minimizing \eqn{\frac12 \|\sum_i T_i - Y\|_F^2}{(1/2) ||sum_i T_i - Y||_F^2}
#' over rank-one tensors \eqn{T_1, \ldots, T_r}{T_1, ..., T_r} with Riemannian
#' gradient descent (RGD) or Riemannian Gauss--Newton (RGN) on the Segre
#' manifold. One iteration updates each component once, in turn, with the
#' others at their current values: the gradient is projected onto the tangent
#' space at the component, and the updated tensor is retracted to a rank-one
#' tensor by the rank-(1, ..., 1) T-HOSVD. For decomposition, RGN is RGD with
#' step size 1.
#'
#' @param x Observed tensor: a numeric array with at least two modes or an
#'   `rTensor::Tensor`.
#' @param rank CP rank.
#' @param method `"rgd"` or `"rgn"`.
#' @param init Start: `"cpca"`, `"hosvd"` or `"random"` (see
#'   [cp_initialize()]), a function called as `init(x, rank, ...)`, or a list
#'   of CP components (see [cp_components()]).
#' @param init_split For `init = "cpca"`, the modes that index the rows of the
#'   unfolding. By default the split that makes the unfolding closest to
#'   square.
#' @param init_args A named list of further arguments for [cp_initialize()] or
#'   for the function `init`.
#' @param seed Seed for `init = "random"`. The random number stream of the
#'   session is left unchanged.
#' @param control Options from [segre_control()]. The default RGD step size is
#'   0.5.
#' @param keep_history If `TRUE`, the components after every iteration are
#'   kept in `component_history`, for example to compute estimation errors with
#'   [convergence_data()]. This stores `rank * sum(dim(x))` numbers per
#'   iteration.
#' @return An object of class `c("segre_cp_decomposition", "segre_cp_fit")`: a
#'   list with elements including
#'   * `components`: the estimated CP components, a list of `list(lambda, U)`
#'     with positive weights `lambda` and unit factors `U`;
#'   * `objective`: the objective at the start and after each iteration;
#'   * `iterations`, `converged`, `stopping_reason`: the number of iterations
#'     run and why the fit stopped;
#'   * `component_history`: if `keep_history = TRUE`, the components at the
#'     start and after each iteration;
#'   * `accepted_steps`: the step size of each component update;
#'   * `initialization`, `control`, `call`: the start and the options used.
#'
#'   [reconstruct_cp()], [coef()], [summary()] and [convergence_data()] apply
#'   to the result. `cp_decompose()` is an alias of `segre_decompose()`.
#' @seealso [segre_regress()] for tensor regression and [segre_control()] for
#'   the options.
#' @examples
#' # A rank-two 8 x 7 x 6 tensor observed with noise
#' set.seed(10)
#' factors <- lapply(c(8, 7, 6), function(p) qr.Q(qr(matrix(rnorm(2 * p), p))))
#' truth <- cp_components(c(5, 3), factors)
#' x <- reconstruct_cp(truth) + array(rnorm(8 * 7 * 6, sd = 0.05), c(8, 7, 6))
#' fit <- segre_decompose(x, rank = 2, method = "rgn",
#'                        control = segre_control(max_iter = 10))
#' fit
#' cp_relative_error(fit, truth)
#' @export
segre_decompose <- function(
    x, rank,
    method = c("rgd", "rgn"),
    init = c("cpca", "hosvd", "random"),
    init_split = NULL,
    init_args = list(),
    seed = NULL,
    control = segre_control(),
    keep_history = FALSE) {
  call <- match.call()
  input_class <- class(x)
  x <- as_segre_array(x)
  if (length(dim(x)) < 2L) {
    stop("CP decomposition requires at least two modes.", call. = FALSE)
  }
  if (sqrt(sum(x * x)) <= .Machine$double.eps) {
    stop("A zero observed tensor has no nonzero Segre initialization.",
         call. = FALSE)
  }
  rank <- .as_whole_number(rank, "rank")
  keep_history <- .as_flag(keep_history, "keep_history")
  method <- match.arg(method)
  control <- .resolve_control(control, "decomposition")
  if (control$step_rule == "exact_quadratic" && method == "rgd") {
    stop("`exact_quadratic` is a regression RGD option; use a fixed step here.",
         call. = FALSE)
  }

  initialization <- if (is.list(init)) {
    .component_start(init, dim(x), rank, init_args)
  } else .resolve_initialization(x, rank, init, init_split, seed, init_args)
  components <- .validate_components(initialization$comps, dim(x), rank)
  component_history <- if (keep_history) list(components) else NULL
  state <- initialize_decomposition_state_cpp(x, components)
  caches <- state$caches
  residual <- state$residual

  objective <- numeric(control$max_iter + 1L)
  objective[[1L]] <- state$objective
  accepted_steps <- matrix(NA_real_, control$max_iter, rank)
  converged <- FALSE
  stopping_reason <- "max_iter"
  iterations <- control$max_iter

  for (iteration in seq_len(control$max_iter)) {
    order_i <- .component_update_order(components, control)
    batch <- NULL
    if (method != "rgd" || control$step_rule != "backtracking") {
      # All components in one compiled sweep (src/kernels.cpp).
      alpha <- if (method == "rgn") 1 else .step_at(control, iteration)
      batch <- rgd_sweep_gs_batch_cpp(
        residual, caches, components, as.integer(order_i), alpha,
        objective[[iteration]]
      )
      components <- batch$components
      accepted_steps[iteration, order_i] <- alpha
    } else for (i in order_i) {
      # RGD with backtracking: one component at a time.
      alpha <- .step_at(control, iteration)
      old_comp <- components[[i]]
      old_cache <- caches[[i]]

      loss_old <- 0.5 * sum(residual * residual)
      projected <- tangent_proj_cpp(residual, old_comp$U)
      gradient_norm_sq <- sum(projected * projected)
      best <- list(comp = old_comp, loss = loss_old, step = 0)
      for (attempt in seq_len(control$max_backtracks)) {
        proposal <- rgd_step_fused_cpp(
          residual, old_comp$U, old_comp$lambda, alpha
        )
        new_comp <- .canonicalize_component(
          list(lambda = proposal$lambda, U = proposal$U)
        )
        new_cache <- .rank1_vector(new_comp$lambda, new_comp$U)
        candidate <- residual + array(new_cache - old_cache, dim = dim(x))
        loss_new <- 0.5 * sum(candidate * candidate)
        if (loss_new < best$loss) {
          best <- list(comp = new_comp, loss = loss_new, step = alpha)
        }
        if (loss_new <= loss_old - control$armijo * alpha * gradient_norm_sq) {
          break
        }
        alpha <- alpha * control$backtrack_factor
      }
      new_comp <- best$comp
      alpha <- best$step

      if (!is.finite(new_comp$lambda) ||
          new_comp$lambda <= .Machine$double.eps) {
        stop(sprintf("Component %d retracted to a zero/nonfinite tensor.", i),
             call. = FALSE)
      }
      update_rank1_cache_inplace_cpp(
        residual, caches[[i]], new_comp$lambda, new_comp$U
      )
      components[[i]] <- new_comp
      accepted_steps[iteration, i] <- alpha
    }

    objective[[iteration + 1L]] <- if (is.null(batch))
      0.5 * sum(residual * residual) else batch$objective
    if (keep_history) {
      component_history[[iteration + 1L]] <- components
    }
    if (control$trace) {
      message(sprintf("iteration %d: objective = %.8g", iteration,
                      objective[[iteration + 1L]]))
    }
    if (control$tolerance > 0) {
      previous <- objective[[iteration]]
      current <- objective[[iteration + 1L]]
      relative <- (previous - current) / (abs(previous) + .Machine$double.eps)
      if (relative >= 0 && relative <= control$tolerance) {
        converged <- TRUE
        stopping_reason <- "relative_objective_change"
        iterations <- iteration
        break
      }
    }
  }

  structure(list(
    call = call,
    task = "decomposition",
    method = method,
    profile = control$profile,
    rank = rank,
    dims = as.integer(dim(x)),
    components = components,
    initialization = initialization[names(initialization) != "comps"],
    objective = objective[seq_len(iterations + 1L)],
    component_history = component_history,
    accepted_steps = accepted_steps[seq_len(iterations), , drop = FALSE],
    iterations = iterations,
    converged = converged,
    stopping_reason = stopping_reason,
    control = control,
    input_class = input_class
  ), class = c("segre_cp_decomposition", "segre_cp_fit"))
}

#' @rdname segre_decompose
#' @export
cp_decompose <- segre_decompose
