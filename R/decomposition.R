#' Fit a CP decomposition with RGD or RGN
#'
#' For decomposition the observation operator is the identity. Consequently,
#' the RGN tangent normal equation is the identity on the current tangent space,
#' and RGN is exactly the unit-step RGD update. Every accepted component update
#' uses the exact rank-(1,...,1) T-HOSVD retraction of the rank-one-plus-tangent
#' proposal.
#'
#' @param x Observed tensor as a base R numeric array or `rTensor::Tensor`.
#' @param rank Positive target CP rank.
#' @param method `"rgd"` or `"rgn"`.
#' @param init `"cpca"`, `"hosvd"`, `"random"`, a custom component list, or a
#'   custom function; see [cp_initialize()].
#' @param init_split Optional CPCA mode split.
#' @param init_args Named list passed to [cp_initialize()].
#' @param seed Optional local seed for random initialization.
#' @param control Controls created by [segre_control()].
#' @param keep_history If `TRUE`, retain the CP components at initialization and
#'   after every completed sweep. This is useful for convergence diagnostics and
#'   costs `rank * sum(dims)` numbers per stored sweep.
#' @return A fitted `segre_cp_decomposition` object. Its `components` field uses
#'   unit factors and nonnegative weights.
#' @examples
#' set.seed(10)
#' u <- lapply(c(5, 6, 4), function(p) {
#'   z <- rnorm(p); z / sqrt(sum(z^2))
#' })
#' signal <- array(3 * Reduce(kronecker, rev(u)), c(5, 6, 4))
#' fit <- segre_decompose(signal, rank = 1, method = "rgn", init = "hosvd",
#'                        control = segre_control(max_iter = 3))
#' reconstruct_cp(fit)
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
  if (length(dim(x)) < 2L) stop("CP decomposition requires at least two modes.")
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

  initialization <- .resolve_initialization(
    x, rank, init, init_split, seed, init_args
  )
  components <- .validate_components(initialization$comps, dim(x), rank)
  component_history <- if (keep_history) list(components) else NULL
  caches <- lapply(components, function(comp) .rank1_vector(comp$lambda, comp$U))
  residual <- array(Reduce(`+`, caches) - as.numeric(x), dim = dim(x))

  objective <- numeric(control$max_iter + 1L)
  objective[[1L]] <- 0.5 * sum(residual * residual)
  accepted_steps <- matrix(NA_real_, control$max_iter, rank)
  converged <- FALSE
  stopping_reason <- "max_iter"
  iterations <- control$max_iter

  for (iteration in seq_len(control$max_iter)) {
    order_i <- .component_update_order(components, control)
    for (i in order_i) {
      alpha <- if (method == "rgn") 1 else .step_at(control, iteration)
      old_comp <- components[[i]]
      old_cache <- caches[[i]]

      if (method == "rgd" && control$step_rule == "backtracking") {
        loss_old <- 0.5 * sum(residual * residual)
        projected <- tangent_proj_cpp(residual, old_comp$U)
        gradient_norm_sq <- sum(projected * projected)
        best <- list(comp = old_comp, loss = loss_old, step = 0)
        for (attempt in seq_len(control$max_backtracks)) {
          proposal <- rgd_step_fused_cpp(
            residual, old_comp$U, old_comp$lambda, alpha, 0L
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
      } else {
        proposal <- rgd_step_fused_cpp(
          residual, old_comp$U, old_comp$lambda, alpha, 0L
        )
        new_comp <- .canonicalize_component(
          list(lambda = proposal$lambda, U = proposal$U)
        )
      }

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

    objective[[iteration + 1L]] <- 0.5 * sum(residual * residual)
    if (keep_history) {
      component_history[[iteration + 1L]] <- components
    }
    if (control$trace) {
      message(sprintf("sweep %d: objective = %.8g", iteration,
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
