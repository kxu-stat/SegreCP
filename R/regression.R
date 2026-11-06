#' Tensor regression with a low-CP-rank coefficient by RGD or RGN
#'
#' Fits scalar-on-tensor or tensor-on-tensor regression
#' \deqn{Y_m = \langle X_m, B \rangle + E_m, \qquad m = 1, \ldots, n,}{Y_m = <X_m, B> + E_m, m = 1, ..., n,}
#' where \eqn{X_m} has dimensions `input_dims`, the response \eqn{Y_m} is a
#' scalar or has dimensions `output_dims`, the coefficient \eqn{B} has
#' dimensions `c(input_dims, output_dims)` and CP rank `rank`, and
#' \eqn{\langle X_m, B \rangle}{<X_m, B>} contracts the input modes. The fit
#' minimizes
#' \eqn{\frac{1}{2n} \sum_m \|Y_m - \langle X_m, B \rangle\|_F^2}{(1/(2n)) sum_m ||Y_m - <X_m, B>||_F^2}
#' (or the same sum without the factor \eqn{1/n} if
#' `normalize_regression = "none"` in [segre_control()]) by Riemannian
#' gradient descent (RGD) or Riemannian Gauss--Newton (RGN) over the rank-one
#' components of \eqn{B}. One iteration updates each component once, in turn,
#' with the others at their current values. RGN solves the least-squares
#' problem restricted to the tangent space at the component. Every update is
#' retracted to a rank-one tensor by the rank-(1, ..., 1) T-HOSVD.
#'
#' The CPCA and HOSVD starts are computed from the moment
#' \eqn{\frac{1}{n \sigma_X^2} \sum_m X_m \otimes Y_m}{sum_m (X_m o Y_m) / (n sigma_X^2)}
#' (without the factor \eqn{1/n} if `normalize_regression = "none"`), which
#' estimates \eqn{B} when the design entries are uncorrelated with variance
#' \eqn{\sigma_X^2} (`design_variance`).
#'
#' @param X Design. The most efficient format is an `n` by `prod(input_dims)`
#'   matrix whose row `m` is `as.vector(X_m)`. Also accepted: an array with
#'   dimensions `c(n, input_dims)`, an `rTensor::Tensor` with these dimensions,
#'   or a list of `n` arrays.
#' @param Y Response: a numeric vector of length `n` (scalar response), an `n`
#'   by `prod(output_dims)` matrix, an array with dimensions
#'   `c(n, output_dims)`, an `rTensor::Tensor`, or a list of `n` arrays.
#' @param rank CP rank of the coefficient.
#' @param input_dims Dimensions of `X_m`. Needed only when `X` is a matrix and
#'   `X_m` has more than one mode.
#' @param output_dims Dimensions of `Y_m`. Needed only when `Y` is a matrix
#'   and `Y_m` has more than one mode.
#' @param method `"rgd"` or `"rgn"`.
#' @param init Start: `"cpca"`, `"hosvd"` or `"random"` (see
#'   [cp_initialize()]), a function called as `init(moment, rank, ...)`, or a
#'   list of CP components of the coefficient (see [cp_components()]).
#' @param init_split For `init = "cpca"`, the coefficient modes that index the
#'   rows of the unfolding (input modes are numbered first). By default the
#'   split that makes the unfolding closest to square.
#' @param init_args A named list of further arguments for [cp_initialize()] or
#'   for the function `init`.
#' @param seed Seed for `init = "random"`. The random number stream of the
#'   session is left unchanged.
#' @param design_variance Variance \eqn{\sigma_X^2} of the design entries. It
#'   scales the start and sets the default RGD step \eqn{1/\sigma_X^2}.
#' @param control Options from [segre_control()].
#' @param keep_history If `TRUE`, the components after every iteration are
#'   kept in `component_history`, for example to compute estimation errors with
#'   [convergence_data()].
#' @return An object of class `c("segre_cp_regression", "segre_cp_fit")`: a
#'   list with elements including
#'   * `components`: the CP components of the estimated coefficient, a list of
#'     `list(lambda, U)` with positive weights `lambda` and unit factors `U`
#'     (input modes first, then output modes);
#'   * `objective`: the objective at the start and after each iteration;
#'   * `iterations`, `converged`, `stopping_reason`: the number of iterations
#'     run and why the fit stopped;
#'   * `component_history`: if `keep_history = TRUE`, the components at the
#'     start and after each iteration;
#'   * `accepted_steps`: the RGD step of each component update (`NA` for RGN);
#'   * `normal_condition`: an estimate of the condition number of each RGN
#'     normal equation (`NA` for RGD);
#'   * `fitted_values`, `residual_values`: `n` by `prod(output_dims)`
#'     matrices on the scale of the data;
#'   * `initialization`, `control`, `call`: the start and the options used.
#'
#'   [coef()], [predict()], [fitted()], [residuals()], [summary()] and
#'   [convergence_data()] apply to the result. `cp_regression()` is an alias
#'   of `segre_regress()`.
#' @seealso [segre_decompose()] for CP decomposition and [segre_control()] for
#'   the options.
#' @examples
#' # Scalar response, 3 x 4 matrix covariates, rank-one coefficient
#' set.seed(2)
#' n <- 80
#' X <- matrix(rnorm(n * 12), n)
#' beta <- 2 * outer(c(1, 0, 0), c(0, 1, 0, 0))
#' y <- as.numeric(X %*% as.vector(beta)) + rnorm(n, sd = 0.05)
#' fit <- segre_regress(X, y, rank = 1, input_dims = c(3, 4), method = "rgn",
#'                      control = segre_control(max_iter = 5))
#' fit
#' predict(fit, X[1:3, ])
#'
#' # Matrix response: one output mode of dimension 2
#' Y <- cbind(y, -y) + rnorm(2 * n, sd = 0.05)
#' fit2 <- segre_regress(X, Y, rank = 1, input_dims = c(3, 4),
#'                       control = segre_control(max_iter = 5))
#' dim(coef(fit2))
#' @export
segre_regress <- function(
    X, Y, rank,
    input_dims = NULL,
    output_dims = NULL,
    method = c("rgd", "rgn"),
    init = c("cpca", "hosvd", "random"),
    init_split = NULL,
    init_args = list(),
    seed = NULL,
    design_variance = 1,
    control = segre_control(),
    keep_history = FALSE) {
  call <- match.call()
  method <- match.arg(method)
  if (length(design_variance) != 1L || !is.finite(design_variance) ||
      design_variance <= 0) {
    stop("`design_variance` must be one positive finite number.",
         call. = FALSE)
  }
  control <- .resolve_control(control, "regression", design_variance)
  rank <- .as_whole_number(rank, "rank")
  keep_history <- .as_flag(keep_history, "keep_history")

  design <- .as_design_matrix(X, input_dims)
  response <- .as_response_block(Y, design$n, output_dims)
  X_original <- design$X
  Y_original <- response$Y
  input_dims <- design$input_dims
  output_dims <- response$output_dims
  input_order <- length(input_dims)
  output_order <- length(output_dims)
  public_dims <- c(input_dims, output_dims)
  internal_dims <- c(output_dims, input_dims)
  if (length(public_dims) < 2L) {
    stop("The coefficient tensor must have at least two modes.", call. = FALSE)
  }

  normalization <- if (control$normalize_regression == "sqrt_n") {
    1 / sqrt(design$n)
  } else 1
  X_work <- if (normalization == 1 && is.double(X_original)) X_original else
    X_original * normalization
  Y_work <- if (normalization == 1 && is.double(Y_original)) Y_original else
    Y_original * normalization

  # A supplied component list needs only the coefficient dimensions, so the
  # moment X'Y is formed only for the other starts.
  if (is.list(init)) {
    initialization <- .component_start(init, public_dims, rank, init_args)
  } else {
    moment <- .coefficient_moment(
      X_work, Y_work, input_dims, output_dims, design_variance
    )
    initialization <- .resolve_initialization(
      moment, rank, init, init_split, seed, init_args
    )
  }
  public_components <- .validate_components(
    initialization$comps, public_dims, rank
  )
  component_history <- if (keep_history) list(public_components) else NULL
  components <- .public_to_internal_components(
    public_components, input_order, output_order
  )
  components <- .validate_components(components, internal_dims, rank)
  split <- .split_internal_dims(internal_dims, output_order)

  predictions <- lapply(components, function(comp) {
    .apply_rank1_operator(comp, X_work, split)$prediction
  })
  prediction <- Reduce(`+`, predictions)
  residual <- prediction - Y_work
  objective <- numeric(control$max_iter + 1L)
  objective[[1L]] <- 0.5 * sum(residual * residual)
  accepted_steps <- matrix(NA_real_, control$max_iter, rank)
  normal_condition <- matrix(NA_real_, control$max_iter, rank)
  converged <- FALSE
  stopping_reason <- "max_iter"
  iterations <- control$max_iter

  for (iteration in seq_len(control$max_iter)) {
    order_i <- .component_update_order(components, control)
    for (i in order_i) {
      old_component <- components[[i]]
      old_prediction <- predictions[[i]]

      if (method == "rgd") {
        need_gram <- control$step_rule == "exact_quadratic"
        gradient <- .regression_gradient_stats(
          X_work, residual, old_component, internal_dims, output_order,
          need_gram = need_gram
        )
        base <- c(old_component$lambda,
                  numeric(length(gradient$rhs) - 1L))
        context <- .prepare_tangent_retraction(
          old_component$U, base, -gradient$rhs
        )
        alpha <- .step_at(control, iteration)
        if (control$step_rule == "exact_quadratic") {
          numerator <- sum(gradient$rhs * gradient$rhs)
          denominator <- as.numeric(crossprod(
            gradient$rhs, gradient$gram %*% gradient$rhs
          ))
          if (is.finite(denominator) && denominator > .Machine$double.eps) {
            alpha <- numerator / denominator
          }
          alpha <- min(max(alpha, 1e-10), 1e2)
        }

        if (control$step_rule == "backtracking") {
          loss_old <- 0.5 * sum(residual * residual)
          gradient_norm_sq <- sum(gradient$rhs * gradient$rhs)
          partial <- residual - old_prediction
          best <- list(component = old_component,
                       prediction = old_prediction,
                       residual = residual,
                       loss = loss_old, step = 0)
          for (attempt in seq_len(control$max_backtracks)) {
            candidate_component <- .apply_prepared_retraction(context, alpha)
            candidate_prediction <- .apply_rank1_operator(
              candidate_component, X_work, split
            )$prediction
            candidate_residual <- partial + candidate_prediction
            loss_new <- 0.5 * sum(candidate_residual * candidate_residual)
            if (loss_new < best$loss) {
              best <- list(component = candidate_component,
                           prediction = candidate_prediction,
                           residual = candidate_residual,
                           loss = loss_new, step = alpha)
            }
            if (loss_new <= loss_old - control$armijo * alpha *
                gradient_norm_sq) break
            alpha <- alpha * control$backtrack_factor
          }
          new_component <- best$component
          new_prediction <- best$prediction
          residual <- best$residual
          alpha <- best$step
        } else {
          new_component <- .apply_prepared_retraction(context, alpha)
          new_prediction <- .apply_rank1_operator(
            new_component, X_work, split
          )$prediction
          residual <- residual + new_prediction - old_prediction
        }
        accepted_steps[iteration, i] <- alpha
      } else {
        # RGN: least squares over the tangent space at T_i with the other
        # components fixed. Its solution is T_i - Delta_i, where Delta_i solves
        # the RGN normal equation of the paper.
        target <- Y_work - (prediction - old_prediction)
        if (ncol(Y_work) == 1L) {
          tangent_design <- .compute_tangent_design(X_work, old_component$U)
          solved <- .solve_scalar_tangent(
            tangent_design, target, control$ridge
          )
          coordinates <- solved$solution
          normal_condition[iteration, i] <- solved$condition
        } else {
          statistics <- .tensor_tangent_statistics(
            X_work, target, old_component$U, internal_dims, output_order,
            include_gram = TRUE, full_gram = FALSE
          )
          coordinates <- .solve_tensor_tangent(statistics, control$ridge)
          joint_gram <- statistics$gram_in
          if (control$ridge > 0) joint_gram <- joint_gram +
            control$ridge * diag(statistics$df_in)
          normal_condition[iteration, i] <- kappa(joint_gram, exact = FALSE)
        }
        new_component <- .retract_tangent_coordinates(
          old_component$U, coordinates
        )
        new_prediction <- .apply_rank1_operator(
          new_component, X_work, split
        )$prediction
        residual <- residual + new_prediction - old_prediction
      }

      if (!is.finite(new_component$lambda) ||
          new_component$lambda <= .Machine$double.eps) {
        stop(sprintf("Component %d retracted to a zero/nonfinite tensor.", i),
             call. = FALSE)
      }
      components[[i]] <- new_component
      predictions[[i]] <- new_prediction
      prediction <- prediction + new_prediction - old_prediction
    }

    objective[[iteration + 1L]] <- 0.5 * sum(residual * residual)
    if (keep_history) {
      component_history[[iteration + 1L]] <-
        .internal_to_public_components(components, input_order, output_order)
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

  public_components <- .internal_to_public_components(
    components, input_order, output_order
  )
  public_components <- .validate_components(
    public_components, public_dims, rank
  )
  # Fitted values on the scale of the data. Without normalization the stored
  # component predictions are reused.
  fitted_values <- if (normalization == 1) Reduce(`+`, predictions) else
    Reduce(`+`, lapply(components, function(comp) {
      .apply_rank1_operator(comp, X_original, split)$prediction
    }))
  residual_values <- Y_original - fitted_values

  structure(list(
    call = call,
    task = "regression",
    method = method,
    profile = control$profile,
    rank = rank,
    dims = public_dims,
    input_dims = input_dims,
    output_dims = output_dims,
    components = public_components,
    initialization = initialization[names(initialization) != "comps"],
    objective = objective[seq_len(iterations + 1L)],
    component_history = component_history,
    accepted_steps = accepted_steps[seq_len(iterations), , drop = FALSE],
    normal_condition = normal_condition[seq_len(iterations), , drop = FALSE],
    iterations = iterations,
    converged = converged,
    stopping_reason = stopping_reason,
    control = control,
    normalization = normalization,
    design_variance = design_variance,
    fitted_values = fitted_values,
    residual_values = residual_values,
    response_source = response$source,
    design_source = design$source
  ), class = c("segre_cp_regression", "segre_cp_fit"))
}

#' @rdname segre_regress
#' @export
cp_regression <- segre_regress

# Internal helpers ---------------------------------------------------------

# Public order: input modes, then output modes. Internal order: output modes,
# then input modes.
.public_to_internal_components <- function(components, input_order, output_order) {
  if (!output_order) return(components)
  permutation <- c(input_order + seq_len(output_order), seq_len(input_order))
  lapply(components, function(comp) {
    list(lambda = comp$lambda, U = comp$U[permutation])
  })
}

.internal_to_public_components <- function(components, input_order, output_order) {
  if (!output_order) return(components)
  permutation <- c(output_order + seq_len(input_order), seq_len(output_order))
  lapply(components, function(comp) {
    list(lambda = comp$lambda, U = comp$U[permutation])
  })
}

.regression_gradient_stats <- function(X, residual, component, dims,
                                       output_order, need_gram = FALSE) {
  if (ncol(residual) == 1L) {
    design <- .compute_tangent_design(X, component$U)
    ans <- list(rhs = as.numeric(crossprod(design, residual)))
    if (need_gram) ans$gram <- crossprod(design)
    return(ans)
  }
  .tensor_tangent_statistics(
    X, residual, component$U, dims, output_order,
    include_gram = need_gram
  )
}

.solve_scalar_tangent <- function(design, target, ridge) {
  gram <- crossprod(design)
  if (ridge > 0) gram <- gram + ridge * diag(ncol(design))
  rhs <- as.numeric(crossprod(design, target))
  solution <- tryCatch(
    solve(gram, rhs),
    error = function(e) stop(
      "The restricted tangent normal operator is singular; use a better ",
      "start or set a positive `ridge` in segre_control().",
      call. = FALSE
    )
  )
  list(solution = solution, condition = kappa(gram, exact = FALSE))
}

.coefficient_moment <- function(X, Y, input_dims, output_dims,
                                design_variance) {
  moment <- crossprod(X, Y) / design_variance
  dims <- c(input_dims, output_dims)
  array(as.numeric(moment), dim = dims)
}
