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
      "initialization or set a positive `ridge` in the practical profile.",
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

#' Fit CP tensor regression with RGD or RGN
#'
#' Fits scalar-on-tensor or tensor-on-tensor least squares with a CP-structured
#' coefficient. Public coefficient modes always follow the manuscript order:
#' **input modes first, output modes second**. Internally the implementation
#' applies a fixed output-first permutation to exploit the exact separated
#' tangent design; the returned factors are permuted back.
#'
#' `profile = "manuscript"` (the default in [segre_control()]) uses normalized
#' observations, fixed component order, the fixed RGD step
#' \eqn{\alpha = \sigma_X^{-2}} unless `step_size` is given, exact T-HOSVD, no
#' early stopping, and unregularized RGN. The practical profile records every
#' deviation explicitly in the fit object.
#'
#' @param X Tensor covariates. The efficient format is an `n` by `p_in` matrix.
#'   Also accepted: an array with dimensions `c(n, input_dims)`, an
#'   `rTensor::Tensor` with those dimensions, or a length-`n` list of tensors.
#' @param Y A length-`n` vector for scalar response; an `n` by `q` matrix; an
#'   array with dimensions `c(n, output_dims)`; an `rTensor::Tensor`; or a list
#'   of response tensors.
#' @param rank Positive target CP rank.
#' @param input_dims Input tensor mode dimensions. Required when a matrix `X`
#'   should represent more than one input mode.
#' @param output_dims Output tensor mode dimensions. For a matrix `Y`, omission
#'   treats its columns as one output mode. A vector `Y` is scalar.
#' @param method `"rgd"` or `"rgn"`.
#' @param init Initialization accepted by [cp_initialize()]. CPCA and HOSVD are
#'   both applied to the normalized adjoint moment in public input-output order.
#' @param init_split Optional CPCA mode split in public coefficient-mode order.
#' @param init_args Named list passed to [cp_initialize()].
#' @param seed Optional local seed for random initialization.
#' @param design_variance Positive design-entry variance \eqn{\sigma_X^2}. It
#'   scales the initialization moment: with `normalize_regression = "sqrt_n"`,
#'   the moment is `crossprod(X, Y) / (n * design_variance)` on the original data
#'   scale. The default RGD step is `1 / design_variance`.
#' @param control Controls from [segre_control()].
#' @param keep_history If `TRUE`, retain the public-order CP components at
#'   initialization and after every completed sweep. This is intended for
#'   convergence diagnostics and simulation reproducibility.
#' @return A `segre_cp_regression` fit. `components` and [coef()] use input modes
#'   followed by output modes.
#' @examples
#' set.seed(2)
#' n <- 80
#' input_dims <- c(3, 4)
#' X <- matrix(rnorm(n * prod(input_dims)), n)
#' u1 <- c(1, 0, 0)
#' u2 <- c(1, 0, 0, 0)
#' beta <- 2 * kronecker(u2, u1)
#' y <- as.numeric(X %*% beta + rnorm(n, sd = 0.05))
#' fit <- segre_regress(
#'   X, y, rank = 1, input_dims = input_dims, method = "rgn",
#'   init = "hosvd", control = segre_control(max_iter = 3)
#' )
#' predict(fit, X[1:3, , drop = FALSE])
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

  if (is.list(init)) {
    if (!is.list(init_args)) stop("`init_args` must be a list.", call. = FALSE)
    # Supplied components use only the coefficient dimensions, not X'Y.
    # Retain cp_initialize's argument matching and first validation pass.
    initialize_components <- function(x, rank, method, split = NULL,
                                      seed = NULL,
                                      project_random_weights = TRUE, ...) {
      list(comps = .validate_components(method, public_dims, rank),
           method = "custom_components")
    }
    initialization <- do.call(initialize_components, c(
      list(x = NULL, rank = rank, method = init, split = init_split,
           seed = seed), init_args
    ))
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
        # The tangent-space least-squares target equals T_i - Delta_i, where
        # Delta_i solves the manuscript RGN normal equation.
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

  public_components <- .internal_to_public_components(
    components, input_order, output_order
  )
  public_components <- .validate_components(
    public_components, public_dims, rank
  )
  # Re-sum component predictions, rather than the incrementally updated sum,
  # to preserve the original final summation order without rescanning X.
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
