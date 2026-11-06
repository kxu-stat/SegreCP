# Tensor-on-tensor regression methods of Section 5 and Appendix C: CP-ALS,
# RGD and RGN (SegreCP) from a shared start, RRR (MultiwayRegression),
# Tucker-RGD, TT-IHT and HOPLS.
#
# X is an n x prod(input_dims) matrix with rows vec(X_m), and Y has dimensions
# c(n, output_dims). Coefficients have the input modes first, then the output
# modes. Internally the output modes come first.
#
# A fit returns a list with
#   coefficient          the final coefficient array;
#   coefficient_history  the coefficient at the start and after each iteration
#                        (for RGD and RGN: component_history, the CP components);
#   objective            half the residual sum of squares at the start and
#                        after each iteration (CP-ALS, RGD, RGN, Tucker-RGD,
#                        TT-IHT).
# The paths are recorded only with record_path = TRUE (not in the timing
# studies).

.regression_problem <- function(X, Y, input_dims, output_dims) {
  input_dims <- as.integer(input_dims)
  output_dims <- as.integer(output_dims)
  n <- nrow(X)
  q <- if (length(output_dims)) prod(output_dims) else 1L
  list(
    X = unname(X), Y = matrix(as.numeric(Y), nrow = n, ncol = q), n = n,
    input_dims = input_dims, output_dims = output_dims,
    public_dims = c(input_dims, output_dims),
    internal_dims = c(output_dims, input_dims),
    input_order = length(input_dims), output_order = length(output_dims),
    p_in = prod(input_dims), q_out = q
  )
}

.response_from_matrix <- function(Y, output_dims) {
  if (!length(output_dims)) return(as.numeric(Y[, 1L]))
  array(as.numeric(Y), dim = c(nrow(Y), output_dims))
}

.predict_coefficient_matrix <- function(coefficient, X, input_dims, output_dims) {
  q <- if (length(output_dims)) prod(output_dims) else 1L
  B <- matrix(as.numeric(coefficient), nrow = prod(input_dims), ncol = q)
  unname(X %*% B)
}

.kron_vector <- function(factors) {
  if (!length(factors)) return(1)
  if (length(factors) == 1L) return(as.numeric(factors[[1L]]))
  Reduce(kronecker, rev(lapply(factors, as.numeric)))
}

.components_to_array <- function(components, dims) {
  value <- numeric(prod(dims))
  for (component in components) {
    value <- value + component$lambda * .kron_vector(component$U)
  }
  array(value, dim = dims)
}

.public_to_internal_components <- function(components, problem) {
  if (!problem$output_order) return(components)
  permutation <- c(problem$input_order + seq_len(problem$output_order),
                   seq_len(problem$input_order))
  lapply(components, function(component) {
    list(lambda = as.numeric(component$lambda), U = component$U[permutation])
  })
}

.internal_to_public_components <- function(components, problem) {
  if (!problem$output_order) return(components)
  permutation <- c(problem$output_order + seq_len(problem$input_order),
                   seq_len(problem$output_order))
  lapply(components, function(component) {
    list(lambda = as.numeric(component$lambda), U = component$U[permutation])
  })
}

.internal_to_public_array <- function(A, problem) {
  if (!problem$output_order) return(array(A, dim = problem$public_dims))
  permutation <- c(problem$output_order + seq_len(problem$input_order),
                   seq_len(problem$output_order))
  aperm(array(A, dim = problem$internal_dims), permutation)
}

.apply_internal_operator <- function(A, X, problem) {
  A_matrix <- matrix(as.numeric(A), nrow = problem$q_out,
                     ncol = problem$p_in)
  unname(X %*% t(A_matrix))
}

.internal_adjoint <- function(Z, X, problem) {
  array(as.numeric(crossprod(Z, X)), dim = problem$internal_dims)
}

.mode_unfold <- function(A, mode) {
  dimensions <- dim(A)
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  matrix(as.numeric(aperm(A, permutation)), nrow = dimensions[[mode]])
}

.mode_fold <- function(M, mode, dimensions) {
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  aperm(array(as.numeric(M), dim = dimensions[permutation]), order(permutation))
}

.mode_product <- function(A, M, mode) {
  dimensions <- dim(A)
  output_dimensions <- dimensions
  output_dimensions[[mode]] <- nrow(M)
  .mode_fold(M %*% .mode_unfold(A, mode), mode, output_dimensions)
}

.objective_path_from_components <- function(history, problem) {
  vapply(history, function(components) {
    coefficient <- .components_to_array(components, problem$public_dims)
    residual <- .predict_coefficient_matrix(
      coefficient, problem$X, problem$input_dims, problem$output_dims
    ) - problem$Y
    0.5 * sum(residual^2)
  }, numeric(1L))
}

# Shared start of CP-ALS, RGD and RGN: cp_initialize() applied to the moment
# X'Y / design_variance, with method "cpca", "hosvd" or "random".
regression_start <- function(X, Y, input_dims, output_dims, rank, split = NULL,
                             design_variance = 1, method = "cpca",
                             seed = NULL) {
  problem <- .regression_problem(X, Y, input_dims, output_dims)
  moment <- crossprod(problem$X, problem$Y) / design_variance
  moment <- array(as.numeric(moment), dim = problem$public_dims)
  SegreCP::cp_initialize(moment, rank, method = method, split = split,
                         seed = seed)$comps
}

# RGD and RGN: SegreCP with the components updated in decreasing order of
# their weights, a fixed RGD step and no early stopping.
.fit_segrecp <- function(method, problem, rank, start, iterations, step_size,
                         record_path) {
  control <- SegreCP::segre_control(
    profile = "manuscript", max_iter = iterations, tolerance = 0,
    step_size = step_size, step_rule = "fixed",
    component_order = "decreasing_weight", normalize_regression = "none"
  )
  fit <- SegreCP::segre_regress(
    X = problem$X, Y = .response_from_matrix(problem$Y, problem$output_dims),
    rank = rank, input_dims = problem$input_dims,
    output_dims = if (length(problem$output_dims)) problem$output_dims else NULL,
    method = tolower(method), init = start, control = control,
    keep_history = record_path
  )
  list(
    coefficient = stats::coef(fit), components = fit$components,
    component_history = fit$component_history,
    objective = if (record_path)
      .objective_path_from_components(fit$component_history, problem)
  )
}

# CP-ALS: alternating least squares over the factors of each mode.
.fit_cp_als <- function(problem, rank, start, iterations, record_path) {
  components <- .public_to_internal_components(start, problem)
  dimensions <- problem$internal_dims
  output_order <- problem$output_order
  output_modes <- if (output_order) seq_len(output_order) else integer()
  input_modes <- output_order + seq_len(problem$input_order)
  factors <- lapply(seq_along(dimensions), function(mode) {
    do.call(cbind, lapply(components, function(component) component$U[[mode]]))
  })
  weights <- vapply(components, `[[`, numeric(1L), "lambda")
  current_coefficient <- function() {
    current <- lapply(seq_len(rank), function(component) {
      list(lambda = weights[[component]], U = lapply(factors, function(A) A[, component]))
    })
    .components_to_array(.internal_to_public_components(current, problem),
                         problem$public_dims)
  }
  prediction <- function() {
    out <- matrix(0, problem$n, problem$q_out)
    for (component in seq_len(rank)) {
      input_vector <- .kron_vector(lapply(input_modes, function(mode) {
        factors[[mode]][, component]
      }))
      score <- as.numeric(problem$X %*% input_vector)
      output_vector <- if (output_order) {
        .kron_vector(lapply(output_modes, function(mode) {
          factors[[mode]][, component]
        }))
      } else 1
      out <- out + weights[[component]] * tcrossprod(score, output_vector)
    }
    out
  }
  coefficient_history <- NULL
  objective <- NULL
  if (record_path) {
    coefficient_history <- list(current_coefficient())
    objective <- 0.5 * sum((problem$Y - prediction())^2)
  }

  for (iteration in seq_len(iterations)) {
    for (mode in seq_along(dimensions)) {
      mode_size <- dimensions[[mode]]
      if (mode %in% input_modes) {
        output_gram <- matrix(1, rank, rank)
        for (output_mode in output_modes) {
          output_gram <- output_gram * crossprod(factors[[output_mode]])
        }
        input_position <- match(mode, input_modes)
        fibre_stack <- matrix(0, problem$n, mode_size * rank)
        for (component in seq_len(rank)) {
          input_factors <- lapply(input_modes, function(index) {
            factors[[index]][, component]
          })
          fibres <- SegreCP:::modewise_design_fibres_cpp(problem$X, input_factors)
          columns <- ((component - 1L) * mode_size + 1L):
            (component * mode_size)
          fibre_stack[, columns] <- fibres[[input_position]]
        }
        cross_gram <- crossprod(fibre_stack)
        scale <- (weights %o% weights) * output_gram
        normal <- cross_gram * kronecker(scale, matrix(1, mode_size, mode_size))
        output_matrix <- matrix(1, 1L, rank)
        if (output_order) {
          output_matrix <- vapply(seq_len(rank), function(component) {
            .kron_vector(lapply(output_modes, function(index) {
              factors[[index]][, component]
            }))
          }, numeric(problem$q_out))
          output_matrix <- matrix(output_matrix, nrow = problem$q_out,
                                  ncol = rank)
        }
        response_projection <- problem$Y %*% output_matrix
        rhs_matrix <- matrix(0, mode_size, rank)
        for (component in seq_len(rank)) {
          columns <- ((component - 1L) * mode_size + 1L):
            (component * mode_size)
          rhs_matrix[, component] <- weights[[component]] *
            crossprod(fibre_stack[, columns, drop = FALSE],
                      response_projection[, component])
        }
        updated <- matrix(solve(normal, as.numeric(rhs_matrix)),
                          nrow = mode_size, ncol = rank)
      } else {
        output_gram <- matrix(1, rank, rank)
        for (other in setdiff(output_modes, mode)) {
          output_gram <- output_gram * crossprod(factors[[other]])
        }
        input_kronecker <- vapply(seq_len(rank), function(component) {
          .kron_vector(lapply(input_modes, function(index) {
            factors[[index]][, component]
          }))
        }, numeric(problem$p_in))
        input_kronecker <- matrix(input_kronecker,
                                  nrow = problem$p_in, ncol = rank)
        scores <- problem$X %*% input_kronecker
        scale <- (weights %o% weights) * crossprod(scores) * output_gram
        rhs_matrix <- matrix(0, mode_size, rank)
        for (component in seq_len(rank)) {
          response_vector <- as.numeric(crossprod(
            problem$Y, scores[, component]
          ))
          if (output_order == 1L) {
            contracted <- response_vector
          } else {
            response_array <- array(response_vector,
                                    dim = dimensions[output_modes])
            for (other in setdiff(output_modes, mode)) {
              response_array <- .mode_product(
                response_array,
                matrix(factors[[other]][, component], nrow = 1L),
                other
              )
            }
            contracted <- as.numeric(response_array)
          }
          rhs_matrix[, component] <- weights[[component]] * contracted
        }
        updated <- t(solve(scale, t(rhs_matrix)))
      }
      norms <- sqrt(colSums(updated^2))
      factors[[mode]] <- sweep(updated, 2L, norms, "/")
      weights <- weights * norms
    }
    if (record_path) {
      objective[[iteration + 1L]] <- 0.5 * sum((problem$Y - prediction())^2)
      coefficient_history[[iteration + 1L]] <- current_coefficient()
    }
  }
  components <- lapply(seq_len(rank), function(component) {
    list(lambda = weights[[component]], U = lapply(factors, function(A) {
      A[, component]
    }))
  })
  public_components <- .internal_to_public_components(components, problem)
  list(
    coefficient = .components_to_array(public_components, problem$public_dims),
    components = public_components,
    coefficient_history = coefficient_history, objective = objective
  )
}

.tucker_hosvd <- function(A, ranks) {
  factors <- lapply(seq_along(dim(A)), function(mode) {
    current_rank <- min(ranks[[mode]], dim(A)[[mode]])
    svd(.mode_unfold(A, mode), nu = current_rank, nv = 0L)$u[
      , seq_len(current_rank), drop = FALSE
    ]
  })
  core <- A
  for (mode in seq_along(factors)) {
    core <- .mode_product(core, t(factors[[mode]]), mode)
  }
  reconstruction <- core
  for (mode in seq_along(factors)) {
    reconstruction <- .mode_product(reconstruction, factors[[mode]], mode)
  }
  list(core = core, factors = factors, recon = reconstruction)
}

.tucker_reconstruct <- function(core, factors) {
  value <- core
  for (mode in seq_along(factors)) {
    value <- .mode_product(value, factors[[mode]], mode)
  }
  value
}

# Projection of Z onto the tangent space of the Tucker manifold at the
# current core and factors.
.tucker_tangent_projection <- function(Z, state) {
  core <- state$core
  factors <- state$factors
  order <- length(factors)
  dot_core <- Z
  for (mode in seq_len(order)) {
    dot_core <- .mode_product(dot_core, t(factors[[mode]]), mode)
  }
  projected <- .tucker_reconstruct(dot_core, factors)
  for (mode in seq_len(order)) {
    contracted <- Z
    for (other in setdiff(seq_len(order), mode)) {
      contracted <- .mode_product(
        contracted, t(factors[[other]]), other
      )
    }
    Z_mode <- .mode_unfold(contracted, mode)
    G_mode <- .mode_unfold(core, mode)
    gram <- G_mode %*% t(G_mode)
    factor <- factors[[mode]]
    dot_factor <- Z_mode %*% t(G_mode) %*% solve(gram)
    dot_factor <- dot_factor - factor %*% crossprod(factor, dot_factor)
    tangent_factors <- factors
    tangent_factors[[mode]] <- dot_factor
    projected <- projected + .tucker_reconstruct(core, tangent_factors)
  }
  projected
}

# Tucker-RGD: projected gradient step with exact line search, halved up to 20
# times until the objective does not increase, and T-HOSVD retraction. The
# start is the T-HOSVD of X'Y (simulations) or of the least-squares estimate
# (trade data).
.fit_tucker_rgd <- function(problem, rank, iterations, real_data, record_path) {
  ranks <- rep(rank, length(problem$internal_dims))
  if (real_data) {
    beta <- solve(crossprod(problem$X), crossprod(problem$X, problem$Y))
    start <- array(as.numeric(t(beta)), dim = problem$internal_dims)
  } else {
    start <- .internal_adjoint(problem$Y, problem$X, problem)
  }
  state <- .tucker_hosvd(start, ranks)
  coefficient_internal <- state$recon
  loss <- function(A) {
    0.5 * sum((.apply_internal_operator(A, problem$X, problem) - problem$Y)^2)
  }
  objective <- numeric(iterations + 1L)
  objective[[1L]] <- loss(coefficient_internal)
  coefficient_history <- if (record_path)
    list(.internal_to_public_array(coefficient_internal, problem))
  for (iteration in seq_len(iterations)) {
    residual <- .apply_internal_operator(
      coefficient_internal, problem$X, problem
    ) - problem$Y
    ambient <- .internal_adjoint(residual, problem$X, problem)
    gradient <- .tucker_tangent_projection(ambient, state)
    measured <- .apply_internal_operator(gradient, problem$X, problem)
    step <- sum(gradient^2) / max(sum(measured^2), .Machine$double.eps)
    candidate_loss <- objective[[iteration]]
    for (backtrack in seq_len(20L)) {
      trial_state <- .tucker_hosvd(
        coefficient_internal - step * gradient, ranks
      )
      trial_loss <- loss(trial_state$recon)
      if (is.finite(trial_loss) && trial_loss <= objective[[iteration]]) {
        coefficient_internal <- trial_state$recon
        state <- trial_state
        candidate_loss <- trial_loss
        break
      }
      step <- step / 2
    }
    objective[[iteration + 1L]] <- candidate_loss
    if (record_path) coefficient_history[[iteration + 1L]] <-
      .internal_to_public_array(coefficient_internal, problem)
  }
  list(
    coefficient = .internal_to_public_array(coefficient_internal, problem),
    coefficient_history = coefficient_history,
    objective = if (record_path) objective
  )
}

# TT-SVD projection of A onto tensors of TT rank tt_ranks.
.tt_svd_project <- function(A, tt_ranks) {
  dimensions <- dim(A)
  order <- length(dimensions)
  cores <- vector("list", order)
  remainder <- A
  left_rank <- 1L
  for (mode in seq_len(order - 1L)) {
    unfolding <- matrix(remainder, nrow = left_rank * dimensions[[mode]])
    current_rank <- min(tt_ranks[[mode]], nrow(unfolding), ncol(unfolding))
    decomposition <- svd(unfolding, nu = current_rank, nv = current_rank)
    cores[[mode]] <- array(
      decomposition$u[, seq_len(current_rank), drop = FALSE],
      dim = c(left_rank, dimensions[[mode]], current_rank)
    )
    remainder <- array(
      diag(decomposition$d[seq_len(current_rank)], nrow = current_rank) %*%
        t(decomposition$v[, seq_len(current_rank), drop = FALSE]),
      dim = c(current_rank, dimensions[(mode + 1L):order])
    )
    left_rank <- current_rank
  }
  cores[[order]] <- array(remainder,
                          dim = c(left_rank, dimensions[[order]], 1L))
  current <- array(cores[[1L]], dim = c(dimensions[[1L]],
                                       dim(cores[[1L]])[[3L]]))
  for (mode in 2:order) {
    left_rank <- dim(cores[[mode]])[[1L]]
    right_rank <- dim(cores[[mode]])[[3L]]
    current <- array(
      matrix(current, nrow = prod(dimensions[seq_len(mode - 1L)]),
             ncol = left_rank) %*%
        matrix(cores[[mode]], nrow = left_rank,
               ncol = dimensions[[mode]] * right_rank),
      dim = c(dimensions[seq_len(mode)], right_rank)
    )
  }
  array(current, dim = dimensions)
}

# TT-IHT: gradient step followed by TT-SVD, with step 1 (simulations, where
# sigma_X = 1) or 1 / lambda_max(X'X) (trade data).
.fit_tt_iht <- function(problem, rank, iterations, real_data, record_path) {
  tt_ranks <- rep(rank, length(problem$internal_dims) - 1L)
  adjoint <- .internal_adjoint(problem$Y, problem$X, problem)
  step <- if (real_data) {
    1 / max(eigen(crossprod(problem$X), symmetric = TRUE,
                  only.values = TRUE)$values)
  } else 1
  coefficient_internal <- .tt_svd_project(step * adjoint, tt_ranks)
  loss <- function(A) {
    0.5 * sum((.apply_internal_operator(A, problem$X, problem) - problem$Y)^2)
  }
  coefficient_history <- NULL
  objective <- NULL
  if (record_path) {
    coefficient_history <- list(.internal_to_public_array(coefficient_internal, problem))
    objective <- loss(coefficient_internal)
  }
  for (iteration in seq_len(iterations)) {
    residual <- .apply_internal_operator(
      coefficient_internal, problem$X, problem
    ) - problem$Y
    gradient <- .internal_adjoint(residual, problem$X, problem)
    coefficient_internal <- .tt_svd_project(
      coefficient_internal - step * gradient, tt_ranks
    )
    if (record_path) {
      objective[[iteration + 1L]] <- loss(coefficient_internal)
      coefficient_history[[iteration + 1L]] <-
        .internal_to_public_array(coefficient_internal, problem)
    }
  }
  list(
    coefficient = .internal_to_public_array(coefficient_internal, problem),
    coefficient_history = coefficient_history, objective = objective
  )
}

# Leading rank-one approximation of A by HOOI from the HOSVD start.
.rank_one_hooi <- function(A, max_iter = 100L, tolerance = 1e-10) {
  dimensions <- dim(A)
  factors <- lapply(seq_along(dimensions), function(mode) {
    svd(.mode_unfold(A, mode), nu = 1L, nv = 0L)$u[, 1L]
  })
  contract_except <- function(mode) {
    value <- A
    for (other in setdiff(seq_along(dimensions), mode)) {
      value <- .mode_product(value, matrix(factors[[other]], nrow = 1L), other)
    }
    as.numeric(value)
  }
  previous <- 0
  for (iteration in seq_len(max_iter)) {
    for (mode in seq_along(dimensions)) {
      update <- contract_except(mode)
      update_norm <- sqrt(sum(update^2))
      if (update_norm > .Machine$double.eps) factors[[mode]] <- update / update_norm
    }
    value <- sum(as.numeric(A) * .kron_vector(factors))
    if (value < 0) {
      factors[[1L]] <- -factors[[1L]]
      value <- -value
    }
    if (abs(value - previous) <= tolerance * max(1, abs(previous))) break
    previous <- value
  }
  factors
}

# HOPLS with rank-one loadings: each latent block takes the leading rank-one
# approximation of the cross-covariance of the residual X and Y, and both are
# deflated.
.fit_hopls <- function(problem, blocks) {
  X_residual <- problem$X
  Y_residual <- problem$Y
  input_loadings <- vector("list", blocks)
  output_loadings <- vector("list", blocks)
  response_core <- numeric(blocks)
  effective_weights <- vector("list", blocks)
  output_modes <- if (problem$output_order) seq_len(problem$output_order) else integer()
  input_modes <- problem$output_order + seq_len(problem$input_order)

  for (component in seq_len(blocks)) {
    cross_covariance <- array(
      as.numeric(crossprod(Y_residual, X_residual)),
      dim = problem$internal_dims
    )
    factors <- .rank_one_hooi(cross_covariance)
    output_loading <- if (problem$output_order) {
      .kron_vector(factors[output_modes])
    } else 1
    input_loading <- .kron_vector(factors[input_modes])
    raw_score <- as.numeric(X_residual %*% input_loading)
    score_norm <- sqrt(sum(raw_score^2))
    score <- raw_score / score_norm
    response_core[[component]] <- sum(
      score * as.numeric(Y_residual %*% output_loading)
    )
    effective_loading <- input_loading
    if (component > 1L) {
      for (previous in rev(seq_len(component - 1L))) {
        direction <- input_loadings[[previous]]
        effective_loading <- effective_loading -
          direction * sum(direction * effective_loading)
      }
    }
    effective_loading <- effective_loading / score_norm
    input_loadings[[component]] <- input_loading
    output_loadings[[component]] <- output_loading
    effective_weights[[component]] <- effective_loading
    block_fit <- tcrossprod(
      score * response_core[[component]], output_loading
    )
    X_residual <- X_residual - tcrossprod(raw_score, input_loading)
    Y_residual <- Y_residual - block_fit
  }
  coefficient_matrix <- Reduce(`+`, lapply(seq_len(blocks), function(component) {
    tcrossprod(
      effective_weights[[component]],
      response_core[[component]] * output_loadings[[component]]
    )
  }))
  list(coefficient = array(as.numeric(coefficient_matrix),
                           dim = problem$public_dims))
}

# MultiwayRegression::rrr with an iteration limit on its outer loop and, if
# coefficient_callback is given, a call with the coefficient at the start and
# after each outer iteration.
.make_rrr_fixed <- function(max_iterations, coefficient_callback = NULL) {
  fit_function <- MultiwayRegression::rrr
  function_formals <- formals(fit_function)
  function_formals$maxIter <- max_iterations
  formals(fit_function) <- function_formals
  patch_loop <- function(expression) {
    if (is.call(expression) && identical(expression[[1L]], as.name("while")) &&
        identical(paste(deparse(expression[[2L]]), collapse = ""),
                  "conv > convThresh")) {
      old_condition <- expression[[2L]]
      expression[[2L]] <- substitute(
        OLD && j < maxIter, list(OLD = old_condition)
      )
      if (!is.null(coefficient_callback)) {
        outer_body <- expression[[3L]]
        expression[[3L]] <- substitute({
          ORIGINAL
          .rrr_coefficient_callback(j, array(t(Bmat) %*% t(Vmat), dim = c(P, Q)))
        }, list(ORIGINAL = outer_body))
        return(substitute({
          .rrr_coefficient_callback(0L, array(t(Bmat) %*% t(Vmat), dim = c(P, Q)))
          LOOP
        }, list(LOOP = expression)))
      }
    }
    if (is.recursive(expression)) {
      for (index in seq_along(expression)) {
        expression[[index]] <- patch_loop(expression[[index]])
      }
    }
    expression
  }
  body(fit_function) <- patch_loop(body(fit_function))
  if (!is.null(coefficient_callback)) {
    callback_environment <- new.env(parent = environment(fit_function))
    callback_environment$.rrr_coefficient_callback <- coefficient_callback
    environment(fit_function) <- callback_environment
  }
  fit_function
}

# RRR: random start with the given seed, no ridge, no annealing.
.fit_rrr <- function(problem, rank, iterations, seed, record_path) {
  coefficient_history <- if (record_path) list()
  callback <- if (record_path) function(iteration, coefficient) {
    coefficient_history[[as.integer(iteration) + 1L]] <<-
      array(as.numeric(coefficient), dim = problem$public_dims)
  }
  fit_function <- .make_rrr_fixed(iterations, coefficient_callback = callback)
  # A single input mode is passed to rrr() as a p x 1 array.
  rrr_input_dims <- if (length(problem$input_dims) == 1L) {
    c(problem$input_dims, 1L)
  } else problem$input_dims
  X_array <- array(as.numeric(problem$X), dim = c(problem$n, rrr_input_dims))
  rrr_output_dims <- if (length(problem$output_dims)) problem$output_dims else 1L
  Y_array <- array(as.numeric(problem$Y), dim = c(problem$n, rrr_output_dims))
  fit <- fit_function(
    X = X_array, Y = Y_array, R = rank, lambda = 0,
    annealIter = 0L, convThresh = 0, seed = seed, maxIter = iterations
  )
  list(
    coefficient = array(as.numeric(fit$B), dim = problem$public_dims),
    coefficient_history = coefficient_history
  )
}

# The data of one fit. In the timing studies only run_regression_fit() is
# timed.
prepare_regression_fit <- function(method, X, Y, input_dims, output_dims,
                                   rank, start = NULL, iterations,
                                   seed = 1L, step_size = 1, blocks = rank,
                                   real_data = FALSE, record_path = TRUE) {
  list(method = method,
       problem = .regression_problem(X, Y, input_dims, output_dims),
       rank = rank, start = start, iterations = iterations, seed = seed,
       step_size = step_size, blocks = blocks, real_data = real_data,
       record_path = record_path)
}

run_regression_fit <- function(prepared) {
  with(prepared, switch(
    method,
    RGD = , RGN = .fit_segrecp(method, problem, rank, start, iterations,
                               step_size, record_path),
    `CP-ALS` = .fit_cp_als(problem, rank, start, iterations, record_path),
    RRR = .fit_rrr(problem, rank, iterations, seed, record_path),
    `Tucker-RGD` = .fit_tucker_rgd(problem, rank, iterations, real_data,
                                   record_path),
    `TT-IHT` = .fit_tt_iht(problem, rank, iterations, real_data, record_path),
    HOPLS = .fit_hopls(problem, blocks)
  ))
}

fit_regression_method <- function(method, ...) {
  run_regression_fit(prepare_regression_fit(method, ...))
}
