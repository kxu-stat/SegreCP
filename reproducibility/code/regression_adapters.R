# Fitting routines for the tensor-on-tensor regression methods of Section 5
# and Appendix C: CP-ALS, RGD, RGN (SegreCP), RRR (MultiwayRegression),
# Tucker-RGD, TT-IHT and HOPLS, and the shared CPCA initialization of CP-ALS,
# RGD and RGN. Coefficients use input modes followed by output modes; X is an
# n x prod(input_dims) matrix and Y is sample-first. Sourcing this file defines
# functions only.

.regression_whole_number <- function(x, name, allow_zero = FALSE) {
  lower <- if (allow_zero) 0 else 1
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < lower) {
    stop(sprintf("`%s` must be one %s integer.", name,
                 if (allow_zero) "nonnegative" else "positive"),
         call. = FALSE)
  }
  as.integer(x)
}

.regression_dims <- function(x, name, allow_empty = FALSE) {
  if (allow_empty && !length(x)) return(integer())
  if (!length(x) || !is.numeric(x) || any(!is.finite(x)) ||
      any(x != floor(x)) || any(x < 1)) {
    stop(sprintf("`%s` must contain positive integers.", name), call. = FALSE)
  }
  as.integer(x)
}

.regression_problem <- function(X, Y, input_dims, output_dims) {
  input_dims <- .regression_dims(input_dims, "input_dims")
  output_dims <- .regression_dims(output_dims, "output_dims", allow_empty = TRUE)
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("`X` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(X) != prod(input_dims)) {
    stop("The number of columns in `X` must equal prod(input_dims).",
         call. = FALSE)
  }
  n <- nrow(X)
  q <- if (length(output_dims)) prod(output_dims) else 1L
  if (is.null(dim(Y))) {
    if (q != 1L || !is.numeric(Y) || length(Y) != n) {
      stop("A vector `Y` is valid only for a scalar response of length n.",
           call. = FALSE)
    }
    Y_matrix <- matrix(as.numeric(Y), nrow = n, ncol = 1L)
  } else {
    expected <- if (length(output_dims)) c(n, output_dims) else c(n, 1L)
    if (!is.numeric(Y) || any(!is.finite(Y)) ||
        (!identical(as.integer(dim(Y)), as.integer(expected)) &&
         !(is.matrix(Y) && identical(dim(Y), c(n, q))))) {
      stop("`Y` must be sample-first with dimensions c(n, output_dims).",
           call. = FALSE)
    }
    Y_matrix <- matrix(as.numeric(Y), nrow = n, ncol = q)
  }
  if (any(!is.finite(Y_matrix))) stop("`Y` must be finite.", call. = FALSE)
  list(
    X = unname(X), Y = unname(Y_matrix), n = n,
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
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) ||
      ncol(X) != prod(input_dims)) {
    stop("`newX` must be finite with prod(input_dims) columns.", call. = FALSE)
  }
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

.option_value <- function(options, name, default = NULL, method = NULL) {
  value <- options[[name]]
  if (is.null(value)) return(default)
  if (!is.null(method) && is.list(value) && !is.null(value[[method]])) {
    return(value[[method]])
  }
  if (!is.null(method) && !is.list(value) && length(value) > 1L &&
      !is.null(names(value)) && method %in% names(value)) {
    return(unname(value[[method]]))
  }
  value
}

.regression_diagnostics <- function(options, method) {
  value <- .option_value(options, "diagnostics", TRUE, method)
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`options$diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.unwrap_shared_initialization <- function(shared_init, problem, rank) {
  components <- if (inherits(shared_init, "regression_cpca_initialization")) {
    shared_init$components
  } else if (is.list(shared_init) && !is.null(shared_init$components)) {
    shared_init$components
  } else shared_init
  if (!is.list(components) || length(components) != rank) {
    stop("`shared_init` must contain exactly `rank` CP components.",
         call. = FALSE)
  }
  valid <- vapply(components, function(component) {
    is.list(component) && length(component$lambda) == 1L &&
      is.finite(component$lambda) && length(component$U) ==
      length(problem$public_dims) &&
      all(vapply(seq_along(component$U), function(mode) {
        is.numeric(component$U[[mode]]) &&
          length(component$U[[mode]]) == problem$public_dims[[mode]] &&
          all(is.finite(component$U[[mode]]))
      }, logical(1L)))
  }, logical(1L))
  if (!all(valid)) stop("`shared_init` has invalid component dimensions.",
                        call. = FALSE)
  components
}

make_shared_regression_initialization <- function(
    X, Y, input_dims, output_dims = integer(), rank,
    split = NULL, design_variance = 1,
    normalize_regression = c("none", "sqrt_n"), seed = NULL,
    init_args = list(), method = c("CPCA", "HOSVD", "Random")) {
  method <- match.arg(method)
  if (!requireNamespace("SegreCP", quietly = TRUE)) {
    stop("SegreCP is required to compute the shared CPCA initialization.",
         call. = FALSE)
  }
  problem <- .regression_problem(X, Y, input_dims, output_dims)
  rank <- .regression_whole_number(rank, "rank")
  normalize_regression <- match.arg(normalize_regression)
  if (length(design_variance) != 1L || !is.finite(design_variance) ||
      design_variance <= 0) {
    stop("`design_variance` must be positive and finite.", call. = FALSE)
  }
  scale <- if (normalize_regression == "sqrt_n") 1 / sqrt(problem$n) else 1
  moment <- crossprod(problem$X * scale, problem$Y * scale) / design_variance
  moment <- array(as.numeric(moment), dim = problem$public_dims)
  arguments <- c(
    list(x = moment, rank = rank, method = tolower(method), split = split, seed = seed),
    init_args
  )
  initialization <- do.call(SegreCP::cp_initialize, arguments)
  structure(list(
    components = initialization$comps,
    initialization = initialization[names(initialization) != "comps"],
    moment = moment,
    metadata = list(
      method = method, split = split, design_variance = design_variance,
      normalization = normalize_regression, scale = scale,
      mode_order = "input_then_output"
    )
  ), class = if (method == "CPCA") "regression_cpca_initialization" else
    "regression_hosvd_initialization")
}

# Preserve the original public experiment helper and its return class.
make_shared_cpca_initialization <- function(
    X, Y, input_dims, output_dims = integer(), rank,
    split = NULL, design_variance = 1,
    normalize_regression = c("none", "sqrt_n"), seed = NULL,
    init_args = list()) {
  make_shared_regression_initialization(X, Y, input_dims, output_dims, rank,
    split, design_variance, match.arg(normalize_regression), seed, init_args,
    method = "CPCA")
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

.modewise_design_fibres_r <- function(X, factors) {
  n <- nrow(X)
  input_dims <- vapply(factors, length, integer(1L))
  X_array <- array(as.numeric(X), dim = c(n, input_dims))
  lapply(seq_along(factors), function(mode) {
    array_mode <- mode + 1L
    permutation <- c(1L, array_mode,
                     setdiff(seq_len(length(input_dims) + 1L), c(1L, array_mode)))
    unfolded <- matrix(aperm(X_array, permutation),
                       nrow = n * input_dims[[mode]])
    contraction <- .kron_vector(factors[-mode])
    matrix(unfolded %*% contraction, nrow = n,
           ncol = input_dims[[mode]])
  })
}

.modewise_design_fibres <- function(X, factors) {
  if (requireNamespace("SegreCP", quietly = TRUE)) {
    namespace <- asNamespace("SegreCP")
    if (exists("modewise_design_fibres_cpp", envir = namespace,
               inherits = FALSE)) {
      return(get("modewise_design_fibres_cpp", envir = namespace,
                 inherits = FALSE)(X, factors))
    }
  }
  .modewise_design_fibres_r(X, factors)
}

.objective_path_from_components <- function(history, problem) {
  objective <- vapply(history, function(components) {
    coefficient <- .components_to_array(components, problem$public_dims)
    residual <- .predict_coefficient_matrix(
      coefficient, problem$X, problem$input_dims, problem$output_dims
    ) - problem$Y
    0.5 * sum(residual^2)
  }, numeric(1L))
  data.frame(iteration = seq_along(objective) - 1L, objective = objective)
}

.fit_segrecp_method <- function(method, problem, rank, shared_init, iterations,
                                ridge, options, seed) {
  if (!requireNamespace("SegreCP", quietly = TRUE)) {
    stop("SegreCP is required for RGD and RGN.", call. = FALSE)
  }
  components <- .unwrap_shared_initialization(shared_init, problem, rank)
  diagnostics <- .regression_diagnostics(options, method)
  profile <- .option_value(options, "profile", "manuscript", method)
  tolerance <- .option_value(options, "tolerance", 0, method)
  step_size <- .option_value(options, "step_size", 0.5, method)
  default_rule <- if (method == "RGD") "exact_quadratic" else "fixed"
  step_rule <- .option_value(options, "step_rule", default_rule, method)
  component_order <- .option_value(
    options, "component_order", "decreasing_weight", method
  )
  normalization <- .option_value(
    options, "normalize_regression", "none", method
  )
  design_variance <- .option_value(options, "design_variance", 1, method)
  control <- SegreCP::segre_control(
    profile = profile, max_iter = iterations, tolerance = tolerance,
    step_size = step_size, step_rule = step_rule,
    component_order = component_order, ridge = ridge,
    normalize_regression = normalization,
    trace = isTRUE(.option_value(options, "trace", FALSE, method)),
    armijo = .option_value(options, "armijo", 1e-4, method),
    backtrack_factor = .option_value(options, "backtrack_factor", 0.5, method),
    max_backtracks = .option_value(options, "max_backtracks", 20L, method)
  )
  response <- .response_from_matrix(problem$Y, problem$output_dims)
  fit <- SegreCP::segre_regress(
    X = problem$X, Y = response, rank = rank,
    input_dims = problem$input_dims,
    output_dims = if (length(problem$output_dims)) problem$output_dims else NULL,
    method = tolower(method), init = components, seed = seed,
    design_variance = design_variance, control = control,
    keep_history = diagnostics
  )
  coefficient <- stats::coef(fit)
  native_objective <- data.frame(iteration = seq_along(fit$objective) - 1L,
                                 objective = fit$objective)
  if (diagnostics) {
    objective <- .objective_path_from_components(fit$component_history, problem)
  } else {
    objective <- native_objective
    objective$objective <- objective$objective / fit$normalization^2
  }
  if (method == "RGD") {
    accepted_component <- is.finite(fit$accepted_steps) &
      fit$accepted_steps > 0
    accepted_iterations <- sum(rowSums(accepted_component) == rank)
    accepted_component_updates <- sum(accepted_component)
    accepted_definition <- paste(
      "outer iterations with a positive accepted step for every component;",
      "component updates are also recorded"
    )
  } else {
    accepted_iterations <- as.integer(fit$iterations)
    accepted_component_updates <- as.integer(fit$iterations * rank)
    accepted_definition <- paste(
      "RGN has no line-search acceptance filter; each completed outer",
      "iteration updates every component"
    )
  }
  list(
    coefficient = coefficient, components = fit$components,
    component_history = fit$component_history,
    normal_condition = fit$normal_condition,
    objective = objective, native_objective = native_objective,
    actual_iterations = as.integer(fit$iterations),
    accepted_iterations = as.integer(accepted_iterations),
    stop_reason = as.character(fit$stopping_reason), stopped_on_loss = FALSE,
    metadata = list(
      engine = "SegreCP::segre_regress", profile = profile,
      initialization = "shared_CPCA", normalization = normalization,
      objective_data_scale = "original_X_Y",
      ridge = ridge,
      ridge_role = if (method == "RGN") "tangent_normal_system" else "unused",
      regularization = if (method == "RGN") ridge else 0,
      control = unclass(fit$control), design_variance = design_variance,
      accepted_component_steps = fit$accepted_steps,
      accepted_component_updates = accepted_component_updates,
      accepted_count_definition = accepted_definition,
      native_objective_scale = if (normalization == "sqrt_n")
        "X_and_Y_divided_by_sqrt_n" else "original_X_Y",
      objective_increased = any(diff(objective$objective) >
        1e-10 * pmax(1, head(abs(objective$objective), -1L)))
    )
  )
}

.regression_collect_coefficients <- function(options, method) {
  value <- .option_value(options, "collect_coefficients", FALSE, method)
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop("`options$collect_coefficients` must be TRUE or FALSE.", call. = FALSE)
  }
  value
}

.fit_cp_als <- function(problem, rank, shared_init, iterations, ridge,
                        diagnostics = TRUE, collect_coefficients = FALSE) {
  public_components <- .unwrap_shared_initialization(shared_init, problem, rank)
  components <- .public_to_internal_components(public_components, problem)
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
  coefficient_history <- if (collect_coefficients) list(current_coefficient()) else NULL

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
  objective <- if (diagnostics) numeric(iterations + 1L) else NULL
  if (diagnostics) objective[[1L]] <- 0.5 * sum((problem$Y - prediction())^2)

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
          fibres <- .modewise_design_fibres(problem$X, input_factors)
          columns <- ((component - 1L) * mode_size + 1L):
            (component * mode_size)
          fibre_stack[, columns] <- fibres[[input_position]]
        }
        cross_gram <- crossprod(fibre_stack)
        scale <- (weights %o% weights) * output_gram
        normal <- cross_gram * kronecker(scale, matrix(1, mode_size, mode_size))
        if (ridge > 0) normal <- normal + ridge * diag(mode_size * rank)
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
        if (ridge > 0) scale <- scale + ridge * diag(rank)
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
      if (any(!is.finite(norms)) || any(norms <= .Machine$double.eps)) {
        stop("CP-ALS produced a zero or nonfinite factor column.", call. = FALSE)
      }
      factors[[mode]] <- sweep(updated, 2L, norms, "/")
      weights <- weights * norms
    }
    if (diagnostics) {
      objective[[iteration + 1L]] <- 0.5 * sum((problem$Y - prediction())^2)
    }
    if (collect_coefficients) coefficient_history[[iteration + 1L]] <- current_coefficient()
  }
  components <- lapply(seq_len(rank), function(component) {
    list(lambda = weights[[component]], U = lapply(factors, function(A) {
      A[, component]
    }))
  })
  public_components <- .internal_to_public_components(components, problem)
  coefficient <- .components_to_array(public_components, problem$public_dims)
  list(
    coefficient = coefficient, components = public_components,
    coefficient_history = coefficient_history,
    objective = if (diagnostics)
      data.frame(iteration = 0:iterations, objective = objective) else NULL,
    native_objective = NULL, actual_iterations = iterations,
    accepted_iterations = iterations, stop_reason = "iteration_budget",
    stopped_on_loss = FALSE,
    metadata = list(
      engine = "local_CP_ALS_block_updates", initialization = "shared_CPCA",
      normalization = "none", objective_data_scale = "original_X_Y",
      ridge = ridge, ridge_role = "block_normal_system",
      regularization = ridge,
      control = list(max_iter = iterations, tolerance = 0, ridge = ridge,
                     update_scheme = "Gauss-Seidel"),
      objective_increased = if (diagnostics) any(diff(objective) >
        1e-10 * pmax(1, head(abs(objective), -1L))) else NA
    )
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

.tucker_projection_stabilizers <- function(state, projection_ridge) {
  vapply(seq_along(state$factors), function(mode) {
    G_mode <- .mode_unfold(state$core, mode)
    gram <- G_mode %*% t(G_mode)
    projection_ridge * max(sum(diag(gram)) / nrow(gram), 1)
  }, numeric(1L))
}

.tucker_tangent_projection <- function(Z, state, projection_ridge = 0) {
  core <- state$core
  factors <- state$factors
  order <- length(factors)
  stabilizers <- .tucker_projection_stabilizers(state, projection_ridge)
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
    stabilizer <- stabilizers[[mode]]
    factor <- factors[[mode]]
    dot_factor <- Z_mode %*% t(G_mode) %*%
      solve(gram + stabilizer * diag(nrow(gram)))
    dot_factor <- dot_factor - factor %*% crossprod(factor, dot_factor)
    tangent_factors <- factors
    tangent_factors[[mode]] <- dot_factor
    projected <- projected + .tucker_reconstruct(core, tangent_factors)
  }
  attr(projected, "projection_stabilizers") <- stabilizers
  projected
}

.fit_tucker_rgd <- function(problem, rank, iterations, ridge, options) {
  diagnostics <- .regression_diagnostics(options, "Tucker-RGD")
  collect_coefficients <- .regression_collect_coefficients(options, "Tucker-RGD")
  protocol <- .option_value(options, "protocol", "synthetic_grid",
                            "Tucker-RGD")
  if (!protocol %in% c("synthetic_grid", "real_trade")) {
    stop("`options$protocol` must be `synthetic_grid` or `real_trade`.",
         call. = FALSE)
  }
  projection_ridge <- .option_value(
    options, "tucker_projection_ridge",
    0, "Tucker-RGD"
  )
  if (length(projection_ridge) != 1L || !is.numeric(projection_ridge) ||
      !is.finite(projection_ridge) || projection_ridge < 0) {
    stop("`options$tucker_projection_ridge` must be one nonnegative finite number.",
         call. = FALSE)
  }
  ranks <- .option_value(
    options, "tucker_ranks", rep(rank, length(problem$internal_dims)),
    "Tucker-RGD"
  )
  ranks <- .regression_dims(ranks, "options$tucker_ranks")
  if (length(ranks) != length(problem$internal_dims)) {
    stop("Tucker ranks must have one entry per coefficient mode.", call. = FALSE)
  }
  if (protocol == "real_trade") {
    normal <- crossprod(problem$X) + ridge * diag(problem$p_in)
    beta <- solve(normal, crossprod(problem$X, problem$Y))
    initialization <- array(as.numeric(t(beta)), dim = problem$internal_dims)
  } else {
    initialization <- .internal_adjoint(problem$Y, problem$X, problem)
  }
  state <- .tucker_hosvd(initialization, ranks)
  coefficient_internal <- state$recon
  coefficient_history <- if (collect_coefficients)
    list(.internal_to_public_array(coefficient_internal, problem)) else NULL
  data_loss <- function(A) {
    0.5 * sum((.apply_internal_operator(A, problem$X, problem) - problem$Y)^2)
  }
  native_loss <- function(A) data_loss(A) +
    if (protocol == "real_trade") 0.5 * ridge * sum(A^2) else 0
  objective <- if (diagnostics) numeric(iterations + 1L) else NULL
  native_objective <- numeric(iterations + 1L)
  projection_stabilizers <- if (diagnostics) matrix(
    NA_real_, nrow = iterations, ncol = length(problem$internal_dims)
  ) else NULL
  if (diagnostics) objective[[1L]] <- data_loss(coefficient_internal)
  native_objective[[1L]] <- native_loss(coefficient_internal)
  accepted <- 0L
  max_backtracks <- .regression_whole_number(
    .option_value(options, "max_backtracks", 20L, "Tucker-RGD"),
    "options$max_backtracks"
  )
  for (iteration in seq_len(iterations)) {
    residual <- .apply_internal_operator(
      coefficient_internal, problem$X, problem
    ) - problem$Y
    ambient <- .internal_adjoint(residual, problem$X, problem)
    if (protocol == "real_trade") {
      ambient <- ambient + ridge * coefficient_internal
    }
    gradient <- .tucker_tangent_projection(
      ambient, state, projection_ridge = projection_ridge
    )
    if (diagnostics) projection_stabilizers[iteration, ] <-
      attr(gradient, "projection_stabilizers")
    measured <- .apply_internal_operator(gradient, problem$X, problem)
    denominator <- sum(measured^2)
    if (protocol == "real_trade") {
      denominator <- denominator + ridge * sum(gradient^2)
    }
    step <- sum(gradient^2) / max(denominator, .Machine$double.eps)
    candidate <- coefficient_internal
    candidate_state <- NULL
    candidate_loss <- native_objective[[iteration]]
    for (backtrack in seq_len(max_backtracks)) {
      trial_state <- .tucker_hosvd(
        coefficient_internal - step * gradient, ranks
      )
      trial_loss <- native_loss(trial_state$recon)
      if (is.finite(trial_loss) &&
          trial_loss <= native_objective[[iteration]]) {
        candidate <- trial_state$recon
        candidate_state <- trial_state
        candidate_loss <- trial_loss
        accepted <- accepted + 1L
        break
      }
      step <- step / 2
    }
    coefficient_internal <- candidate
    if (!is.null(candidate_state)) state <- candidate_state
    if (diagnostics) objective[[iteration + 1L]] <- data_loss(coefficient_internal)
    native_objective[[iteration + 1L]] <- candidate_loss
    if (collect_coefficients) coefficient_history[[iteration + 1L]] <-
      .internal_to_public_array(coefficient_internal, problem)
  }
  coefficient <- .internal_to_public_array(coefficient_internal, problem)
  list(
    coefficient = coefficient, components = NULL,
    coefficient_history = coefficient_history,
    objective = if (diagnostics)
      data.frame(iteration = 0:iterations, objective = objective) else NULL,
    native_objective = data.frame(
      iteration = 0:iterations, objective = native_objective
    ), actual_iterations = iterations,
    accepted_iterations = accepted, stop_reason = "iteration_budget",
    stopped_on_loss = FALSE,
    metadata = list(
      engine = "local_Tucker_RGD",
      initialization = if (protocol == "real_trade" && ridge > 0)
        "native_ridge_OLS_T_HOSVD" else if (protocol == "real_trade")
        "native_OLS_T_HOSVD" else "native_adjoint_T_HOSVD",
      normalization = "none", objective_data_scale = "original_X_Y",
      ridge = ridge,
      ridge_role = if (protocol == "real_trade")
        "coefficient_penalty" else "unused",
      regularization = if (protocol == "real_trade") ridge else 0,
      control = list(max_iter = iterations, multilinear_rank = ranks,
                     step_rule = "exact_quadratic_with_backtracking",
                     max_backtracks = max_backtracks, protocol = protocol,
                     tucker_projection_ridge = projection_ridge),
      tucker_projection_ridge = projection_ridge,
      tucker_projection_ridge_scaling = paste(
        "per-mode diagonal addition equals tucker_projection_ridge times",
        "max(trace(G_mode %*% t(G_mode)) / nrow(G_mode), 1)"
      ),
      tucker_projection_stabilizers = projection_stabilizers,
      objective_increased = if (diagnostics) any(diff(objective) >
        1e-10 * pmax(1, head(abs(objective), -1L))) else NA
    )
  )
}

.tt_svd_project <- function(A, tt_ranks) {
  dimensions <- dim(A)
  order <- length(dimensions)
  if (length(tt_ranks) != order - 1L) {
    stop("TT ranks must have order - 1 entries.", call. = FALSE)
  }
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

.fit_tt_iht <- function(problem, rank, iterations, ridge, options) {
  diagnostics <- .regression_diagnostics(options, "TT-IHT")
  collect_coefficients <- .regression_collect_coefficients(options, "TT-IHT")
  protocol <- .option_value(options, "protocol", "synthetic_grid", "TT-IHT")
  if (!protocol %in% c("synthetic_grid", "real_trade")) {
    stop("`options$protocol` must be `synthetic_grid` or `real_trade`.",
         call. = FALSE)
  }
  tt_ranks <- .option_value(
    options, "tt_ranks", rep(rank, length(problem$internal_dims) - 1L),
    "TT-IHT"
  )
  tt_ranks <- .regression_dims(tt_ranks, "options$tt_ranks")
  adjoint <- .internal_adjoint(problem$Y, problem$X, problem)
  step <- if (protocol == "real_trade") {
    largest <- max(eigen(crossprod(problem$X), symmetric = TRUE,
                         only.values = TRUE)$values)
    1 / max(largest, 1e-12)
  } else 1
  coefficient_internal <- .tt_svd_project(step * adjoint, tt_ranks)
  coefficient_history <- if (collect_coefficients)
    list(.internal_to_public_array(coefficient_internal, problem)) else NULL
  objective <- if (diagnostics) numeric(iterations + 1L) else NULL
  if (diagnostics) objective[[1L]] <- 0.5 * sum((.apply_internal_operator(
    coefficient_internal, problem$X, problem
  ) - problem$Y)^2)
  for (iteration in seq_len(iterations)) {
    residual <- .apply_internal_operator(
      coefficient_internal, problem$X, problem
    ) - problem$Y
    gradient <- .internal_adjoint(residual, problem$X, problem)
    coefficient_internal <- .tt_svd_project(
      coefficient_internal - step * gradient, tt_ranks
    )
    if (diagnostics) objective[[iteration + 1L]] <- 0.5 * sum((.apply_internal_operator(
      coefficient_internal, problem$X, problem
    ) - problem$Y)^2)
    if (collect_coefficients) coefficient_history[[iteration + 1L]] <-
      .internal_to_public_array(coefficient_internal, problem)
  }
  coefficient <- .internal_to_public_array(coefficient_internal, problem)
  list(
    coefficient = coefficient, components = NULL,
    coefficient_history = coefficient_history,
    objective = if (diagnostics)
      data.frame(iteration = 0:iterations, objective = objective) else NULL,
    native_objective = NULL, actual_iterations = iterations,
    accepted_iterations = iterations, stop_reason = "iteration_budget",
    stopped_on_loss = FALSE,
    metadata = list(
      engine = "local_TT_IHT", initialization = "native_TT_SVD_adjoint",
      normalization = "none", objective_data_scale = "original_X_Y",
      ridge = ridge, ridge_role = "unused", regularization = 0,
      control = list(max_iter = iterations, tt_ranks = tt_ranks,
                     step_size = step, protocol = protocol),
      objective_increased = if (diagnostics) any(diff(objective) >
        1e-10 * pmax(1, head(abs(objective), -1L))) else NA
    )
  )
}

.rank_one_hooi <- function(A, max_iter, tolerance) {
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
  value <- NA_real_
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
  list(factors = factors, value = value, iterations = iteration)
}

.fit_hopls <- function(problem, rank, options) {
  diagnostics <- .regression_diagnostics(options, "HOPLS")
  blocks <- .regression_whole_number(
    .option_value(options, "hopls_blocks", rank, "HOPLS"),
    "options$hopls_blocks"
  )
  hooi_max <- .regression_whole_number(
    .option_value(options, "hopls_hooi_max_iter", 100L, "HOPLS"),
    "options$hopls_hooi_max_iter"
  )
  hooi_tolerance <- .option_value(
    options, "hopls_hooi_tolerance", 1e-10, "HOPLS"
  )
  X_residual <- problem$X
  Y_residual <- problem$Y
  input_loadings <- vector("list", blocks)
  output_loadings <- vector("list", blocks)
  response_core <- numeric(blocks)
  effective_weights <- vector("list", blocks)
  hooi_iterations <- integer(blocks)
  fitted_from_scores <- if (diagnostics) matrix(0, problem$n, problem$q_out) else NULL
  objective <- if (diagnostics) numeric(blocks + 1L) else NULL
  if (diagnostics) objective[[1L]] <- 0.5 * sum(problem$Y^2)
  output_modes <- if (problem$output_order) seq_len(problem$output_order) else integer()
  input_modes <- problem$output_order + seq_len(problem$input_order)

  for (component in seq_len(blocks)) {
    cross_covariance <- array(
      as.numeric(crossprod(Y_residual, X_residual)),
      dim = problem$internal_dims
    )
    hooi <- .rank_one_hooi(cross_covariance, hooi_max, hooi_tolerance)
    output_loading <- if (problem$output_order) {
      .kron_vector(hooi$factors[output_modes])
    } else 1
    input_loading <- .kron_vector(hooi$factors[input_modes])
    raw_score <- as.numeric(X_residual %*% input_loading)
    score_norm <- sqrt(sum(raw_score^2))
    if (!is.finite(score_norm) || score_norm <= .Machine$double.eps) {
      stop("HOPLS encountered a zero or nonfinite predictor score.",
           call. = FALSE)
    }
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
    hooi_iterations[[component]] <- hooi$iterations
    block_fit <- tcrossprod(
      score * response_core[[component]], output_loading
    )
    if (diagnostics) {
      fitted_from_scores <- fitted_from_scores + block_fit
      objective[[component + 1L]] <- 0.5 * sum(
        (problem$Y - fitted_from_scores)^2
      )
    }
    X_residual <- X_residual - tcrossprod(raw_score, input_loading)
    Y_residual <- Y_residual - block_fit
  }
  coefficient_matrix <- Reduce(`+`, lapply(seq_len(blocks), function(component) {
    tcrossprod(
      effective_weights[[component]],
      response_core[[component]] * output_loadings[[component]]
    )
  }))
  coefficient <- array(as.numeric(coefficient_matrix), dim = problem$public_dims)
  representation_error <- if (diagnostics) sqrt(sum((
    .predict_coefficient_matrix(
      coefficient, problem$X, problem$input_dims, problem$output_dims
    ) - fitted_from_scores
  )^2)) / max(sqrt(sum(fitted_from_scores^2)), .Machine$double.eps) else NA_real_
  list(
    coefficient = coefficient, components = NULL,
    objective = if (diagnostics)
      data.frame(iteration = 0:blocks, objective = objective) else NULL,
    native_objective = NULL, actual_iterations = NA_integer_,
    accepted_iterations = NA_integer_, stop_reason = "latent_block_budget",
    stopped_on_loss = FALSE,
    metadata = list(
      engine = "local_HOPLS_rank_one_loadings", initialization = "native_HOSVD_HOOI",
      normalization = "none", objective_data_scale = "original_X_Y",
      ridge = 0, ridge_role = "unused", regularization = 0,
      latent_blocks = blocks, hooi_iterations = hooi_iterations,
      representation_error = representation_error,
      objective_index = "latent_block",
      control = list(latent_blocks = blocks, hooi_max_iter = hooi_max,
                     hooi_tolerance = hooi_tolerance),
      objective_increased = if (diagnostics) any(diff(objective) >
        1e-10 * pmax(1, head(abs(objective), -1L))) else NA
    )
  )
}

.make_rrr_fixed <- function(max_iterations, coefficient_callback = NULL) {
  if (!requireNamespace("MultiwayRegression", quietly = TRUE)) {
    stop("The RRR adapter requires the MultiwayRegression package.",
         call. = FALSE)
  }
  fit_function <- MultiwayRegression::rrr
  function_formals <- formals(fit_function)
  function_formals$maxIter <- max_iterations
  formals(fit_function) <- function_formals
  patched <- 0L
  patch_loop <- function(expression) {
    if (is.call(expression) && identical(expression[[1L]], as.name("while")) &&
        identical(paste(deparse(expression[[2L]]), collapse = ""),
                  "conv > convThresh")) {
      old_condition <- expression[[2L]]
      expression[[2L]] <- substitute(
        OLD && j < maxIter, list(OLD = old_condition)
      )
      patched <<- patched + 1L
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
  if (patched != 1L) stop("Could not patch the unique RRR outer loop.",
                          call. = FALSE)
  if (!is.null(coefficient_callback)) {
    if (!is.function(coefficient_callback)) stop("RRR coefficient callback must be a function.")
    callback_environment <- new.env(parent = environment(fit_function))
    callback_environment$.rrr_coefficient_callback <- coefficient_callback
    environment(fit_function) <- callback_environment
  }
  fit_function
}

.fit_rrr <- function(problem, rank, iterations, ridge, seed, options) {
  collect_coefficients <- .regression_collect_coefficients(options, "RRR")
  coefficient_history <- if (collect_coefficients) list() else NULL
  callback <- if (collect_coefficients) function(iteration, coefficient) {
    coefficient_history[[as.integer(iteration) + 1L]] <<-
      array(as.numeric(coefficient), dim = problem$public_dims)
  } else NULL
  fit_function <- .make_rrr_fixed(iterations, coefficient_callback = callback)
  protocol <- .option_value(options, "protocol", "synthetic_grid", "RRR")
  if (!protocol %in% c("synthetic_grid", "real_trade")) {
    stop("`options$protocol` must be `synthetic_grid` or `real_trade`.",
         call. = FALSE)
  }
  rrr_input_dims <- if (protocol == "real_trade" &&
                        length(problem$input_dims) == 1L) {
    c(problem$input_dims, 1L)
  } else problem$input_dims
  X_array <- array(as.numeric(problem$X), dim = c(problem$n, rrr_input_dims))
  rrr_output_dims <- if (length(problem$output_dims)) problem$output_dims else 1L
  Y_array <- array(as.numeric(problem$Y), dim = c(problem$n, rrr_output_dims))
  anneal <- .regression_whole_number(
    .option_value(options, "rrr_anneal", 0L, "RRR"),
    "options$rrr_anneal", allow_zero = TRUE
  )
  rrr_seed <- .regression_whole_number(
    .option_value(options, "rrr_seed", seed, "RRR"), "options$rrr_seed"
  )
  convergence_threshold <- .option_value(
    options, "rrr_tolerance", 0, "RRR"
  )
  fit <- fit_function(
    X = X_array, Y = Y_array, R = rank, lambda = ridge,
    annealIter = anneal, convThresh = convergence_threshold, seed = rrr_seed,
    maxIter = iterations
  )
  blocks_per_iteration <- length(rrr_input_dims) + length(rrr_output_dims)
  actual <- length(fit$sse) / blocks_per_iteration
  if (!is.finite(actual) || actual != floor(actual)) {
    stop("RRR returned an incomplete outer iteration.", call. = FALSE)
  }
  actual <- as.integer(actual)
  cycle_end <- seq.int(blocks_per_iteration, length(fit$sse),
                       by = blocks_per_iteration)
  objective <- data.frame(
    iteration = seq_len(actual),
    objective = 0.5 * as.numeric(fit$sse[cycle_end])
  )
  coefficient_with_native_dims <- array(
    as.numeric(fit$B), dim = c(rrr_input_dims, rrr_output_dims)
  )
  coefficient <- array(as.numeric(coefficient_with_native_dims),
                       dim = c(problem$input_dims, rrr_output_dims))
  if (!length(problem$output_dims)) {
    coefficient <- array(as.numeric(coefficient), dim = problem$input_dims)
  }
  list(
    coefficient = coefficient, components = NULL,
    coefficient_history = coefficient_history,
    objective = objective, native_objective = data.frame(
      block_update = seq_along(fit$sse), sse = as.numeric(fit$sse),
      penalized_sse = as.numeric(fit$sseR)
    ),
    actual_iterations = actual, accepted_iterations = actual,
    stop_reason = if (actual < iterations) "native_convergence" else
      "iteration_budget",
    stopped_on_loss = FALSE,
    metadata = list(
      engine = "MultiwayRegression::rrr_with_max_iteration_guard",
      initialization = "native_random", normalization = "none",
      objective_data_scale = "original_X_Y",
      ridge = ridge, ridge_role = "coefficient_penalty",
      regularization = ridge,
      control = list(
        max_iter = iterations,
        convergence_threshold = convergence_threshold,
        anneal_iterations = anneal, lambda = ridge, protocol = protocol,
        native_input_dims = rrr_input_dims, native_seed = rrr_seed
      ),
      package_version = as.character(
        utils::packageVersion("MultiwayRegression")
      ),
      objective_increased = if (nrow(objective) > 1L) any(
        diff(objective$objective) >
          1e-10 * pmax(1, head(abs(objective$objective), -1L))
      ) else FALSE
    )
  )
}

.finalize_regression_fit <- function(result, method, problem, rank, iterations,
                                     seed) {
  expected_dimensions <- problem$public_dims
  if (!identical(as.integer(dim(result$coefficient)),
                 as.integer(expected_dimensions))) {
    stop(method, " returned a coefficient with incorrect dimensions.",
         call. = FALSE)
  }
  fitted_matrix <- .predict_coefficient_matrix(
    result$coefficient, problem$X, problem$input_dims, problem$output_dims
  )
  final_objective <- 0.5 * sum((fitted_matrix - problem$Y)^2)
  if (is.null(result$objective)) {
    result$objective <- data.frame(
      iteration = if (method == "HOPLS") result$metadata$latent_blocks else
        result$actual_iterations,
      objective = final_objective
    )
  }
  result$metadata$input_dims <- problem$input_dims
  result$metadata$output_dims <- problem$output_dims
  result$metadata$rank <- rank
  result$metadata$seed <- seed
  result$metadata$requested_iterations <- iterations
  result$metadata$common_final_objective <- final_objective
  result$metadata$mode_order <- "input_then_output"
  result$method <- method
  result$fitted <- .response_from_matrix(fitted_matrix, problem$output_dims)
  result$requested_iterations <- iterations
  result$finite <- all(is.finite(result$coefficient)) &&
    all(is.finite(fitted_matrix)) &&
    is.data.frame(result$objective) &&
    nrow(result$objective) > 0L &&
    all(is.finite(result$objective$objective))
  coefficient <- result$coefficient
  input_dims <- problem$input_dims
  output_dims <- problem$output_dims
  result$predict <- local({
    coefficient_local <- coefficient
    input_dims_local <- input_dims
    output_dims_local <- output_dims
    function(newX) {
      matrix_prediction <- .predict_coefficient_matrix(
        coefficient_local, newX, input_dims_local, output_dims_local
      )
      .response_from_matrix(matrix_prediction, output_dims_local)
    }
  })
  class(result) <- c("current_regression_fit", "list")
  result
}

# Prepare outside a wall-clock timer. Shared initialization is supplied by the
# caller; each baseline's native initialization remains inside the kernel.
# The timing default disables diagnostics. fit_regression_method() keeps its
# historical diagnostics-on default unless the caller explicitly changes it.
prepare_regression_kernel <- function(
    method, X, Y, input_dims, output_dims = integer(), rank,
    shared_init = NULL, iterations = 10L, seed = 1L, ridge = 0,
    options = list()) {
  methods <- c("RGD", "RGN", "CP-ALS", "RRR", "Tucker-RGD", "TT-IHT",
               "HOPLS")
  if (length(method) != 1L || !method %in% methods) {
    stop("`method` must be one of: ", paste(methods, collapse = ", "),
         call. = FALSE)
  }
  problem <- .regression_problem(X, Y, input_dims, output_dims)
  rank <- .regression_whole_number(rank, "rank")
  iterations <- .regression_whole_number(iterations, "iterations")
  seed <- .regression_whole_number(seed, "seed")
  if (length(ridge) != 1L || !is.numeric(ridge) || !is.finite(ridge) ||
      ridge < 0) {
    stop("`ridge` must be one nonnegative finite number.", call. = FALSE)
  }
  if (!is.list(options)) stop("`options` must be a named list.", call. = FALSE)
  if (is.null(options$diagnostics)) options$diagnostics <- FALSE
  diagnostics <- .regression_diagnostics(options, method)
  if (method %in% c("RGD", "RGN", "CP-ALS") && is.null(shared_init)) {
    stop(method, " requires the realized common `shared_init` object.",
         call. = FALSE)
  }

  initialization <- NULL
  shared_hash <- NA_character_
  if (method %in% c("CP-ALS", "RGD", "RGN")) {
    initialization <- if (!is.null(shared_init$metadata$method))
      paste0("shared_", shared_init$metadata$method) else "shared_custom"
    if (!requireNamespace("digest", quietly = TRUE)) {
      stop("The digest package is required to record the shared start.",
           call. = FALSE)
    }
    shared_hash <- digest::digest(
      .unwrap_shared_initialization(shared_init, problem, rank),
      algo = "sha256", serialize = TRUE
    )
  }
  structure(list(
    method = method, problem = problem, rank = rank, shared_init = shared_init,
    iterations = iterations, seed = seed, ridge = ridge, options = options,
    diagnostics = diagnostics, initialization = initialization,
    shared_init_sha256 = shared_hash
  ), class = "prepared_regression_kernel")
}

# Time this call only. It invokes SegreCP directly for RGD/RGN, with no copied
# SegreCP solver, and omits adapter endpoint prediction/metric calculations.
fit_regression_kernel <- function(prepared) {
  if (!inherits(prepared, "prepared_regression_kernel")) {
    stop("Call prepare_regression_kernel() before fit_regression_kernel().",
         call. = FALSE)
  }
  method <- prepared$method
  problem <- prepared$problem
  rank <- prepared$rank
  shared_init <- prepared$shared_init
  iterations <- prepared$iterations
  seed <- prepared$seed
  ridge <- prepared$ridge
  options <- prepared$options

  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  result <- switch(
    method,
    RGD = .fit_segrecp_method(
      method, problem, rank, shared_init, iterations, ridge, options, seed
    ),
    RGN = .fit_segrecp_method(
      method, problem, rank, shared_init, iterations, ridge, options, seed
    ),
    `CP-ALS` = .fit_cp_als(
      problem, rank, shared_init, iterations, ridge, prepared$diagnostics,
      .regression_collect_coefficients(options, method)
    ),
    RRR = .fit_rrr(problem, rank, iterations, ridge, seed, options),
    `Tucker-RGD` = .fit_tucker_rgd(
      problem, rank, iterations, ridge, options
    ),
    `TT-IHT` = .fit_tt_iht(problem, rank, iterations, ridge, options),
    HOPLS = .fit_hopls(problem, rank, options)
  )
  if (!is.null(prepared$initialization)) {
    result$metadata$initialization <- prepared$initialization
  }
  result$metadata$shared_init_sha256 <- prepared$shared_init_sha256
  result$metadata$diagnostics <- prepared$diagnostics
  result$metadata$objective_recording <- if (prepared$diagnostics)
    "diagnostic_trajectory" else if (is.null(result$objective))
      "endpoint_only_after_timer" else "native_solver_trajectory"
  result$method <- method
  result$requested_iterations <- iterations
  class(result) <- c("current_regression_kernel", "list")
  result
}

# Call after stopping the timer. This independently evaluates the common
# original-data-scale endpoint and returns the usual prediction/metric API.
finalize_regression_kernel <- function(result, prepared) {
  if (!inherits(result, "current_regression_kernel") ||
      !inherits(prepared, "prepared_regression_kernel") ||
      !identical(result$method, prepared$method) ||
      !identical(result$requested_iterations, prepared$iterations)) {
    stop("The kernel result and preparation must identify the same fit.",
         call. = FALSE)
  }
  .finalize_regression_fit(result, prepared$method, prepared$problem,
                           prepared$rank, prepared$iterations, prepared$seed)
}

fit_regression_method <- function(
    method, X, Y, input_dims, output_dims = integer(), rank,
    shared_init = NULL, iterations = 10L, seed = 1L, ridge = 0,
    options = list()) {
  if (!is.list(options)) stop("`options` must be a named list.", call. = FALSE)
  if (is.null(options$diagnostics)) options$diagnostics <- TRUE
  prepared <- prepare_regression_kernel(
    method, X, Y, input_dims, output_dims, rank, shared_init, iterations,
    seed, ridge, options
  )
  finalize_regression_kernel(fit_regression_kernel(prepared), prepared)
}

# Options of the regression simulations (Section 5 and Appendix C); the
# real-data script changes the protocol, ranks and RGD step.
regression_trajectory_options <- function() {
  list(protocol = "synthetic_grid", profile = "manuscript", tolerance = 0,
    # RGD: fixed step alpha = 1, i.e. alpha * sigma_X^2 = 1 for the simulation
    # design (entries with variance 1/n), as in the regression corollaries.
    # RGN does not use a step size.
    step_size = 1, step_rule = c(RGD = "fixed", RGN = "fixed"),
    component_order = "decreasing_weight", normalize_regression = "none", design_variance = 1,
    tucker_ranks = rep(3L, 5L), tt_ranks = rep(3L, 4L), hopls_blocks = 3L,
    tucker_projection_ridge = 0, rrr_anneal = 0L, rrr_tolerance = 0,
    diagnostics = TRUE, collect_coefficients = TRUE)
}

