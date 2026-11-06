# CP-decomposition comparators of Appendix C: Orth-ALS, ICO, NLS-GN, GN-CG,
# CP-ALS with enhanced line search and hybrid AMDM. Each is an R
# implementation of the cited algorithm (references below); all start from the
# CPCA components supplied by the caller. Sourcing this file defines functions
# only.

.dc_references <- c(
  "Orth-ALS" = "Sharan and Valiant (2017), orthogonalized ALS",
  "ICO" = "Han and Zhang (2023), iterative concurrent orthogonalization",
  "NLS-GN" = "Sorber, Van Barel and De Lathauwer (2013), NLS Gauss-Newton",
  "GN-CG" = "Singh et al. (2021), implicit Gauss-Newton with CG",
  "CP-ALS+ELS" = "Rajih, Comon and Harshman (2008), enhanced line search",
  "Hybrid AMDM" = "Singh and Solomonik (2023), hybrid AMDM"
)

.dc_whole_number <- function(x, name, positive = TRUE) {
  lower <- if (positive) 1 else 0
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < lower) {
    stop("`", name, "` must be one ",
         if (positive) "positive" else "nonnegative", " integer.",
         call. = FALSE)
  }
  as.integer(x)
}

.dc_scalar <- function(x, name, lower = -Inf, strict_lower = FALSE) {
  invalid_lower <- if (strict_lower) x <= lower else x < lower
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || invalid_lower) {
    stop("`", name, "` has an invalid value.", call. = FALSE)
  }
  as.numeric(x)
}

.dc_option <- function(options, name, default) {
  value <- options[[name]]
  if (is.null(value)) default else value
}

.dc_observed_array <- function(observed) {
  if (inherits(observed, "Tensor")) {
    out <- array(as.numeric(observed@data), dim = dim(observed))
  } else {
    out <- observed
  }
  if (!is.array(out) || length(dim(out)) < 2L || any(dim(out) < 1L) ||
      !is.numeric(out) || any(!is.finite(out))) {
    stop("`observed` must be a finite numeric tensor of order at least two.",
         call. = FALSE)
  }
  out
}

.dc_validate_components <- function(components, dims, name = "initialization") {
  if (!is.list(components) || !length(components)) {
    stop("`", name, "` must be a nonempty component list.", call. = FALSE)
  }
  valid <- vapply(components, function(component) {
    is.list(component) && length(component$lambda) == 1L &&
      is.numeric(component$lambda) && is.finite(component$lambda) &&
      is.list(component$U) && length(component$U) == length(dims) &&
      all(vapply(seq_along(dims), function(mode) {
        u <- component$U[[mode]]
        is.numeric(u) && length(u) == dims[[mode]] && all(is.finite(u))
      }, logical(1L)))
  }, logical(1L))
  if (!all(valid)) {
    stop("`", name, "` is not compatible with the observed dimensions.",
         call. = FALSE)
  }
  invisible(components)
}

.dc_problem <- function(observed) {
  if (!requireNamespace("rTensor", quietly = TRUE)) {
    stop("The decomposition comparators require the 'rTensor' package.",
         call. = FALSE)
  }
  observed <- .dc_observed_array(observed)
  tensor <- rTensor::as.tensor(observed)
  dims <- dim(observed)
  structure(list(
    observed = observed,
    tensor = tensor,
    dims = dims,
    order = length(dims),
    y_norm_sq = sum(observed * observed),
    unfoldings = lapply(seq_along(dims), function(mode) {
      rTensor::k_unfold(tensor, m = mode)@data
    })
  ), class = c("decomposition_comparator_problem", "list"))
}

.dc_validate_problem <- function(problem) {
  required <- c("observed", "tensor", "dims", "order", "y_norm_sq",
                "unfoldings")
  if (!inherits(problem, "decomposition_comparator_problem") ||
      !all(required %in% names(problem)) ||
      length(problem$unfoldings) != problem$order ||
      !identical(as.integer(dim(problem$observed)),
                 as.integer(problem$dims))) {
    stop("The prepared decomposition problem is invalid.", call. = FALSE)
  }
  problem
}

.dc_raw_factors <- function(components) {
  modes <- length(components[[1L]]$U)
  list(
    factors = lapply(seq_len(modes), function(mode) {
      do.call(cbind, lapply(components, function(z) as.numeric(z$U[[mode]])))
    }),
    weights = vapply(components, `[[`, numeric(1L), "lambda")
  )
}

.dc_balanced_factors <- function(components) {
  raw <- .dc_raw_factors(components)
  modes <- length(raw$factors)
  scales <- abs(raw$weights)^(1 / modes)
  factors <- lapply(raw$factors, function(A) sweep(A, 2L, scales, "*"))
  factors[[1L]] <- sweep(
    factors[[1L]], 2L, ifelse(raw$weights < 0, -1, 1), "*"
  )
  factors
}

.dc_components_from_factors <- function(factors, floor = 1e-14) {
  rank <- ncol(factors[[1L]])
  lapply(seq_len(rank), function(component) {
    norms <- vapply(factors, function(A) {
      sqrt(sum(A[, component] * A[, component]))
    }, numeric(1L))
    if (any(!is.finite(norms)) || any(norms <= floor)) {
      return(list(
        lambda = 0,
        U = lapply(factors, function(A) {
          u <- A[, component]
          norm <- sqrt(sum(u * u))
          if (!is.finite(norm) || norm <= floor) {
            u <- numeric(nrow(A)); u[[1L]] <- 1
            u
          } else u / norm
        })
      ))
    }
    list(
      lambda = prod(norms),
      U = lapply(seq_along(factors), function(mode) {
        factors[[mode]][, component] / norms[[mode]]
      })
    )
  })
}

.dc_components_from_weighted_factors <- function(factors, weights,
                                                  floor = 1e-14) {
  rank <- ncol(factors[[1L]])
  lapply(seq_len(rank), function(component) {
    norms <- vapply(factors, function(A) {
      sqrt(sum(A[, component] * A[, component]))
    }, numeric(1L))
    if (any(!is.finite(norms)) || any(norms <= floor) ||
        !is.finite(weights[[component]])) {
      return(list(
        lambda = 0,
        U = lapply(factors, function(A) {
          u <- A[, component]
          norm <- sqrt(sum(u * u))
          if (!is.finite(norm) || norm <= floor) {
            u <- numeric(nrow(A)); u[[1L]] <- 1
            u
          } else u / norm
        })
      ))
    }
    list(
      lambda = as.numeric(weights[[component]] * prod(norms)),
      U = lapply(seq_along(factors), function(mode) {
        factors[[mode]][, component] / norms[[mode]]
      })
    )
  })
}

.dc_reconstruct_components <- function(components) {
  dims <- vapply(components[[1L]]$U, length, integer(1L))
  reconstructed <- array(0, dim = dims)
  for (component in components) {
    value <- component$lambda * as.numeric(component$U[[1L]])
    if (length(dims) >= 2L) {
      for (mode in 2:length(dims)) {
        value <- outer(value, as.numeric(component$U[[mode]]))
      }
    }
    dim(value) <- dims
    reconstructed <- reconstructed + value
  }
  reconstructed
}

.dc_relative_error <- function(estimate, truth) {
  estimate_dense <- .dc_reconstruct_components(estimate)
  truth_dense <- .dc_reconstruct_components(truth)
  truth_norm_sq <- sum(truth_dense * truth_dense)
  if (!is.finite(truth_norm_sq) || truth_norm_sq <= 0) {
    stop("`truth` has zero or nonfinite Frobenius norm.", call. = FALSE)
  }
  sqrt(sum((estimate_dense - truth_dense)^2) / truth_norm_sq)
}

.dc_mttkrp <- function(problem, factors, mode) {
  khatri_rao <- rTensor::khatri_rao_list(factors[-mode], reverse = TRUE)
  problem$unfoldings[[mode]] %*% khatri_rao
}

.dc_objective <- function(problem, factors) {
  gram <- Reduce("*", lapply(factors, crossprod))
  inner <- sum(factors[[1L]] * .dc_mttkrp(problem, factors, 1L))
  value <- 0.5 * (problem$y_norm_sq - 2 * inner + sum(gram))
  if (is.finite(value) && value < 0 &&
      value > -1e-8 * max(problem$y_norm_sq, 1)) 0 else value
}

.dc_objective_components <- function(problem, components) {
  residual <- .dc_reconstruct_components(components) - problem$observed
  0.5 * sum(residual * residual)
}

.dc_solver_state <- function() {
  state <- new.env(parent = emptyenv())
  state$pseudoinverse_fallbacks <- 0L
  state
}

.dc_right_solve <- function(B, G, diagonal_shift = 0, state = NULL) {
  rank <- ncol(G)
  shifted <- (G + t(G)) / 2 + diag(diagonal_shift, rank)
  answer <- tryCatch(t(solve(shifted, t(B))), error = function(e) NULL)
  if (!is.null(answer) && all(is.finite(answer))) return(answer)
  if (!is.null(state)) {
    state$pseudoinverse_fallbacks <- state$pseudoinverse_fallbacks + 1L
  }
  eigen_decomposition <- eigen(shifted, symmetric = TRUE)
  tolerance <- max(dim(shifted)) * .Machine$double.eps *
    max(abs(eigen_decomposition$values), 1)
  inverse <- ifelse(abs(eigen_decomposition$values) > tolerance,
                    1 / eigen_decomposition$values, 0)
  B %*% (eigen_decomposition$vectors %*%
           (inverse * t(eigen_decomposition$vectors)))
}

.dc_record <- function(iteration, components, problem, truth, diagnostics,
                       objective_rows, error_rows) {
  if (!diagnostics) {
    return(list(objective = objective_rows, error = error_rows))
  }
  objective_rows[[length(objective_rows) + 1L]] <- data.frame(
    iteration = as.integer(iteration),
    objective = .dc_objective_components(problem, components),
    stringsAsFactors = FALSE
  )
  if (!is.null(truth)) {
    error_rows[[length(error_rows) + 1L]] <- data.frame(
      iteration = as.integer(iteration),
      relative_frobenius_error = .dc_relative_error(components, truth),
      stringsAsFactors = FALSE
    )
  }
  list(objective = objective_rows, error = error_rows)
}

.dc_orthogonalize <- function(factors) {
  rank <- ncol(factors[[1L]])
  for (mode in seq_along(factors)) {
    Q <- factors[[mode]]
    for (component in seq_len(rank)) {
      norm <- sqrt(sum(Q[, component] * Q[, component]))
      if (!is.finite(norm) || norm <= .Machine$double.eps) {
        stop("Orth-ALS encountered a zero Gram-Schmidt column.",
             call. = FALSE)
      }
      Q[, component] <- Q[, component] / norm
      if (component < rank) {
        for (later in seq.int(component + 1L, rank)) {
          Q[, later] <- Q[, later] -
            sum(Q[, component] * Q[, later]) * Q[, component]
        }
      }
    }
    factors[[mode]] <- Q
  }
  factors
}

.dc_fit_orth_als <- function(problem, initialization, iterations, truth,
                             diagnostics, options, state) {
  raw <- .dc_raw_factors(initialization)
  factors <- raw$factors
  weights <- raw$weights
  objective_rows <- error_rows <- list()
  recorded <- .dc_record(0L, initialization, problem, truth, diagnostics,
                         objective_rows, error_rows)
  objective_rows <- recorded$objective; error_rows <- recorded$error
  for (iteration in seq_len(iterations)) {
    factors <- .dc_orthogonalize(factors)
    grams <- lapply(factors, crossprod)
    for (mode in seq_along(factors)) {
      normal <- Reduce("*", grams[-mode])
      updated <- .dc_right_solve(
        .dc_mttkrp(problem, factors, mode), normal,
        diagonal_shift = 0, state = state
      )
      norms <- sqrt(colSums(updated * updated))
      if (any(!is.finite(norms)) || any(norms <= .Machine$double.eps)) {
        stop("Orth-ALS encountered a zero mode-update column.",
             call. = FALSE)
      }
      factors[[mode]] <- sweep(updated, 2L, norms, "/")
    }
    contractions <- colSums(
      factors[[1L]] * .dc_mttkrp(problem, factors, 1L)
    )
    weights <- as.numeric(contractions)
    components <- .dc_components_from_weighted_factors(factors, weights)
    recorded <- .dc_record(
      iteration, components, problem, truth, diagnostics,
      objective_rows, error_rows
    )
    objective_rows <- recorded$objective; error_rows <- recorded$error
  }
  list(
    components = .dc_components_from_weighted_factors(factors, weights),
    objective_rows = objective_rows, error_rows = error_rows,
    accepted = iterations,
    extras = list(orthogonalization = "classical_Gram_Schmidt_before_each_sweep")
  )
}

.dc_right_inverse <- function(A) {
  decomposition <- svd(A)
  tolerance <- max(dim(A)) * .Machine$double.eps * decomposition$d[[1L]]
  inverse <- ifelse(decomposition$d > tolerance, 1 / decomposition$d, 0)
  decomposition$v %*% (inverse * t(decomposition$u))
}

.dc_fit_ico <- function(problem, initialization, iterations, truth,
                        diagnostics, options, state) {
  raw <- .dc_raw_factors(initialization)
  factors <- raw$factors
  weights <- raw$weights
  rank <- length(initialization)
  objective_rows <- error_rows <- list()
  recorded <- .dc_record(0L, initialization, problem, truth, diagnostics,
                         objective_rows, error_rows)
  objective_rows <- recorded$objective; error_rows <- recorded$error
  zero_contractions <- 0L
  for (iteration in seq_len(iterations)) {
    for (mode in seq_along(factors)) {
      inverses <- lapply(factors, .dc_right_inverse)
      dual_factors <- lapply(inverses, t)
      updated <- .dc_mttkrp(problem, dual_factors, mode)
      norms <- sqrt(colSums(updated * updated))
      zero <- !is.finite(norms) | norms <= .Machine$double.eps
      zero_contractions <- zero_contractions + sum(zero)
      if (any(zero)) updated[, zero] <- 0
      if (any(!zero)) {
        updated[, !zero] <- sweep(updated[, !zero, drop = FALSE], 2L,
                                  norms[!zero], "/")
      }
      factors[[mode]] <- updated
    }
    inverses <- lapply(factors, .dc_right_inverse)
    dual_factors <- lapply(inverses, t)
    weights <- abs(colSums(
      dual_factors[[1L]] * .dc_mttkrp(problem, dual_factors, 1L)
    ))
    components <- .dc_components_from_weighted_factors(factors, weights)
    recorded <- .dc_record(
      iteration, components, problem, truth, diagnostics,
      objective_rows, error_rows
    )
    objective_rows <- recorded$objective; error_rows <- recorded$error
  }
  list(
    components = .dc_components_from_weighted_factors(factors, weights),
    objective_rows = objective_rows, error_rows = error_rows,
    accepted = iterations,
    extras = list(
      pseudoinverse = "SVD_Moore_Penrose_MATLAB_pinv_tolerance",
      zero_contractions = zero_contractions
    )
  )
}

.dc_balance_factors <- function(factors, floor = 1e-14) {
  rank <- ncol(factors[[1L]])
  for (component in seq_len(rank)) {
    norms <- vapply(factors, function(A) {
      sqrt(sum(A[, component] * A[, component]))
    }, numeric(1L))
    if (any(!is.finite(norms)) || any(norms <= floor)) next
    common <- exp(mean(log(norms)))
    for (mode in seq_along(factors)) {
      factors[[mode]][, component] <- factors[[mode]][, component] *
        common / norms[[mode]]
    }
  }
  factors
}

.dc_normalize_for_mode <- function(factors, mode, floor = 1e-14) {
  rank <- ncol(factors[[1L]])
  other_modes <- setdiff(seq_along(factors), mode)
  for (component in seq_len(rank)) {
    multiplier <- 1
    for (other in other_modes) {
      norm <- sqrt(sum(factors[[other]][, component]^2))
      if (is.finite(norm) && norm > floor) {
        factors[[other]][, component] <-
          factors[[other]][, component] / norm
        multiplier <- multiplier * norm
      }
    }
    factors[[mode]][, component] <-
      factors[[mode]][, component] * multiplier
  }
  factors
}

.dc_gram_products <- function(factors) {
  grams <- lapply(factors, crossprod)
  list(
    grams = grams,
    leave_one = lapply(seq_along(factors), function(mode) {
      Reduce("*", grams[-mode])
    })
  )
}

.dc_gradient <- function(problem, factors) {
  products <- .dc_gram_products(factors)
  gradient <- lapply(seq_along(factors), function(mode) {
    factors[[mode]] %*% products$leave_one[[mode]] -
      .dc_mttkrp(problem, factors, mode)
  })
  list(gradient = gradient, leave_one = products$leave_one)
}

.dc_pack <- function(factors) {
  unlist(lapply(factors, as.numeric), use.names = FALSE)
}

.dc_unpack <- function(vector, dims, rank) {
  factors <- vector("list", length(dims))
  offset <- 0L
  for (mode in seq_along(dims)) {
    count <- dims[[mode]] * rank
    factors[[mode]] <- matrix(
      vector[offset + seq_len(count)], nrow = dims[[mode]], ncol = rank
    )
    offset <- offset + count
  }
  factors
}

.dc_jtj_apply <- function(factors, direction, damping = 0) {
  modes <- length(factors)
  rank <- ncol(factors[[1L]])
  grams <- lapply(factors, crossprod)
  out <- lapply(seq_len(modes), function(mode) {
    direction[[mode]] %*% Reduce("*", grams[-mode])
  })
  for (mode in seq_len(modes)) {
    for (other in setdiff(seq_len(modes), mode)) {
      remaining <- setdiff(seq_len(modes), c(mode, other))
      product <- if (length(remaining)) {
        Reduce("*", grams[remaining])
      } else matrix(1, nrow = rank, ncol = rank)
      cross <- crossprod(factors[[other]], direction[[other]])
      out[[mode]] <- out[[mode]] + factors[[mode]] %*% t(product * cross)
    }
    if (damping > 0) {
      out[[mode]] <- out[[mode]] + damping * direction[[mode]]
    }
  }
  out
}

.dc_block_precondition <- function(factors, residual, damping, state) {
  grams <- lapply(factors, crossprod)
  lapply(seq_along(factors), function(mode) {
    normal <- Reduce("*", grams[-mode])
    .dc_right_solve(residual[[mode]], normal, damping, state)
  })
}

.dc_pcg <- function(apply_operator, right_hand_side, apply_preconditioner,
                    maxit, tolerance) {
  solution <- numeric(length(right_hand_side))
  residual <- right_hand_side
  preconditioned <- apply_preconditioner(residual)
  direction <- preconditioned
  residual_inner <- sum(residual * preconditioned)
  scale <- max(sqrt(sum(right_hand_side^2)), .Machine$double.eps)
  if (!is.finite(residual_inner) || residual_inner <= 0) {
    return(list(solution = solution, iterations = 0L, converged = FALSE))
  }
  completed <- 0L
  for (iteration in seq_len(maxit)) {
    completed <- iteration
    operator_direction <- apply_operator(direction)
    denominator <- sum(direction * operator_direction)
    if (!is.finite(denominator) || denominator <= 0) break
    step <- residual_inner / denominator
    solution <- solution + step * direction
    residual <- residual - step * operator_direction
    if (sqrt(sum(residual^2)) <= tolerance * scale) {
      return(list(solution = solution, iterations = iteration,
                  converged = TRUE))
    }
    preconditioned <- apply_preconditioner(residual)
    next_inner <- sum(residual * preconditioned)
    if (!is.finite(next_inner) || next_inner <= 0) break
    direction <- preconditioned + (next_inner / residual_inner) * direction
    residual_inner <- next_inner
  }
  list(solution = solution, iterations = completed, converged = FALSE)
}

.dc_boundary_tau <- function(point, direction, radius) {
  quadratic <- sum(direction * direction)
  linear <- 2 * sum(point * direction)
  constant <- sum(point * point) - radius^2
  discriminant <- max(linear^2 - 4 * quadratic * constant, 0)
  roots <- c(
    (-linear + sqrt(discriminant)) / (2 * quadratic),
    (-linear - sqrt(discriminant)) / (2 * quadratic)
  )
  roots <- roots[is.finite(roots) & roots >= 0]
  if (!length(roots)) 0 else min(roots)
}

.dc_steihaug <- function(apply_operator, right_hand_side,
                         apply_preconditioner, radius, maxit, tolerance) {
  solution <- numeric(length(right_hand_side))
  residual <- right_hand_side
  preconditioned <- apply_preconditioner(residual)
  direction <- preconditioned
  residual_inner <- sum(residual * preconditioned)
  scale <- max(sqrt(sum(right_hand_side^2)), .Machine$double.eps)
  if (!is.finite(residual_inner) || residual_inner <= 0) {
    return(list(solution = solution, iterations = 0L, boundary = FALSE,
                converged = FALSE))
  }
  completed <- 0L
  for (iteration in seq_len(maxit)) {
    completed <- iteration
    operator_direction <- apply_operator(direction)
    curvature <- sum(direction * operator_direction)
    if (!is.finite(curvature) || curvature <= 0) {
      tau <- .dc_boundary_tau(solution, direction, radius)
      return(list(solution = solution + tau * direction,
                  iterations = iteration, boundary = TRUE, converged = FALSE))
    }
    step <- residual_inner / curvature
    trial <- solution + step * direction
    if (sqrt(sum(trial^2)) >= radius) {
      tau <- .dc_boundary_tau(solution, direction, radius)
      return(list(solution = solution + tau * direction,
                  iterations = iteration, boundary = TRUE, converged = FALSE))
    }
    solution <- trial
    residual <- residual - step * operator_direction
    if (sqrt(sum(residual^2)) <= tolerance * scale) {
      return(list(solution = solution, iterations = iteration,
                  boundary = FALSE, converged = TRUE))
    }
    preconditioned <- apply_preconditioner(residual)
    next_inner <- sum(residual * preconditioned)
    if (!is.finite(next_inner) || next_inner <= 0) break
    direction <- preconditioned + (next_inner / residual_inner) * direction
    residual_inner <- next_inner
  }
  list(solution = solution, iterations = completed, boundary = FALSE,
       converged = FALSE)
}

.dc_generic_history <- function(initialization, problem, truth, diagnostics) {
  objective_rows <- error_rows <- list()
  recorded <- .dc_record(0L, initialization, problem, truth, diagnostics,
                         objective_rows, error_rows)
  list(objective_rows = recorded$objective, error_rows = recorded$error)
}

.dc_append_generic_history <- function(history, iteration, factors, problem,
                                       truth, diagnostics) {
  components <- .dc_components_from_factors(factors)
  recorded <- .dc_record(
    iteration, components, problem, truth, diagnostics,
    history$objective_rows, history$error_rows
  )
  history$objective_rows <- recorded$objective
  history$error_rows <- recorded$error
  history
}

.dc_fit_nls_gn <- function(problem, initialization, iterations, truth,
                           diagnostics, options, state) {
  factors <- .dc_balanced_factors(initialization)
  dims <- problem$dims
  rank <- length(initialization)
  cg_maxit <- .dc_whole_number(.dc_option(options, "cg_maxit", 20L),
                               "options$cg_maxit")
  cg_tolerance <- .dc_scalar(.dc_option(options, "cg_tolerance", 1e-3),
                             "options$cg_tolerance", 0, strict_lower = TRUE)
  damping_factor <- .dc_scalar(
    .dc_option(options, "nls_damping_factor", 1e-10),
    "options$nls_damping_factor", 0
  )
  norm <- sqrt(sum(.dc_pack(factors)^2))
  radius <- .dc_option(options, "trust_radius", NULL)
  if (is.null(radius)) radius <- max(0.25 * norm, 1)
  radius <- .dc_scalar(radius, "options$trust_radius", 0, strict_lower = TRUE)
  maximum_radius <- .dc_option(options, "trust_max", NULL)
  if (is.null(maximum_radius)) maximum_radius <- max(4 * norm, radius)
  maximum_radius <- .dc_scalar(maximum_radius, "options$trust_max", radius)
  objective <- .dc_objective(problem, factors)
  history <- .dc_generic_history(initialization, problem, truth, diagnostics)
  inner_iterations <- integer(iterations)
  accepted <- logical(iterations)
  damping_path <- numeric(iterations)
  for (iteration in seq_len(iterations)) {
    gradient <- .dc_gradient(problem, factors)
    gradient_vector <- .dc_pack(gradient$gradient)
    damping <- damping_factor * max(
      1,
      max(vapply(gradient$leave_one, function(normal) max(diag(normal)),
                 numeric(1L)))
    )
    damping_path[[iteration]] <- damping
    apply_operator <- function(vector) {
      direction <- .dc_unpack(vector, dims, rank)
      .dc_pack(.dc_jtj_apply(factors, direction, damping))
    }
    apply_preconditioner <- function(vector) {
      residual <- .dc_unpack(vector, dims, rank)
      .dc_pack(.dc_block_precondition(
        factors, residual, damping = damping, state = state
      ))
    }
    step <- .dc_steihaug(
      apply_operator, -gradient_vector, apply_preconditioner,
      radius = radius, maxit = cg_maxit, tolerance = cg_tolerance
    )
    inner_iterations[[iteration]] <- step$iterations
    operator_step <- apply_operator(step$solution)
    predicted <- -sum(gradient_vector * step$solution) -
      0.5 * sum(step$solution * operator_step)
    trial_factors <- .dc_unpack(
      .dc_pack(factors) + step$solution, dims, rank
    )
    trial_objective <- .dc_objective(problem, trial_factors)
    actual <- objective - trial_objective
    ratio <- if (is.finite(predicted) && predicted > 0) {
      actual / predicted
    } else -Inf
    if (!is.finite(ratio) || ratio < 0.25) radius <- max(0.25 * radius, 1e-10)
    if (is.finite(ratio) && ratio > 0.75 && step$boundary) {
      radius <- min(2 * radius, maximum_radius)
    }
    if (is.finite(ratio) && ratio > 1e-4 && is.finite(trial_objective)) {
      factors <- .dc_balance_factors(trial_factors)
      objective <- .dc_objective(problem, factors)
      accepted[[iteration]] <- TRUE
    }
    history <- .dc_append_generic_history(
      history, iteration, factors, problem, truth, diagnostics
    )
  }
  list(
    components = .dc_components_from_factors(factors),
    objective_rows = history$objective_rows, error_rows = history$error_rows,
    accepted = sum(accepted),
    extras = list(
      accepted_path = accepted, inner_iterations = inner_iterations,
      inner_converged_count = NA_integer_, trust_radius_final = radius,
      trust_radius_max = maximum_radius, damping_path = damping_path,
      damping_role = "intrinsic_damped_Gauss_Newton_normal_operator"
    )
  )
}

.dc_fit_gn_cg <- function(problem, initialization, iterations, truth,
                          diagnostics, options, state) {
  factors <- .dc_balanced_factors(initialization)
  dims <- problem$dims
  rank <- length(initialization)
  cg_maxit <- .dc_whole_number(.dc_option(options, "cg_maxit", 20L),
                               "options$cg_maxit")
  cg_tolerance <- .dc_scalar(.dc_option(options, "cg_tolerance", 1e-3),
                             "options$cg_tolerance", 0, strict_lower = TRUE)
  armijo <- .dc_scalar(.dc_option(options, "armijo", 1e-4),
                       "options$armijo", 0)
  max_backtracks <- .dc_whole_number(
    .dc_option(options, "max_backtracks", 12L),
    "options$max_backtracks", positive = FALSE
  )
  initial_damping <- .dc_option(options, "gn_damping", NULL)
  if (!is.null(initial_damping)) {
    initial_damping <- .dc_scalar(initial_damping, "options$gn_damping", 0,
                                  strict_lower = TRUE)
  }
  objective <- .dc_objective(problem, factors)
  history <- .dc_generic_history(initialization, problem, truth, diagnostics)
  inner_iterations <- integer(iterations)
  inner_converged <- logical(iterations)
  accepted <- logical(iterations)
  step_lengths <- numeric(iterations)
  damping_path <- numeric(iterations)
  damping <- initial_damping
  for (iteration in seq_len(iterations)) {
    gradient <- .dc_gradient(problem, factors)
    if (is.null(damping)) {
      damping <- 1e-4 * max(
        1,
        max(vapply(gradient$leave_one, function(normal) max(diag(normal)),
                   numeric(1L)))
      )
    }
    damping_path[[iteration]] <- damping
    gradient_vector <- .dc_pack(gradient$gradient)
    apply_operator <- function(vector) {
      direction <- .dc_unpack(vector, dims, rank)
      .dc_pack(.dc_jtj_apply(factors, direction, damping))
    }
    apply_preconditioner <- function(vector) {
      residual <- .dc_unpack(vector, dims, rank)
      .dc_pack(.dc_block_precondition(
        factors, residual, damping = damping, state = state
      ))
    }
    step <- .dc_pcg(
      apply_operator, -gradient_vector, apply_preconditioner,
      maxit = cg_maxit, tolerance = cg_tolerance
    )
    inner_iterations[[iteration]] <- step$iterations
    inner_converged[[iteration]] <- step$converged
    direction <- step$solution
    directional_derivative <- sum(gradient_vector * direction)
    if (!is.finite(directional_derivative) || directional_derivative >= 0) {
      direction <- -apply_preconditioner(gradient_vector)
      directional_derivative <- sum(gradient_vector * direction)
    }
    step_length <- 1
    trial_objective <- Inf
    trial_factors <- factors
    for (backtrack in 0:max_backtracks) {
      trial_factors <- .dc_unpack(
        .dc_pack(factors) + step_length * direction, dims, rank
      )
      trial_objective <- .dc_objective(problem, trial_factors)
      if (is.finite(trial_objective) &&
          trial_objective <= objective +
            armijo * step_length * directional_derivative) break
      step_length <- step_length / 2
    }
    if (is.finite(trial_objective) && trial_objective < objective) {
      factors <- .dc_balance_factors(trial_factors)
      objective <- .dc_objective(problem, factors)
      damping <- max(damping / 3, 1e-12)
      accepted[[iteration]] <- TRUE
      step_lengths[[iteration]] <- step_length
    } else {
      damping <- min(damping * 10, 1e12)
      step_lengths[[iteration]] <- 0
    }
    history <- .dc_append_generic_history(
      history, iteration, factors, problem, truth, diagnostics
    )
  }
  list(
    components = .dc_components_from_factors(factors),
    objective_rows = history$objective_rows, error_rows = history$error_rows,
    accepted = sum(accepted),
    extras = list(
      accepted_path = accepted, inner_iterations = inner_iterations,
      inner_converged = inner_converged, step_lengths = step_lengths,
      damping_path = damping_path, damping_final = damping,
      damping_role = "intrinsic_adaptive_Gauss_Newton_damping"
    )
  )
}

.dc_als_sweep <- function(problem, factors, state, stabilizer = 0) {
  for (mode in seq_along(factors)) {
    normal <- Reduce("*", lapply(factors[-mode], crossprod))
    scale <- max(1, max(abs(diag(normal))))
    factors[[mode]] <- .dc_right_solve(
      .dc_mttkrp(problem, factors, mode), normal,
      diagonal_shift = stabilizer * scale, state = state
    )
  }
  .dc_balance_factors(factors)
}

.dc_polynomial_product <- function(left, right) {
  out <- numeric(length(left) + length(right) - 1L)
  for (i in seq_along(left)) {
    for (j in seq_along(right)) {
      out[[i + j - 1L]] <- out[[i + j - 1L]] + left[[i]] * right[[j]]
    }
  }
  out
}

.dc_rank_one_inner <- function(problem, vectors) {
  matrices <- lapply(vectors, function(vector) matrix(vector, nrow = 1L))
  as.numeric(rTensor::ttl(
    problem$tensor, matrices, ms = seq_len(problem$order)
  )@data)
}

# Exact coefficients of ||Y - CP(A + rho D)||_F^2. This implements the same
# common-relaxation polynomial as the legacy Rcpp entry point without compiling
# code while the adapter is sourced.
.dc_els_polynomial <- function(problem, base, direction) {
  order <- problem$order
  rank <- ncol(base[[1L]])
  cross_coefficients <- numeric(order + 1L)
  for (component in seq_len(rank)) {
    for (mask in 0:(2^order - 1L)) {
      bits <- as.integer(intToBits(mask))[seq_len(order)]
      degree <- sum(bits)
      vectors <- lapply(seq_len(order), function(mode) {
        if (bits[[mode]]) direction[[mode]][, component] else
          base[[mode]][, component]
      })
      cross_coefficients[[degree + 1L]] <-
        cross_coefficients[[degree + 1L]] +
        .dc_rank_one_inner(problem, vectors)
    }
  }
  norm_coefficients <- numeric(2L * order + 1L)
  for (left in seq_len(rank)) {
    for (right in seq_len(rank)) {
      pair_polynomial <- 1
      for (mode in seq_len(order)) {
        a_left <- base[[mode]][, left]
        a_right <- base[[mode]][, right]
        d_left <- direction[[mode]][, left]
        d_right <- direction[[mode]][, right]
        mode_polynomial <- c(
          sum(a_left * a_right),
          sum(d_left * a_right) + sum(a_left * d_right),
          sum(d_left * d_right)
        )
        pair_polynomial <- .dc_polynomial_product(
          pair_polynomial, mode_polynomial
        )
      }
      norm_coefficients[seq_along(pair_polynomial)] <-
        norm_coefficients[seq_along(pair_polynomial)] + pair_polynomial
    }
  }
  cross_padded <- c(cross_coefficients,
                    numeric(length(norm_coefficients) -
                              length(cross_coefficients)))
  polynomial <- norm_coefficients - 2 * cross_padded
  polynomial[[1L]] <- polynomial[[1L]] + problem$y_norm_sq
  polynomial
}

.dc_polynomial_value <- function(coefficients, point) {
  sum(coefficients * point^(seq_along(coefficients) - 1L))
}

.dc_els_relaxation <- function(problem, base, direction, max_abs) {
  polynomial <- .dc_els_polynomial(problem, base, direction)
  derivative <- (seq_len(length(polynomial) - 1L)) * polynomial[-1L]
  tolerance <- .Machine$double.eps * max(1, max(abs(derivative)))
  while (length(derivative) > 1L &&
         abs(utils::tail(derivative, 1L)) <= tolerance) {
    derivative <- utils::head(derivative, -1L)
  }
  roots <- if (length(derivative) <= 1L) complex() else polyroot(derivative)
  real_roots <- Re(roots)[
    abs(Im(roots)) <= 1e-7 * (1 + abs(Re(roots)))
  ]
  candidates <- unique(c(
    0, 1,
    real_roots[is.finite(real_roots) & abs(real_roots) <= max_abs]
  ))
  values <- vapply(candidates, function(point) {
    .dc_polynomial_value(polynomial, point)
  }, numeric(1L))
  list(relaxation = candidates[[which.min(values)]],
       polynomial = polynomial)
}

.dc_fit_als_els <- function(problem, initialization, iterations, truth,
                            diagnostics, options, state) {
  stabilizer <- .dc_scalar(
    .dc_option(options, "linear_solve_stabilizer", 0),
    "options$linear_solve_stabilizer", 0
  )
  maximum_relaxation <- .dc_scalar(
    .dc_option(options, "els_max_abs", 100), "options$els_max_abs", 0,
    strict_lower = TRUE
  )
  base <- .dc_balanced_factors(initialization)
  latest <- base
  history <- .dc_generic_history(initialization, problem, truth, diagnostics)
  relaxations <- rep(NA_real_, iterations)
  polynomial_degrees <- integer(iterations)
  for (iteration in seq_len(iterations)) {
    if (iteration == 1L) {
      next_factors <- .dc_als_sweep(problem, latest, state, stabilizer)
      relaxations[[iteration]] <- 1
    } else {
      direction <- Map("-", latest, base)
      relaxation <- .dc_els_relaxation(
        problem, base, direction, maximum_relaxation
      )
      polynomial_degrees[[iteration]] <- length(relaxation$polynomial) - 1L
      extrapolated <- Map(
        function(A, D) A + relaxation$relaxation * D, base, direction
      )
      rho <- relaxation$relaxation
      if (!is.finite(.dc_objective(problem, extrapolated)) ||
          .dc_objective(problem, extrapolated) >
            .dc_objective(problem, latest)) {
        rho <- 1
        extrapolated <- latest
      }
      next_factors <- .dc_als_sweep(
        problem, extrapolated, state, stabilizer
      )
      base <- extrapolated
      relaxations[[iteration]] <- rho
    }
    latest <- next_factors
    history <- .dc_append_generic_history(
      history, iteration, latest, problem, truth, diagnostics
    )
  }
  list(
    components = .dc_components_from_factors(latest),
    objective_rows = history$objective_rows, error_rows = history$error_rows,
    accepted = iterations,
    extras = list(
      relaxations = relaxations,
      polynomial_degrees = polynomial_degrees,
      linear_solve_stabilizer = stabilizer,
      stabilizer_role = if (stabilizer == 0) "none" else
        "explicit_numerical_diagonal_shift_for_legacy_replay"
    )
  )
}

.dc_amdm_mode_update <- function(problem, factors, mode, threshold,
                                 sv_floor, stabilizer, state) {
  factors <- .dc_normalize_for_mode(factors, mode)
  other_modes <- setdiff(seq_along(factors), mode)
  transformed_factors <- vector("list", length(other_modes))
  transformed_grams <- vector("list", length(other_modes))
  for (index in seq_along(other_modes)) {
    other <- other_modes[[index]]
    decomposition <- svd(
      factors[[other]],
      nu = ncol(factors[[other]]), nv = ncol(factors[[other]])
    )
    values <- pmax(
      decomposition$d,
      sv_floor * max(decomposition$d, 1)
    )
    transformed <- values
    if (threshold > 0L) {
      take <- seq_len(min(threshold, length(values)))
      transformed[take] <- 1 / values[take]
    }
    transformed_factors[[index]] <- decomposition$u %*%
      (transformed * t(decomposition$v))
    transformed_grams[[index]] <- decomposition$v %*%
      ((transformed * values) * t(decomposition$v))
  }
  khatri_rao <- rTensor::khatri_rao_list(
    transformed_factors, reverse = TRUE
  )
  right_hand_side <- problem$unfoldings[[mode]] %*% khatri_rao
  normal <- Reduce("*", transformed_grams)
  scale <- max(1, max(abs(diag(normal))))
  factors[[mode]] <- .dc_right_solve(
    right_hand_side, normal,
    diagonal_shift = stabilizer * scale, state = state
  )
  factors
}

.dc_fit_hybrid_amdm <- function(problem, initialization, iterations, truth,
                                diagnostics, options, state) {
  hybrid_every <- .dc_whole_number(
    .dc_option(options, "hybrid_every", 10L), "options$hybrid_every"
  )
  sv_floor <- .dc_scalar(.dc_option(options, "sv_floor", 1e-10),
                         "options$sv_floor", 0)
  stabilizer <- .dc_scalar(
    .dc_option(options, "linear_solve_stabilizer", 0),
    "options$linear_solve_stabilizer", 0
  )
  rank <- length(initialization)
  factors <- .dc_balanced_factors(initialization)
  history <- .dc_generic_history(initialization, problem, truth, diagnostics)
  thresholds <- integer(iterations)
  for (iteration in seq_len(iterations)) {
    threshold <- max(rank - floor((iteration - 1L) / hybrid_every), 0L)
    thresholds[[iteration]] <- threshold
    for (mode in seq_along(factors)) {
      factors <- .dc_amdm_mode_update(
        problem, factors, mode, threshold, sv_floor, stabilizer, state
      )
    }
    factors <- .dc_balance_factors(factors)
    history <- .dc_append_generic_history(
      history, iteration, factors, problem, truth, diagnostics
    )
  }
  list(
    components = .dc_components_from_factors(factors),
    objective_rows = history$objective_rows, error_rows = history$error_rows,
    accepted = iterations,
    extras = list(
      thresholds = thresholds, hybrid_every = hybrid_every,
      hybrid_transition_count = sum(diff(thresholds) != 0),
      sv_floor = sv_floor,
      sv_floor_role = "numerical_SVD_inverse_threshold",
      linear_solve_stabilizer = stabilizer,
      stabilizer_role = if (stabilizer == 0) "none" else
        "explicit_numerical_diagonal_shift_for_legacy_replay"
    )
  )
}

#' Fit a repository-only CP-decomposition comparator
#'
#' @param method One of Orth-ALS, ICO, NLS-GN, GN-CG, CP-ALS+ELS, or
#'   Hybrid AMDM.
#' @param observed Numeric observed tensor or `rTensor::Tensor`.
#' @param initialization Realized common CP component initialization.
#' @param iterations Number of complete outer updates.
#' @param truth Optional CP truth used only for oracle diagnostics.
#' @param diagnostics Whether to retain iteration-wise objective/error records.
#' @param options Method controls. A statistical ridge is not permitted.
#' @return Common comparator fit object with final components, trajectory,
#'   update counts, stopping information, method diagnostics, and provenance.
fit_decomposition_comparator <- function(
    method,
    observed,
    initialization,
    iterations,
    truth = NULL,
    diagnostics = TRUE,
    options = list()) {
  methods <- c("Orth-ALS", "ICO", "NLS-GN", "GN-CG",
               "CP-ALS+ELS", "Hybrid AMDM")
  method <- match.arg(method, methods)
  iterations <- .dc_whole_number(iterations, "iterations")
  if (length(diagnostics) != 1L || !is.logical(diagnostics) ||
      is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(options)) stop("`options` must be a list.", call. = FALSE)
  for (ridge_name in c("ridge", "statistical_ridge")) {
    if (!is.null(options[[ridge_name]]) &&
        (!is.numeric(options[[ridge_name]]) ||
         length(options[[ridge_name]]) != 1L ||
         !is.finite(options[[ridge_name]]) || options[[ridge_name]] != 0)) {
      stop("Statistical ridge regularization is fixed at zero for these ",
           "decomposition comparisons.", call. = FALSE)
    }
  }
  prepared_by_caller <- inherits(observed, "decomposition_comparator_problem")
  problem <- if (prepared_by_caller) {
    .dc_validate_problem(observed)
  } else .dc_problem(observed)
  .dc_validate_components(initialization, problem$dims)
  if (length(initialization) > min(problem$dims)) {
    stop("The comparator rank cannot exceed a tensor mode dimension.",
         call. = FALSE)
  }
  if (!is.null(truth)) .dc_validate_components(truth, problem$dims, "truth")
  state <- .dc_solver_state()
  fitted <- switch(
    method,
    `Orth-ALS` = .dc_fit_orth_als(
      problem, initialization, iterations, truth, diagnostics, options, state
    ),
    ICO = .dc_fit_ico(
      problem, initialization, iterations, truth, diagnostics, options, state
    ),
    `NLS-GN` = .dc_fit_nls_gn(
      problem, initialization, iterations, truth, diagnostics, options, state
    ),
    `GN-CG` = .dc_fit_gn_cg(
      problem, initialization, iterations, truth, diagnostics, options, state
    ),
    `CP-ALS+ELS` = .dc_fit_als_els(
      problem, initialization, iterations, truth, diagnostics, options, state
    ),
    `Hybrid AMDM` = .dc_fit_hybrid_amdm(
      problem, initialization, iterations, truth, diagnostics, options, state
    )
  )
  objective <- if (length(fitted$objective_rows)) {
    do.call(rbind, fitted$objective_rows)
  } else data.frame(iteration = integer(), objective = numeric())
  if (length(fitted$error_rows)) {
    errors <- do.call(rbind, fitted$error_rows)
    objective$relative_frobenius_error <- errors$relative_frobenius_error[
      match(objective$iteration, errors$iteration)
    ]
  }
  if (diagnostics &&
      (!identical(as.integer(objective$iteration), 0:iterations) ||
       any(!is.finite(objective$objective)))) {
    stop("The comparator did not produce a complete finite objective path.",
         call. = FALSE)
  }
  accepted <- as.integer(fitted$accepted)
  list(
    method = method,
    components = fitted$components,
    objective = objective,
    actual_iterations = iterations,
    accepted_iterations = accepted,
    stop_reason = if (accepted == iterations) "iteration_budget" else
      "iteration_budget_with_rejected_outer_updates",
    diagnostics = fitted$extras,
    metadata = list(
      engine = "decomposition_comparators.R",
      reference = .dc_references[[method]],
      initialization = "caller_supplied_common_components",
      problem_prepared_by_caller = prepared_by_caller,
      iteration_definition = paste(
        "one full mode sweep for alternating methods or one outer",
        "trust-region/Gauss-Newton update"
      ),
      diagnostic_objective = paste(
        "0.5 * sum((direct_dense_reconstruction - observed)^2); evaluated",
        "only when diagnostics=TRUE"
      ),
      diagnostic_relative_frobenius_error = paste(
        "sqrt(sum((direct_dense_reconstruction - direct_dense_truth)^2) /",
        "sum(direct_dense_truth^2)); evaluated only when diagnostics=TRUE"
      ),
      diagnostic_metric_provenance = paste(
        "direct dense residual evaluation; avoids cancellation",
        "in ||estimate||^2 + ||truth||^2 - 2<estimate,truth>"
      ),
      requested_iterations = iterations,
      statistical_ridge = 0,
      intrinsic_damping = method %in% c("NLS-GN", "GN-CG"),
      pseudoinverse_fallbacks = state$pseudoinverse_fallbacks,
      diagnostics_retained = diagnostics,
      control = options
    )
  )
}
