manual_decomposition_sweep <- function(observed, init, alpha, update_order,
                                       refresh_residual = TRUE) {
  components <- SegreCP:::.validate_components(
    init, dim(observed), length(init)
  )
  caches <- lapply(components, function(comp) {
    SegreCP:::.rank1_vector(comp$lambda, comp$U)
  })
  initial_residual <- array(
    Reduce(`+`, caches) - as.numeric(observed), dim = dim(observed)
  )
  residual <- initial_residual

  for (i in update_order) {
    working_residual <- if (refresh_residual) residual else initial_residual
    old_cache <- caches[[i]]
    proposal <- SegreCP:::rgd_step_fused_cpp(
      working_residual, components[[i]]$U, components[[i]]$lambda,
      alpha, 0L
    )
    new_component <- SegreCP:::.canonicalize_component(
      list(lambda = proposal$lambda, U = proposal$U)
    )
    new_cache <- SegreCP:::.rank1_vector(
      new_component$lambda, new_component$U
    )
    if (refresh_residual) {
      residual <- residual + array(new_cache - old_cache, dim = dim(observed))
    }
    components[[i]] <- new_component
    caches[[i]] <- new_cache
  }
  components
}

test_that("both profiles default to unregularized RGN", {
  expect_identical(segre_control()$ridge, 0)
  expect_identical(segre_control(profile = "practical")$ridge, 0)
  expect_identical(segre_control(profile = "practical", ridge = 0.01)$ridge,
                   0.01)
})

component_arrays <- function(components) {
  lapply(components, function(comp) rank1_array(comp$lambda, comp$U))
}

test_that("rank-two decomposition uses immediate Gauss--Seidel residual refresh", {
  set.seed(31)
  dims <- c(2, 3, 4)
  init <- list(
    list(lambda = 0.9, U = lapply(dims, unit_vector)),
    list(lambda = 1.8, U = lapply(dims, unit_vector))
  )
  observed <- rank1_array(1.4, lapply(dims, unit_vector)) +
    rank1_array(0.8, lapply(dims, unit_vector)) +
    array(rnorm(prod(dims), sd = 0.04), dims)
  alpha <- 0.3

  expected <- manual_decomposition_sweep(
    observed, init, alpha, update_order = 1:2, refresh_residual = TRUE
  )
  stale_residual_update <- manual_decomposition_sweep(
    observed, init, alpha, update_order = 1:2, refresh_residual = FALSE
  )
  fit <- segre_decompose(
    observed, rank = 2, method = "rgd", init = init,
    control = segre_control(
      max_iter = 1, step_size = alpha, component_order = "fixed"
    )
  )

  expect_equal(component_arrays(fit$components), component_arrays(expected),
               tolerance = 2e-12)
  expect_gt(
    sqrt(sum((reconstruct_cp(expected) -
      reconstruct_cp(stale_residual_update))^2)),
    1e-6
  )
  expect_identical(fit$control$scheme, "sequential")
})

test_that("the practical profile follows decreasing initial component weight", {
  set.seed(32)
  dims <- c(2, 3, 4)
  init <- list(
    list(lambda = 0.65, U = lapply(dims, unit_vector)),
    list(lambda = 2.1, U = lapply(dims, unit_vector))
  )
  observed <- rank1_array(1.6, lapply(dims, unit_vector)) +
    rank1_array(0.7, lapply(dims, unit_vector))
  alpha <- 0.2

  expected <- manual_decomposition_sweep(
    observed, init, alpha, update_order = c(2L, 1L),
    refresh_residual = TRUE
  )
  fit <- segre_decompose(
    observed, rank = 2, method = "rgd", init = init,
    control = segre_control(
      profile = "practical", max_iter = 1, tolerance = 0,
      step_size = alpha
    )
  )
  fixed_order_fit <- segre_decompose(
    observed, rank = 2, method = "rgd", init = init,
    control = segre_control(
      profile = "practical", max_iter = 1, tolerance = 0,
      step_size = alpha, component_order = "fixed"
    )
  )

  expect_identical(fit$control$component_order, "decreasing_weight")
  expect_equal(component_arrays(fit$components), component_arrays(expected),
               tolerance = 2e-12)
  expect_gt(
    sqrt(sum((reconstruct_cp(fit) - reconstruct_cp(fixed_order_fit))^2)),
    1e-6
  )
})

test_that("regression exact-quadratic RGD matches its tangent-line oracle", {
  set.seed(33)
  dims <- c(2, 3)
  U0 <- lapply(dims, unit_vector)
  lambda0 <- 1.15
  init <- list(list(lambda = lambda0, U = U0))
  n <- 55
  X <- matrix(rnorm(n * prod(dims)), n)
  beta <- as.numeric(rank1_array(1.7, lapply(dims, unit_vector)))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.08))

  tangent_design <- SegreCP:::.compute_tangent_design(X, U0)
  residual <- as.numeric(X %*% as.numeric(rank1_array(lambda0, U0)) - y)
  rhs <- as.numeric(crossprod(tangent_design, residual))
  gram <- crossprod(tangent_design)
  numerator <- sum(rhs * rhs)
  denominator <- as.numeric(crossprod(rhs, gram %*% rhs))
  alpha <- min(max(numerator / denominator, 1e-10), 1e2)
  base <- c(lambda0, numeric(length(rhs) - 1L))
  expected <- SegreCP:::.retract_tangent_coordinates(
    U0, base - alpha * rhs
  )

  fit <- segre_regress(
    X, y, rank = 1, input_dims = dims, method = "rgd", init = init,
    control = segre_control(
      max_iter = 1, step_rule = "exact_quadratic",
      normalize_regression = "none"
    )
  )

  expect_equal(fit$accepted_steps[1, 1], alpha, tolerance = 2e-13)
  expect_equal(
    as.array(fit), rank1_array(expected$lambda, expected$U),
    tolerance = 2e-11
  )
})

test_that("sqrt-n normalization equals manually pre-scaled regression data", {
  set.seed(34)
  dims <- c(2, 3)
  init <- list(list(lambda = 0.95, U = lapply(dims, unit_vector)))
  n <- 48
  X <- matrix(rnorm(n * prod(dims)), n)
  beta <- as.numeric(rank1_array(1.5, lapply(dims, unit_vector)))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.06))
  alpha <- 0.22

  normalized <- segre_regress(
    X, y, rank = 1, input_dims = dims, method = "rgd", init = init,
    control = segre_control(
      max_iter = 2, step_size = alpha, normalize_regression = "sqrt_n"
    )
  )
  pre_scaled <- segre_regress(
    X / sqrt(n), y / sqrt(n), rank = 1, input_dims = dims,
    method = "rgd", init = init,
    control = segre_control(
      max_iter = 2, step_size = alpha, normalize_regression = "none"
    )
  )

  expect_equal(as.array(normalized), as.array(pre_scaled), tolerance = 2e-13)
  expect_equal(normalized$objective, pre_scaled$objective, tolerance = 2e-13)
  expect_equal(normalized$accepted_steps, pre_scaled$accepted_steps,
               tolerance = 2e-13)
  expect_equal(normalized$fitted_values / sqrt(n), pre_scaled$fitted_values,
               tolerance = 2e-13)
  expect_equal(normalized$residual_values / sqrt(n),
               pre_scaled$residual_values, tolerance = 2e-13)
})
