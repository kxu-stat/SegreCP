test_that("decomposition component history aligns with objective history", {
  set.seed(31)
  truth <- orthogonal_cp_tensor(c(4, 5, 3), c(3, 1.5))
  observed <- truth$array + array(rnorm(length(truth$array), sd = 0.03),
                                  dim(truth$array))
  init <- cp_initialize(observed, 2, "cpca", split = 1)$comps
  fit <- segre_decompose(
    observed, 2, method = "rgd", init = init,
    control = segre_control(max_iter = 3, step_size = 0.5),
    keep_history = TRUE
  )
  fit_without_history <- segre_decompose(
    observed, 2, method = "rgd", init = init,
    control = segre_control(max_iter = 3, step_size = 0.5),
    keep_history = FALSE
  )

  expect_length(fit$component_history, fit$iterations + 1L)
  expect_equal(reconstruct_cp(fit$component_history[[1L]]),
               reconstruct_cp(init), tolerance = 1e-13)
  expect_equal(reconstruct_cp(fit$component_history[[4L]]),
               as.array(fit), tolerance = 1e-13)
  expect_equal(as.array(fit), as.array(fit_without_history), tolerance = 0)
  expect_equal(fit$objective, fit_without_history$objective, tolerance = 0)
  expect_equal(fit$accepted_steps, fit_without_history$accepted_steps,
               tolerance = 0)

  trajectory <- convergence_data(fit, truth$components)
  expect_equal(trajectory$iteration, 0:3)
  expect_equal(trajectory$objective, fit$objective, tolerance = 0)
  expect_equal(trajectory$objective_ratio[[1L]], 1, tolerance = 0)
  expect_equal(
    trajectory$relative_objective_decrease[-1L],
    (fit$objective[-length(fit$objective)] - fit$objective[-1L]) /
      pmax(abs(fit$objective[-length(fit$objective)]),
           .Machine$double.eps),
    tolerance = 0
  )
  expect_true(all(is.finite(trajectory$relative_frobenius_error)))
})

test_that("early stopping trims component history to completed sweeps", {
  truth <- orthogonal_cp_tensor(c(3, 4, 5), 2)
  fit <- segre_decompose(
    truth$array, 1, method = "rgd", init = truth$components,
    control = segre_control(
      profile = "practical", max_iter = 5, tolerance = 1e-8
    ),
    keep_history = TRUE
  )
  expect_true(fit$converged)
  expect_lt(fit$iterations, 5L)
  expect_length(fit$component_history, fit$iterations + 1L)
  expect_equal(reconstruct_cp(tail(fit$component_history, 1L)[[1L]]),
               as.array(fit), tolerance = 0)
})

test_that("convergence reconstruction error requires stored components", {
  truth <- orthogonal_cp_tensor(c(3, 4, 5), c(2, 1))
  fit <- segre_decompose(
    truth$array, 2, method = "rgn", init = "hosvd",
    control = segre_control(max_iter = 1)
  )
  expect_null(fit$component_history)
  expect_s3_class(convergence_data(fit), "data.frame")
  expect_error(convergence_data(fit, truth$components), "keep_history")
  expect_error(convergence_data(list()), "fitted object")
})

test_that("regression history remains in public input-output mode order", {
  set.seed(32)
  input_dims <- c(2, 3)
  output_dims <- c(2, 4)
  public_dims <- c(input_dims, output_dims)
  truth <- list(list(lambda = 1.7, U = lapply(public_dims, unit_vector)))
  beta <- matrix(as.numeric(reconstruct_cp(truth)),
                 nrow = prod(input_dims), ncol = prod(output_dims))
  X <- matrix(rnorm(70 * prod(input_dims)), 70)
  Y <- X %*% beta
  fit <- segre_regress(
    X, Y, 1, input_dims, output_dims, method = "rgn", init = "hosvd",
    control = segre_control(max_iter = 2), keep_history = TRUE
  )
  fit_without_history <- segre_regress(
    X, Y, 1, input_dims, output_dims, method = "rgn", init = "hosvd",
    control = segre_control(max_iter = 2), keep_history = FALSE
  )

  expect_length(fit$component_history, 3L)
  expect_equal(vapply(fit$component_history[[1L]][[1L]]$U,
                      length, integer(1L)), public_dims)
  expect_equal(reconstruct_cp(fit$component_history[[3L]]),
               as.array(fit), tolerance = 1e-12)
  expect_equal(convergence_data(fit, truth)$iteration, 0:2)
  expect_equal(as.array(fit), as.array(fit_without_history), tolerance = 0)
  expect_equal(fit$objective, fit_without_history$objective, tolerance = 0)
  expect_equal(fit$normal_condition, fit_without_history$normal_condition,
               tolerance = 0)
})

test_that("keep_history is a validated logical flag", {
  truth <- orthogonal_cp_tensor(c(3, 3, 3), 2)
  expect_error(
    segre_decompose(truth$array, 1, keep_history = NA),
    "TRUE or FALSE"
  )
})
