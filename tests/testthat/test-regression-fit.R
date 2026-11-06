regression_problem <- function(tensor_response = TRUE) {
  set.seed(7103)
  n <- 64L
  input_dims <- c(3L, 4L)
  output_dims <- if (tensor_response) c(2L, 3L) else integer()
  dims <- c(input_dims, output_dims)
  truth <- lapply(c(2, .8), function(weight) {
    list(lambda = weight, U = lapply(dims, unit_vector))
  })
  initial <- lapply(truth, function(component) {
    list(lambda = .95 * component$lambda,
      U = lapply(component$U, function(u) {
        v <- u + rnorm(length(u), sd = .03)
        v / sqrt(sum(v * v))
      }))
  })
  X <- matrix(rnorm(n * prod(input_dims)), n) / sqrt(n)
  B <- matrix(reconstruct_cp(truth), nrow = prod(input_dims))
  Y <- X %*% B + matrix(rnorm(n * prod(output_dims), sd = .002), n)
  if (!tensor_response) Y <- as.numeric(Y)
  list(X = X, Y = Y, input_dims = input_dims, output_dims = output_dims,
       initial = initial)
}

test_that("a component start skips the moment and matches a function start", {
  original_moment <- SegreCP:::.coefficient_moment
  moment_calls <- 0L
  local_mocked_bindings(.coefficient_moment = function(...) {
    moment_calls <<- moment_calls + 1L
    original_moment(...)
  }, .package = "SegreCP")
  fields <- c("components", "component_history", "objective", "accepted_steps",
    "normal_condition", "iterations", "converged", "stopping_reason",
    "control", "normalization", "fitted_values", "residual_values",
    "response_source", "design_source")
  for (tensor_response in c(FALSE, TRUE)) {
    data <- regression_problem(tensor_response)
    for (normalization in c("none", "sqrt_n")) {
      for (route in c("fixed", "exact_quadratic", "backtracking", "rgn")) {
        method <- if (route == "rgn") "rgn" else "rgd"
        control <- segre_control(profile = "practical", max_iter = 4,
          tolerance = 0, step_size = .4,
          step_rule = if (route == "rgn") "fixed" else route,
          normalize_regression = normalization, ridge = 0)
        args <- list(X = data$X, Y = data$Y, rank = 2,
          input_dims = data$input_dims, output_dims = data$output_dims,
          method = method, init = data$initial, control = control,
          keep_history = TRUE, seed = 39)
        before_calls <- moment_calls
        before_seed <- .Random.seed
        supplied <- do.call(segre_regress, args)
        expect_identical(moment_calls, before_calls)
        expect_identical(.Random.seed, before_seed)
        expect_identical(supplied$initialization$method, "custom_components")
        # A function start forms the moment, then returns the same components.
        args$init <- function(x, rank) data$initial
        functional <- do.call(segre_regress, args)
        expect_identical(moment_calls, before_calls + 1L)
        expect_identical(supplied[fields], functional[fields])
        expect_identical(length(supplied$component_history), 5L)
        expect_identical(supplied$control$ridge, 0)
      }
    }
  }
})

test_that("an initializer function receives the scaled moment", {
  data <- regression_problem(TRUE)
  for (normalization in c("none", "sqrt_n")) {
    scale <- if (normalization == "sqrt_n") 1 / sqrt(nrow(data$X)) else 1
    expected <- array(crossprod(data$X * scale, data$Y * scale) / 2,
      c(data$input_dims, data$output_dims))
    calls <- 0L
    initialize <- function(x, rank, extra) {
      calls <<- calls + 1L
      expect_identical(x, expected)
      expect_identical(rank, 2L)
      expect_identical(extra, "retained")
      data$initial
    }
    fit <- segre_regress(data$X, data$Y, 2, data$input_dims,
      data$output_dims, method = "rgn", init = initialize,
      init_args = list(extra = "retained"), design_variance = 2,
      control = segre_control(max_iter = 1,
        normalize_regression = normalization))
    expect_identical(calls, 1L)
    expect_identical(fit$initialization$method, "custom_function")
  }
})

test_that("a component start takes no init_args and may be wrapped", {
  data <- regression_problem(TRUE)
  args <- list(X = data$X, Y = data$Y, rank = 2,
    input_dims = data$input_dims, output_dims = data$output_dims,
    method = "rgn", init = data$initial,
    control = segre_control(max_iter = 1, normalize_regression = "none"))
  expect_error(do.call(segre_regress, c(args, list(init_args = 1))),
    "`init_args` must be a list", fixed = TRUE)
  expect_error(do.call(segre_regress, c(args, list(init_args = list(seed = 1)))),
    "used only when `init` is a method name or a function", fixed = TRUE)
  expect_error(segre_decompose(reconstruct_cp(data$initial), 2,
    init = data$initial, init_args = list(seed = 1)),
    "used only when `init` is a method name or a function", fixed = TRUE)
  plain <- do.call(segre_regress, args)
  args$init <- list(comps = data$initial)
  wrapped <- do.call(segre_regress, args)
  expect_identical(plain$components, wrapped$components)
  expect_identical(plain$initialization, wrapped$initialization)
  expect_error(do.call(segre_regress, c(args, list(design_variance = 0))),
    "positive finite number", fixed = TRUE)
})

test_that("fitted values are on the scale of the data", {
  for (tensor_response in c(FALSE, TRUE)) {
    data <- regression_problem(tensor_response)
    for (normalization in c("none", "sqrt_n")) {
      for (method in c("rgd", "rgn")) {
        fit <- segre_regress(data$X, data$Y, 2, data$input_dims,
          data$output_dims, method = method, init = data$initial,
          control = segre_control(max_iter = 3, step_size = .4,
            normalize_regression = normalization))
        direct <- matrix(predict(fit, data$X), nrow(data$X))
        expect_equal(fit$fitted_values, direct, tolerance = 1e-12)
        expect_identical(fit$residual_values,
          matrix(data$Y, nrow(data$X)) - fit$fitted_values)
        expect_equal(tail(fit$objective, 1L),
          .5 * sum(fit$residual_values^2) * fit$normalization^2,
          tolerance = 1e-11)
      }
    }
  }
})

test_that("integer inputs give the same fit as double inputs", {
  data <- regression_problem(TRUE)
  X <- matrix(as.integer(round(20 * data$X)), nrow(data$X))
  Y <- matrix(as.integer(round(100 * data$Y)), nrow(data$X))
  ctl <- segre_control(max_iter = 2, normalize_regression = "none")
  a <- segre_regress(X, Y, 2, data$input_dims, data$output_dims,
    method = "rgn", init = data$initial, control = ctl, keep_history = TRUE)
  b <- segre_regress(X * 1, Y * 1, 2, data$input_dims, data$output_dims,
    method = "rgn", init = data$initial, control = ctl, keep_history = TRUE)
  expect_identical(a$component_history, b$component_history)
  expect_identical(a$objective, b$objective)
  expect_identical(a$fitted_values, b$fitted_values)
})

test_that("the block solve of the RGN normal equation equals the full solve", {
  set.seed(7104)
  for (dims in list(c(2L, 3L, 3L, 4L), c(1L, 2L, 1L, 3L),
                   c(2L, 3L, 1L, 1L))) {
    U <- lapply(dims, unit_vector)
    X <- matrix(rnorm(50 * prod(dims[3:4])), 50)
    target <- matrix(rnorm(50 * prod(dims[1:2])), 50)
    full <- SegreCP:::.tensor_tangent_statistics(X, target, U, dims, 2L,
      include_gram = TRUE)
    compact <- SegreCP:::.tensor_tangent_statistics(X, target, U, dims, 2L,
      include_gram = TRUE, full_gram = FALSE)
    expect_null(compact$gram)
    expect_identical(full[c("rhs", "gram_in", "df_in", "df_out")],
      compact[c("rhs", "gram_in", "df_in", "df_out")])
    for (ridge in c(0, .01)) {
      expected <- SegreCP:::.solve_tensor_tangent(full, ridge)
      actual <- SegreCP:::.solve_tensor_tangent(compact, ridge)
      expect_identical(actual, expected)
      expect_equal(actual, as.numeric(solve(full$gram +
        ridge * diag(nrow(full$gram)), full$rhs)), tolerance = 1e-11)
    }
  }
})
