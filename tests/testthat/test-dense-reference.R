# The public API against dense reference implementations: explicit tangent
# projections, unfolding SVDs and least squares, without the compiled sweep or
# the structured regression statistics.
dense_components <- function(dims, rank) {
  lapply(seq_len(rank), function(i) {
    list(lambda = 2 + i / rank, U = lapply(dims, unit_vector))
  })
}

dense_problem <- function(dims, rank, seed) {
  set.seed(seed)
  truth <- dense_components(dims, rank)
  initial <- lapply(truth, function(component) {
    list(lambda = .96 * component$lambda, U = lapply(component$U, function(u) {
      v <- u + rnorm(length(u), sd = .01)
      v / sqrt(sum(v * v))
    }))
  })
  list(initial = initial, observed = reconstruct_cp(truth) +
    array(rnorm(prod(dims), sd = .03), dims))
}

dense_dense_tensor <- function(components) {
  Reduce(`+`, lapply(components, function(z) rank1_array(z$lambda, z$U)))
}

dense_dense_projection <- function(residual, component) {
  U <- component$U
  d <- length(U)
  direction <- rank1_array(1, U)
  full <- sum(residual * direction)
  projected <- full * direction
  for (k in seq_len(d)) {
    unfolded <- matrix(aperm(residual, c(k, setdiff(seq_len(d), k))),
      nrow = length(U[[k]]))
    other <- Reduce(kronecker, rev(U[-k]))
    contraction <- as.numeric(unfolded %*% other)
    factors <- U
    factors[[k]] <- contraction - full * U[[k]]
    projected <- projected + rank1_array(1, factors)
  }
  projected
}

dense_dense_decomposition <- function(observed, initial, method, control) {
  components <- initial
  residual <- dense_dense_tensor(components) - observed
  history <- list(components)
  objective <- .5 * sum(residual^2)
  steps <- matrix(NA_real_, control$max_iter, length(components))
  for (iteration in seq_len(control$max_iter)) {
    order_i <- if (control$component_order == "fixed") seq_along(components) else
      order(vapply(components, `[[`, numeric(1), "lambda"), decreasing = TRUE)
    alpha <- if (method == "rgn") 1 else if (length(control$step_size) > 1L)
      control$step_size[[iteration]] else control$step_size
    if (method == "rgd" && control$step_rule == "diminishing")
      alpha <- alpha / sqrt(iteration)
    for (i in order_i) {
      old <- components[[i]]
      old_array <- rank1_array(old$lambda, old$U)
      projected <- dense_dense_projection(residual, old)
      step <- alpha
      if (method == "rgd" && control$step_rule == "backtracking") {
        loss_old <- .5 * sum(residual^2)
        best <- list(component = old, residual = residual, loss = loss_old, step = 0)
        for (attempt in seq_len(control$max_backtracks)) {
          proposal <- explicit_rank1_hosvd(old_array - step * projected)
          candidate <- residual + rank1_array(proposal$lambda, proposal$U) - old_array
          loss_new <- .5 * sum(candidate^2)
          if (loss_new < best$loss)
            best <- list(component = proposal, residual = candidate, loss = loss_new, step = step)
          if (loss_new <= loss_old - control$armijo * step * sum(projected^2)) break
          step <- step * control$backtrack_factor
        }
        components[[i]] <- best$component
        residual <- best$residual
        step <- best$step
      } else {
        proposal <- explicit_rank1_hosvd(old_array - step * projected)
        components[[i]] <- proposal
        residual <- residual + rank1_array(proposal$lambda, proposal$U) - old_array
      }
      steps[iteration, i] <- step
    }
    history[[iteration + 1L]] <- components
    objective <- c(objective, .5 * sum(residual^2))
  }
  list(history = history, objective = objective, accepted_steps = steps)
}

dense_expect_decomposition <- function(fit, observed, reference = NULL) {
  expect_s3_class(fit, "segre_cp_decomposition")
  expect_identical(fit$control$ridge, 0)
  expect_identical(length(fit$component_history), fit$iterations + 1L)
  expect_identical(dim(fit$accepted_steps), c(fit$iterations, fit$rank))
  expect_true(all(is.finite(fit$objective)))
  expect_true(all(is.finite(fit$accepted_steps)))
  for (iteration in 0:fit$iterations) {
    current <- dense_dense_tensor(fit$component_history[[iteration + 1L]])
    expect_equal(fit$objective[[iteration + 1L]], .5 * sum((current - observed)^2),
      tolerance = 2e-10)
    if (!is.null(reference)) {
      target <- dense_dense_tensor(reference$history[[iteration + 1L]])
      expect_lt(sqrt(sum((current - target)^2)) / max(1, sqrt(sum(target^2))), 2e-10)
    }
    for (component in fit$component_history[[iteration + 1L]]) {
      expect_true(is.finite(component$lambda) && component$lambda > 0)
      expect_equal(vapply(component$U, function(u) sum(u^2), numeric(1)),
        rep(1, length(dim(observed))), tolerance = 2e-13)
    }
  }
  if (!is.null(reference))
    expect_equal(fit$accepted_steps, reference$accepted_steps, tolerance = 1e-13)
}

test_that("decomposition matches the dense reference for orders two to six", {
  for (d in 2:6) for (rank in c(1L, 2L, 3L, 5L, 8L)) {
    dims <- rep(c(2L, 3L, 4L), length.out = d)
    problem <- dense_problem(dims, rank, 9100 + 100 * d + rank)
    initial_before <- serialize(problem$initial, NULL)
    observed_before <- serialize(problem$observed, NULL)
    for (method in c("rgd", "rgn")) {
      ctl <- segre_control(max_iter = 3, step_size = .25, step_rule = "fixed",
        component_order = if (rank %% 2L) "fixed" else "decreasing_weight", ridge = 0)
      fit <- segre_decompose(problem$observed, rank, method,
        init = problem$initial, control = ctl, keep_history = TRUE)
      reference <- dense_dense_decomposition(problem$observed,
        fit$component_history[[1L]], method, ctl)
      dense_expect_decomposition(fit, problem$observed, reference)
      expect_identical(fit$iterations, 3L)
      expect_identical(fit$stopping_reason, "max_iter")
    }
    expect_identical(serialize(problem$initial, NULL), initial_before)
    expect_identical(serialize(problem$observed, NULL), observed_before)
  }
})

test_that("singleton modes and custom overcomplete ranks remain valid", {
  cases <- list(list(dims = c(1L, 3L), rank = 1L),
    list(dims = c(1L, 3L, 2L), rank = 2L),
    list(dims = c(2L, 1L, 3L, 1L), rank = 3L),
    list(dims = c(1L, 2L, 1L, 3L, 2L), rank = 5L),
    list(dims = c(2L, 1L, 3L, 1L, 2L, 1L), rank = 8L))
  for (j in seq_along(cases)) {
    case <- cases[[j]]
    problem <- dense_problem(case$dims, case$rank, 9400 + j)
    for (method in c("rgd", "rgn")) {
      ctl <- segre_control(max_iter = 3, step_size = .2, step_rule = "fixed", ridge = 0)
      fit <- segre_decompose(problem$observed, case$rank, method,
        init = problem$initial, control = ctl, keep_history = TRUE)
      reference <- dense_dense_decomposition(problem$observed,
        fit$component_history[[1L]], method, ctl)
      dense_expect_decomposition(fit, problem$observed, reference)
    }
  }
})

test_that("all initialization routes work across tensor orders", {
  for (d in 2:6) {
    dims <- rep(c(3L, 2L), length.out = d)
    problem <- dense_problem(dims, 2L, 9500 + d)
    custom <- function(x, rank, supplied, tag) {
      expect_equal(x, problem$observed, tolerance = 0)
      expect_identical(rank, 2L)
      expect_identical(tag, "general_api")
      list(comps = supplied)
    }
    routes <- list(cpca = "cpca", hosvd = "hosvd", random = "random",
      components = problem$initial, custom_function = custom)
    for (route in names(routes)) for (method in c("rgd", "rgn")) {
      before_rng <- .Random.seed
      args <- if (route == "custom_function")
        list(supplied = problem$initial, tag = "general_api") else list()
      fit <- segre_decompose(problem$observed, 2, method, init = routes[[route]],
        init_args = args, seed = 19,
        control = segre_control(max_iter = 2, step_size = .2, ridge = 0),
        keep_history = TRUE)
      expect_identical(.Random.seed, before_rng)
      dense_expect_decomposition(fit, problem$observed)
      expected_method <- switch(route, components = "custom_components",
        custom_function = "custom_function", route)
      expect_identical(fit$initialization$method, expected_method)
    }
  }
})

test_that("fixed schedules diminishing steps and backtracking preserve their updates", {
  for (d in c(2L, 4L, 6L)) {
    problem <- dense_problem(rep(c(2L, 3L), length.out = d), 3L, 9600 + d)
    for (rule in c("fixed", "diminishing", "backtracking")) {
      schedule <- if (rule == "backtracking") c(8, 6, 4) else c(.4, .3, .2)
      ctl <- segre_control(max_iter = 3, step_size = schedule, step_rule = rule,
        component_order = "decreasing_weight", ridge = 0, max_backtracks = 10)
      fit <- segre_decompose(problem$observed, 3, "rgd", init = problem$initial,
        control = ctl, keep_history = TRUE)
      reference <- dense_dense_decomposition(problem$observed,
        fit$component_history[[1L]], "rgd", ctl)
      dense_expect_decomposition(fit, problem$observed, reference)
      if (rule == "backtracking") {
        expect_true(all(diff(fit$objective) <= 1e-10))
        expect_true(all(fit$accepted_steps >= 0))
        expect_true(all(fit$accepted_steps <= schedule))
      }
    }
  }
})

test_that("RGN remains unit-step RGD under general dimensions", {
  for (d in 2:6) {
    problem <- dense_problem(rep(c(2L, 3L), length.out = d), 2L, 9700 + d)
    rgd <- segre_decompose(problem$observed, 2, "rgd", init = problem$initial,
      control = segre_control(max_iter = 3, step_size = 1), keep_history = TRUE)
    rgn <- segre_decompose(problem$observed, 2, "rgn", init = problem$initial,
      control = segre_control(max_iter = 3, step_size = c(.1, .2, .3),
        step_rule = "diminishing"), keep_history = TRUE)
    expect_equal(rgd$component_history, rgn$component_history, tolerance = 0)
    expect_identical(rgd$objective, rgn$objective)
    expect_identical(rgd$accepted_steps, rgn$accepted_steps)
  }
})

test_that("aliases, history and inputs are consistent", {
  problem <- dense_problem(c(3L, 2L, 4L), 2L, 9801)
  dimnames(problem$observed) <- lapply(dim(problem$observed), function(n) letters[seq_len(n)])
  before <- serialize(problem, NULL)
  control <- segre_control(max_iter = 3, ridge = 0)
  fit <- segre_decompose(problem$observed, 2, "rgn", init = problem$initial,
    control = control, keep_history = TRUE)
  tensor_fit <- segre_decompose(rTensor::as.tensor(problem$observed), 2, "rgn",
    init = list(comps = problem$initial), control = control, keep_history = TRUE)
  plain <- segre_decompose(problem$observed, 2, "rgn", init = problem$initial,
    control = control, keep_history = FALSE)
  expect_identical(serialize(problem, NULL), before)
  expect_equal(fit$component_history, tensor_fit$component_history, tolerance = 0)
  expect_equal(fit$components, plain$components, tolerance = 0)
  expect_identical(fit$objective, plain$objective)
  expect_null(plain$component_history)
  expect_equal(coef(fit), as.array(fit), tolerance = 0)
  expect_s4_class(as_rtensor(fit), "Tensor")
  integer_x <- array(as.integer(round(20 * problem$observed)), dim(problem$observed))
  integer_fit <- segre_decompose(integer_x, 2, "rgd", init = problem$initial,
    control = control, keep_history = TRUE)
  double_fit <- segre_decompose(integer_x * 1, 2, "rgd", init = problem$initial,
    control = control, keep_history = TRUE)
  expect_equal(integer_fit$component_history, double_fit$component_history, tolerance = 0)
  expect_identical(integer_fit$objective, double_fit$objective)
})

test_that("stopping and invalid-input boundaries are explicit", {
  dims <- c(2L, 1L, 3L)
  initial <- list(list(lambda = 2, U = lapply(dims, function(p) c(1, rep(0, p - 1L)))))
  observed <- reconstruct_cp(initial)
  for (method in c("rgd", "rgn")) {
    stopped <- segre_decompose(observed, 1, method, init = initial,
      control = segre_control(max_iter = 5, tolerance = 1e-8), keep_history = TRUE)
    full <- segre_decompose(observed, 1, method, init = initial,
      control = segre_control(max_iter = 5, tolerance = 0), keep_history = TRUE)
    expect_identical(stopped$iterations, 1L)
    expect_true(stopped$converged)
    expect_identical(stopped$stopping_reason, "relative_objective_change")
    expect_identical(length(stopped$component_history), 2L)
    expect_identical(full$iterations, 5L)
    expect_false(full$converged)
    expect_true(all(full$objective == 0))
    expect_error(segre_decompose(array(0, dims), 1, method, init = initial), "zero observed")
    bad <- observed; bad[1L] <- Inf
    expect_error(segre_decompose(bad, 1, method, init = initial), "finite")
  }
  expect_error(cp_initialize(observed, 2, "hosvd"), "rank")
  expect_error(cp_initialize(observed, 2, "cpca", split = 2L), "rank")
  expect_error(segre_decompose(observed, 1, init = initial,
    control = segre_control(max_iter = 3, step_size = c(.2, .1))), "schedule")
  expect_error(segre_decompose(observed, 1, "rgd", init = initial,
    control = segre_control(max_iter = 1, step_rule = "exact_quadratic")), "regression RGD")
})

dense_dense_regression <- function(X, Y, initial, method, control, input_dims) {
  components <- initial
  dims <- vapply(initial[[1L]]$U, length, integer(1))
  p <- prod(input_dims)
  n <- nrow(X)
  scale <- if (control$normalize_regression == "sqrt_n") 1 / sqrt(n) else 1
  X <- X * scale; Y <- matrix(Y, nrow = n) * scale
  predict_component <- function(component) X %*%
    matrix(rank1_array(component$lambda, component$U), nrow = p)
  predictions <- lapply(components, predict_component)
  residual <- Reduce(`+`, predictions) - Y
  history <- list(components)
  for (iteration in seq_len(control$max_iter)) {
    order_i <- if (control$component_order == "fixed") seq_along(components) else
      order(vapply(components, `[[`, numeric(1), "lambda"), decreasing = TRUE)
    for (i in order_i) {
      old <- components[[i]]
      if (method == "rgn") {
        # Dense tangent design and least squares by QR.
        V <- SegreCP:::.build_tangent_basis(old$U)
        design <- vapply(seq_len(ncol(V)), function(j)
          as.numeric(X %*% matrix(V[, j], nrow = p)), numeric(length(Y)))
        target <- Y - (Reduce(`+`, predictions) - predictions[[i]])
        coordinates <- qr.solve(design, as.numeric(target))
        proposal <- array(as.numeric(V %*% coordinates), dims)
      } else {
        gradient <- array(as.numeric(crossprod(X, residual)), dims)
        projected <- dense_dense_projection(gradient, old)
        numerator <- sum(projected^2)
        denominator <- sum((X %*% matrix(projected, nrow = p))^2)
        alpha <- if (is.finite(denominator) && denominator > .Machine$double.eps)
          numerator / denominator else control$step_size
        alpha <- min(max(alpha, 1e-10), 1e2)
        proposal <- rank1_array(old$lambda, old$U) - alpha * projected
      }
      next_component <- explicit_rank1_hosvd(proposal)
      next_prediction <- predict_component(next_component)
      residual <- residual + next_prediction - predictions[[i]]
      components[[i]] <- next_component
      predictions[[i]] <- next_prediction
    }
    history[[iteration + 1L]] <- components
  }
  history
}

test_that("regression matches the dense reference", {
  cases <- list(list(input = c(2L, 3L), output = integer(), rank = 1L),
    list(input = c(2L, 3L), output = 2L, rank = 2L),
    list(input = c(2L, 1L, 3L), output = 1L, rank = 3L),
    list(input = c(2L, 3L), output = c(2L, 1L, 3L), rank = 5L),
    list(input = c(2L, 1L, 3L), output = c(2L, 3L, 1L), rank = 8L))
  for (j in seq_along(cases)) {
    case <- cases[[j]]; dims <- c(case$input, case$output)
    problem <- dense_problem(dims, case$rank, 9900 + j)
    p <- prod(case$input); n <- p + 20L
    X <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
    Y <- X %*% matrix(problem$observed, nrow = p)
    if (!length(case$output)) Y <- as.numeric(Y)
    before <- serialize(list(X, Y, problem$initial), NULL)
    for (normalization in c("none", "sqrt_n")) for (method in c("rgd", "rgn")) {
      ctl <- segre_control(max_iter = 2, ridge = 0, step_size = .3,
        step_rule = if (method == "rgd") "exact_quadratic" else "fixed",
        normalize_regression = normalization)
      fit <- segre_regress(X, Y, case$rank, case$input, case$output,
        method = method, init = problem$initial, control = ctl, keep_history = TRUE)
      reference <- dense_dense_regression(X, Y, fit$component_history[[1L]],
        method, ctl, case$input)
      expect_identical(fit$iterations, 2L)
      expect_identical(fit$control$ridge, 0)
      expect_identical(fit$dims, dims)
      expect_identical(length(fit$component_history), 3L)
      for (iteration in 0:2) {
        coefficient <- dense_dense_tensor(fit$component_history[[iteration + 1L]])
        expected <- dense_dense_tensor(reference[[iteration + 1L]])
        expect_lt(sqrt(sum((coefficient - expected)^2)) / max(1, sqrt(sum(expected^2))), 3e-10)
        residual <- matrix(Y, n) - X %*% matrix(coefficient, nrow = p)
        expect_equal(fit$objective[[iteration + 1L]],
          .5 * sum(residual^2) * fit$normalization^2, tolerance = 2e-10)
      }
      direct <- X %*% matrix(as.array(fit), nrow = p)
      expect_equal(fit$fitted_values, direct, tolerance = 2e-11)
      expect_equal(matrix(predict(fit, X), n), direct, tolerance = 2e-11)
      expect_identical(fit$residual_values, matrix(Y, n) - fit$fitted_values)
      expect_true(all(is.finite(fit$objective)))
      if (method == "rgn") expect_true(all(is.finite(fit$normal_condition)))
    }
    expect_identical(serialize(list(X, Y, problem$initial), NULL), before)
  }
})
