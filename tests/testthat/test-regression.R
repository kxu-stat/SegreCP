test_that("scalar-response RGN one-step matches explicit tangent least squares", {
  set.seed(7)
  dims <- c(2, 3)
  U0 <- lapply(dims, unit_vector)
  init <- list(list(lambda = 1.2, U = U0))
  n <- 40
  X <- matrix(rnorm(n * prod(dims)), n)
  beta <- as.numeric(rank1_array(1.8, lapply(dims, unit_vector)))
  y <- as.numeric(X %*% beta)

  V <- SegreCP:::.build_tangent_basis(U0)
  A <- X %*% V
  coordinates <- solve(crossprod(A), crossprod(A, y))
  expected <- SegreCP:::.retract_tangent_coordinates(U0, coordinates)
  fit <- segre_regress(
    X, y, 1, input_dims = dims, method = "rgn", init = init,
    control = segre_control(max_iter = 1, normalize_regression = "none")
  )
  expect_equal(as.array(fit), rank1_array(expected$lambda, expected$U),
               tolerance = 2e-11)
})

test_that("tensor-response structured RGN matches dense tangent oracle", {
  set.seed(8)
  input_dims <- c(2, 3)
  output_dims <- c(2, 4)
  public_dims <- c(input_dims, output_dims)
  init_public <- list(list(lambda = 1.1, U = lapply(public_dims, unit_vector)))
  n <- 50
  X <- matrix(rnorm(n * prod(input_dims)), n)
  truth <- list(list(lambda = 1.7, U = lapply(public_dims, unit_vector)))
  truth_internal <- SegreCP:::.public_to_internal_components(
    truth, length(input_dims), length(output_dims)
  )
  split <- SegreCP:::.split_internal_dims(c(output_dims, input_dims),
                                          length(output_dims))
  Y <- SegreCP:::.apply_rank1_operator(truth_internal[[1]], X, split)$prediction

  init_internal <- SegreCP:::.public_to_internal_components(
    init_public, length(input_dims), length(output_dims)
  )[[1]]
  V <- SegreCP:::.build_tangent_basis(init_internal$U)
  q <- prod(output_dims)
  p <- prod(input_dims)
  dense_design <- vapply(seq_len(ncol(V)), function(j) {
    Tmat <- matrix(V[, j], nrow = q, ncol = p)
    as.numeric(X %*% t(Tmat))
  }, numeric(n * q))
  coordinates <- solve(crossprod(dense_design),
                       crossprod(dense_design, as.numeric(Y)))
  expected_internal <- SegreCP:::.retract_tangent_coordinates(
    init_internal$U, coordinates
  )
  expected_public <- SegreCP:::.internal_to_public_components(
    list(expected_internal), length(input_dims), length(output_dims)
  )[[1]]

  fit <- segre_regress(
    X, Y, 1, input_dims, output_dims, method = "rgn", init = init_public,
    control = segre_control(max_iter = 1, normalize_regression = "none")
  )
  expect_equal(as.array(fit),
               rank1_array(expected_public$lambda, expected_public$U),
               tolerance = 5e-10)
  expect_equal(dim(predict(fit, X[1:3, , drop = FALSE])),
               c(3, output_dims))
})

test_that("practical RGD backtracking is finite and non-increasing", {
  set.seed(9)
  dims <- c(2, 3)
  n <- 60
  X <- matrix(rnorm(n * prod(dims)), n)
  beta <- as.numeric(rank1_array(1.5, lapply(dims, unit_vector)))
  y <- as.numeric(X %*% beta + rnorm(n, sd = 0.1))
  fit <- segre_regress(
    X, y, 1, input_dims = dims, method = "rgd", init = "hosvd",
    control = segre_control(profile = "practical", max_iter = 4,
                            normalize_regression = "sqrt_n")
  )
  expect_true(all(is.finite(fit$objective)))
  expect_true(all(diff(fit$objective) <= 1e-10))
  expect_true(all(is.finite(fit$accepted_steps)))
})

test_that("regression array and rTensor sample stacks are compatible", {
  set.seed(10)
  n <- 30
  input_dims <- c(2, 3)
  Xarray <- array(rnorm(n * prod(input_dims)), c(n, input_dims))
  beta <- as.numeric(rank1_array(1.2, lapply(input_dims, unit_vector)))
  y <- matrix(matrix(as.numeric(Xarray), n) %*% beta, ncol = 1)
  ctl <- segre_control(max_iter = 1)
  a <- segre_regress(Xarray, y, 1, method = "rgn", init = "hosvd",
                     control = ctl)
  b <- segre_regress(rTensor::as.tensor(Xarray), rTensor::as.tensor(y), 1,
                     method = "rgn", init = "hosvd", control = ctl)
  expect_equal(as.array(a), as.array(b), tolerance = 1e-12)
})

