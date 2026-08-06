test_that("array and rTensor inputs have identical storage and initialization", {
  set.seed(1)
  x <- array(rnorm(3 * 4 * 5), c(3, 4, 5))
  tx <- rTensor::as.tensor(x)
  expect_identical(as_segre_array(tx), x)

  a <- cp_initialize(x, 2, "cpca", split = c(1, 2))
  b <- cp_initialize(tx, 2, "cpca", split = c(1, 2))
  expect_equal(reconstruct_cp(a$comps), reconstruct_cp(b$comps), tolerance = 1e-12)
  expect_true(inherits(as_rtensor(a$comps |> reconstruct_cp()), "Tensor"))
})

test_that("CPCA and HOSVD recover an orthogonal CP tensor", {
  truth <- orthogonal_cp_tensor(c(3, 4, 5), c(4, 2))
  cpca <- cp_initialize(truth$array, 2, "cpca", split = 1L)
  hosvd <- cp_initialize(truth$array, 2, "hosvd")
  expect_lt(cp_relative_error(cpca$comps, truth$components), 1e-10)
  expect_lt(cp_relative_error(hosvd$comps, truth$components), 1e-10)
  expect_true(all(vapply(cpca$comps, `[[`, numeric(1), "lambda") > 0))
  expect_true(all(vapply(hosvd$comps, `[[`, numeric(1), "lambda") > 0))
})

test_that("automatic CPCA split is a deterministic proper subset", {
  set.seed(2)
  x <- array(rnorm(2^4), rep(2, 4))
  first <- cp_initialize(x, 1, "cpca")
  second <- cp_initialize(x, 1, "cpca")
  expect_identical(first$split, second$split)
  expect_gt(length(first$split), 0)
  expect_lt(length(first$split), length(dim(x)))
})

test_that("automatic CPCA split uses the manuscript product-scale criterion", {
  # Raw product imbalance and log-product imbalance select different subsets
  # for these dimensions.  Candidate enumeration breaks the raw-scale tie
  # deterministically at modes 2 and 3.
  x <- array(seq_len(24), c(3, 2, 2, 2))
  init <- cp_initialize(x, 1, "cpca")
  expect_identical(init$split, c(2L, 3L))
})

test_that("custom initialization accepts functions and components", {
  truth <- orthogonal_cp_tensor(c(3, 3, 3), c(2, 1))
  custom <- function(x, rank, supplied) supplied
  a <- cp_initialize(truth$array, 2, truth$components)
  b <- cp_initialize(truth$array, 2, custom, supplied = truth$components)
  expect_equal(reconstruct_cp(a$comps), truth$array, tolerance = 1e-12)
  expect_equal(reconstruct_cp(b$comps), truth$array, tolerance = 1e-12)
})

test_that("initialization rejects invalid ranks, splits, and zero components", {
  x <- array(rnorm(24), c(2, 3, 4))
  expect_error(cp_initialize(x, 3, "hosvd"), "rank")
  expect_error(cp_initialize(x, 1.5, "hosvd"), "positive integer")
  expect_error(cp_initialize(x, 1, "cpca", split = integer()), "proper subset")
  expect_error(cp_initialize(x, 1, "cpca", split = 1.5), "whole-number")
  expect_error(cp_initialize(x, 1, "cpca", split = 1:3), "proper subset")
  bad <- list(list(lambda = 0, U = list(c(1, 0), c(1, 0, 0), c(1, 0, 0, 0))))
  expect_error(cp_initialize(x, 1, bad), "zero weight")
  x[1] <- NA_real_
  expect_error(as_segre_array(x), "finite")
})

test_that("regression input validation rejects nonfinite data and fractional dimensions", {
  X <- matrix(rnorm(20), 5, 4)
  y <- rnorm(5)
  expect_error(
    segre_regress(X, y, rank = 1.5, input_dims = c(2, 2)),
    "positive integer"
  )
  expect_error(
    segre_regress(matrix(rnorm(25), 5, 5), y, rank = 1,
                  input_dims = c(2.5, 2)),
    "positive integers"
  )
  X[1, 1] <- Inf
  expect_error(
    segre_regress(X, y, rank = 1, input_dims = c(2, 2)),
    "Design entries"
  )
  X[1, 1] <- 0
  y[1] <- NA_real_
  expect_error(
    segre_regress(X, y, rank = 1, input_dims = c(2, 2)),
    "Response entries"
  )
  expect_error(
    segre_control(max_iter = 1.5),
    "positive integer"
  )
})
