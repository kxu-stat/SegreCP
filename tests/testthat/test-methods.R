test_that("fit methods return consistent coefficients and predictions", {
  truth <- orthogonal_cp_tensor(c(3, 4, 5), c(2, 1))
  fit <- segre_decompose(truth$array, 2, "rgn", init = "hosvd",
                         control = segre_control(max_iter = 1))
  expect_equal(coef(fit), as.array(fit), tolerance = 0)
  expect_s4_class(as_rtensor(fit), "Tensor")
  expect_match(paste(capture.output(print(fit)), collapse = " "),
               "CP decomposition")
  expect_s3_class(summary(fit), "summary.segre_cp_fit")
})

test_that("regression fitted and residual methods use original scale", {
  set.seed(11)
  n <- 40
  dims <- c(2, 3)
  X <- matrix(rnorm(n * prod(dims)), n)
  y <- as.numeric(X %*% as.numeric(rank1_array(
    1.3, lapply(dims, unit_vector)
  )))
  fit <- segre_regress(X, y, 1, dims, method = "rgn", init = "hosvd",
                       control = segre_control(max_iter = 2))
  expect_equal(fitted(fit) + residuals(fit), y, tolerance = 1e-12)
  expect_equal(predict(fit, X), fitted(fit), tolerance = 1e-12)
})
