test_that("decomposition RGN equals unit-step RGD", {
  set.seed(5)
  dims <- c(3, 4, 5)
  U0 <- lapply(dims, unit_vector)
  init <- list(list(lambda = 1.7, U = U0))
  observed <- rank1_array(2.2, lapply(dims, unit_vector)) +
    array(rnorm(prod(dims), sd = 0.02), dims)
  ctl <- segre_control(max_iter = 1, step_size = 1)
  rgd <- segre_decompose(observed, 1, "rgd", init = init, control = ctl)
  rgn <- segre_decompose(observed, 1, "rgn", init = init, control = ctl)
  expect_equal(as.array(rgd), as.array(rgn), tolerance = 2e-13)
})

test_that("one decomposition step matches direct projection and full T-HOSVD", {
  set.seed(6)
  dims <- c(2, 3, 4)
  U0 <- lapply(dims, unit_vector)
  init <- list(list(lambda = 1.4, U = U0))
  observed <- array(rnorm(prod(dims)), dims)
  alpha <- 0.35
  residual <- rank1_array(1.4, U0) - observed
  projected <- SegreCP:::tangent_proj_cpp(residual, U0)
  direct <- explicit_rank1_hosvd(rank1_array(1.4, U0) - alpha * projected)
  fit <- segre_decompose(
    observed, 1, "rgd", init = init,
    control = segre_control(max_iter = 1, step_size = alpha)
  )
  expect_equal(as.array(fit), rank1_array(direct$lambda, direct$U),
               tolerance = 2e-11)
})

test_that("decomposition accepts CPCA, HOSVD, and rTensor input", {
  truth <- orthogonal_cp_tensor(c(3, 4, 5), c(3, 1.5))
  ctl <- segre_control(max_iter = 1)
  a <- segre_decompose(truth$array, 2, init = "cpca", control = ctl)
  b <- segre_decompose(rTensor::as.tensor(truth$array), 2,
                       init = "hosvd", control = ctl)
  expect_s3_class(a, "segre_cp_decomposition")
  expect_s3_class(b, "segre_cp_decomposition")
  expect_lt(cp_relative_error(a, truth$components), 1e-9)
  expect_lt(cp_relative_error(b, truth$components), 1e-9)
})

