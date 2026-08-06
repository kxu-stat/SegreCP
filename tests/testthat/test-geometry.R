test_that("compiled tangent projection matches an explicit orthonormal basis", {
  set.seed(3)
  U <- lapply(c(2, 3, 4), unit_vector)
  dims <- vapply(U, length, integer(1L))
  x <- array(rnorm(prod(dims)), dims)
  V <- SegreCP:::.build_tangent_basis(U)
  direct <- as.numeric(V %*% crossprod(V, as.numeric(x)))
  compiled <- as.numeric(SegreCP:::tangent_proj_cpp(x, U))
  expect_equal(compiled, direct, tolerance = 2e-12)
  expect_equal(crossprod(V), diag(ncol(V)), tolerance = 2e-12)
})

test_that("structured T-HOSVD equals explicit unfolding SVD", {
  set.seed(4)
  U <- lapply(c(3, 4, 2, 3), unit_vector)
  df <- 1L + sum(vapply(U, length, integer(1L)) - 1L)
  coords <- c(2, rnorm(df - 1L, sd = 0.15))
  proposal <- array(
    SegreCP:::.tangent_coords_to_vector(U, coords),
    vapply(U, length, integer(1L))
  )
  structured <- SegreCP:::.retract_tangent_coordinates(U, coords)
  explicit <- explicit_rank1_hosvd(proposal)
  expect_equal(rank1_array(structured$lambda, structured$U),
               rank1_array(explicit$lambda, explicit$U), tolerance = 2e-11)

  tiny <- c(2, rep(sqrt(.Machine$double.eps) / 10, df - 1L))
  tiny_fit <- SegreCP:::.retract_tangent_coordinates(U, tiny)
  expect_true(all(is.finite(unlist(tiny_fit))))
})

