# Data generation for the tensor-on-tensor regression simulations: CP
# coefficient with coherent factors, Gaussian design with entry variance 1/n
# and Gaussian noise. Sourcing this file defines functions only.

.simulation_rng <- function(seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
}

.coherent_factors <- function(dimensions, rank, coherence) {
  lapply(dimensions, function(p) {
    gram <- stats::toeplitz(coherence^(0:(rank - 1L)))
    columns <- rbind(diag(1, rank), matrix(0, p - rank, rank)) %*%
      chol(gram)
    columns <- apply(columns, 2L, function(z) z / sqrt(sum(z^2)))
    rotation <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
    rotation %*% columns
  })
}

.tensor_regression_data <- function(
    input_dims, output_dims, rank, kappa, coherence, noise_sd,
    n_train, n_test, seed) {
  .simulation_rng(seed)
  generation_dims <- c(output_dims, input_dims)
  factors <- .coherent_factors(generation_dims, rank, coherence)
  internal_components <- lapply(seq_len(rank), function(component) {
    list(
      lambda = kappa^(0.5 - (component - 1L) / (rank - 1L)),
      U = lapply(factors, function(A) A[, component])
    )
  })
  truth_output_first <- .components_to_array(
    internal_components, generation_dims
  )
  if (length(output_dims)) {
    truth_public <- aperm(
      truth_output_first,
      c(length(output_dims) + seq_along(input_dims), seq_along(output_dims))
    )
    public_components <- lapply(internal_components, function(component) {
      component$U <- component$U[
        c(length(output_dims) + seq_along(input_dims), seq_along(output_dims))
      ]
      component
    })
  } else {
    truth_public <- truth_output_first
    public_components <- internal_components
  }
  p_in <- prod(input_dims)
  q_out <- if (length(output_dims)) prod(output_dims) else 1L
  X_train <- matrix(stats::rnorm(n_train * p_in, sd = 1 / sqrt(n_train)),
                    nrow = n_train, ncol = p_in)
  X_test <- matrix(stats::rnorm(n_test * p_in, sd = 1 / sqrt(n_train)),
                   nrow = n_test, ncol = p_in)
  coefficient_matrix <- matrix(as.numeric(truth_public), p_in, q_out)
  Y_train_clean <- unname(X_train %*% coefficient_matrix)
  Y_test_clean <- unname(X_test %*% coefficient_matrix)
  Z_train <- matrix(stats::rnorm(n_train * q_out), n_train, q_out)
  Z_test <- matrix(stats::rnorm(n_test * q_out), n_test, q_out)
  Y_train <- Y_train_clean + noise_sd * Z_train / sqrt(n_train)
  Y_test <- Y_test_clean + noise_sd * Z_test / sqrt(n_train)
  response <- function(Y) {
    if (length(output_dims)) array(as.numeric(Y), c(nrow(Y), output_dims))
    else as.numeric(Y)
  }
  list(
    X_train = X_train, X_test = X_test,
    Y_train = response(Y_train), Y_test = response(Y_test),
    Y_train_matrix = Y_train, Y_test_matrix = Y_test,
    Y_train_clean = response(Y_train_clean),
    Y_test_clean = response(Y_test_clean),
    Y_test_clean_matrix = Y_test_clean,
    truth = truth_public, components = public_components
  )
}

