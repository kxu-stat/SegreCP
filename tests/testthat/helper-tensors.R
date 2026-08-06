unit_vector <- function(p) {
  z <- stats::rnorm(p)
  z / sqrt(sum(z * z))
}

rank1_array <- function(lambda, U) {
  array(lambda * Reduce(kronecker, rev(U)),
        vapply(U, length, integer(1L)))
}

orthogonal_cp_tensor <- function(dims, weights) {
  rank <- length(weights)
  factors <- lapply(dims, function(p) diag(p)[, seq_len(rank), drop = FALSE])
  components <- cp_components(weights, factors)
  list(array = reconstruct_cp(components), components = components)
}

explicit_rank1_hosvd <- function(x) {
  dims <- dim(x)
  U <- lapply(seq_along(dims), function(k) {
    perm <- c(k, setdiff(seq_along(dims), k))
    mat <- matrix(aperm(x, perm), nrow = dims[[k]])
    base::svd(mat, nu = 1L, nv = 0L)$u[, 1L]
  })
  direction <- as.numeric(Reduce(kronecker, rev(U)))
  lambda <- sum(as.numeric(x) * direction)
  if (lambda < 0) {
    lambda <- -lambda
    U[[1L]] <- -U[[1L]]
  }
  list(lambda = lambda, U = U)
}

