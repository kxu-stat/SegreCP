.orth_complement_householder <- function(u) {
  u <- .normalize_vector(u, "factor")
  p <- length(u)
  if (p == 1L) return(list(v = u, coef = numeric()))
  sign1 <- if (u[[1L]] >= 0) 1 else -1
  v <- u
  v[[1L]] <- v[[1L]] + sign1
  vv <- sum(v * v)
  list(v = v, coef = 2 * v[-1L] / vv)
}

.orth_complement <- function(u) {
  p <- length(u)
  if (p == 1L) return(matrix(nrow = 1L, ncol = 0L))
  hh <- .orth_complement_householder(u)
  diag(p)[, -1L, drop = FALSE] - outer(hh$v, hh$coef)
}

.right_multiply_orth_complement <- function(F, u) {
  p <- length(u)
  if (p == 1L) return(matrix(nrow = nrow(F), ncol = 0L))
  hh <- .orth_complement_householder(u)
  F[, -1L, drop = FALSE] - outer(as.numeric(F %*% hh$v), hh$coef)
}

.lift_from_orth_complement <- function(u, z) {
  p <- length(u)
  if (length(z) != p - 1L) {
    stop("Tangent coordinate block has the wrong length.", call. = FALSE)
  }
  if (p == 1L) return(0)
  hh <- .orth_complement_householder(u)
  c(0, as.numeric(z)) - hh$v * sum(hh$coef * z)
}

.build_tangent_basis <- function(U) {
  U <- lapply(U, .normalize_vector)
  dims <- vapply(U, length, integer(1L))
  basis <- matrix(.rank1_vector(1, U), ncol = 1L)
  for (k in seq_along(U)) {
    Q <- .orth_complement(U[[k]])
    if (!ncol(Q)) next
    block <- vapply(seq_len(ncol(Q)), function(j) {
      Uk <- U
      Uk[[k]] <- Q[, j]
      .rank1_vector(1, Uk)
    }, numeric(prod(dims)))
    basis <- cbind(basis, block)
  }
  basis
}

.compute_tangent_design <- function(X, U) {
  dims <- vapply(U, length, integer(1L))
  df <- 1L + sum(dims - 1L)
  fibres <- modewise_design_fibres_cpp(X, U)
  A <- matrix(0, nrow(X), df)
  A[, 1L] <- as.numeric(fibres[[1L]] %*% U[[1L]])
  col <- 2L
  for (k in seq_along(U)) {
    pk <- dims[[k]]
    if (pk == 1L) next
    idx <- col:(col + pk - 2L)
    A[, idx] <- .right_multiply_orth_complement(fibres[[k]], U[[k]])
    col <- col + pk - 1L
  }
  A
}

.tangent_coords_to_vector <- function(U, coordinates) {
  dims <- vapply(U, length, integer(1L))
  expected <- 1L + sum(dims - 1L)
  if (length(coordinates) != expected) {
    stop(sprintf("Expected %d tangent coordinates; received %d.",
                 expected, length(coordinates)), call. = FALSE)
  }
  out <- coordinates[[1L]] * .rank1_vector(1, U)
  col <- 2L
  for (k in seq_along(U)) {
    block <- dims[[k]] - 1L
    if (!block) next
    idx <- col:(col + block - 1L)
    Uk <- U
    Uk[[k]] <- .lift_from_orth_complement(U[[k]], coordinates[idx])
    out <- out + .rank1_vector(1, Uk)
    col <- col + block
  }
  out
}

.retract_tangent_vectors <- function(U, a, w) {
  b <- vapply(w, function(wk) {
    scale <- max(abs(wk))
    if (!is.finite(scale) || scale == 0) 0 else
      scale * sqrt(sum((wk / scale)^2))
  }, numeric(1L))
  gram_scale <- max(c(abs(a), b))
  if (!is.finite(gram_scale) || gram_scale <= .Machine$double.eps) {
    stop("The T-HOSVD proposal is numerically zero.", call. = FALSE)
  }
  as <- a / gram_scale
  bs <- b / gram_scale
  total_b2 <- sum(bs * bs)
  ccoef <- scoef <- numeric(length(U))
  Unew <- vector("list", length(U))
  for (k in seq_along(U)) {
    gram <- matrix(c(as^2 + total_b2 - bs[[k]]^2, as * bs[[k]],
                     as * bs[[k]], bs[[k]]^2), 2L, 2L)
    z <- base::eigen(gram, symmetric = TRUE)$vectors[, 1L]
    if (z[[1L]] < 0 || (z[[1L]] == 0 && z[[2L]] < 0)) z <- -z
    ccoef[[k]] <- z[[1L]]
    scoef[[k]] <- z[[2L]]
    if (b[[k]] > 0) {
      Unew[[k]] <- ccoef[[k]] * U[[k]] + scoef[[k]] * w[[k]] / b[[k]]
    } else {
      Unew[[k]] <- U[[k]]
      ccoef[[k]] <- 1
      scoef[[k]] <- 0
    }
    Unew[[k]] <- .normalize_vector(Unew[[k]])
  }
  lambda <- a * prod(ccoef)
  for (k in seq_along(U)) {
    other <- if (length(U) == 1L) 1 else prod(ccoef[-k])
    lambda <- lambda + b[[k]] * scoef[[k]] * other
  }
  .canonicalize_component(list(lambda = as.numeric(lambda), U = Unew))
}

.unpack_tangent_coordinates <- function(U, coordinates) {
  dims <- vapply(U, length, integer(1L))
  expected <- 1L + sum(dims - 1L)
  if (length(coordinates) != expected) {
    stop(sprintf("Expected %d tangent coordinates; received %d.",
                 expected, length(coordinates)), call. = FALSE)
  }
  w <- vector("list", length(U))
  col <- 2L
  for (k in seq_along(U)) {
    block <- dims[[k]] - 1L
    if (block) {
      idx <- col:(col + block - 1L)
      w[[k]] <- .lift_from_orth_complement(U[[k]], coordinates[idx])
      col <- col + block
    } else {
      w[[k]] <- numeric(dims[[k]])
    }
  }
  list(a = as.numeric(coordinates[[1L]]), w = w)
}

.retract_tangent_coordinates <- function(U, coordinates) {
  tangent <- .unpack_tangent_coordinates(U, coordinates)
  .retract_tangent_vectors(U, tangent$a, tangent$w)
}

.prepare_tangent_retraction <- function(U, base, direction) {
  b <- .unpack_tangent_coordinates(U, base)
  d <- .unpack_tangent_coordinates(U, direction)
  list(U = U, a_base = b$a, a_direction = d$a,
       w_base = b$w, w_direction = d$w)
}

.apply_prepared_retraction <- function(context, alpha) {
  w <- Map(function(w0, wd) w0 + alpha * wd,
           context$w_base, context$w_direction)
  .retract_tangent_vectors(
    context$U,
    context$a_base + alpha * context$a_direction,
    w
  )
}

.split_internal_dims <- function(dims, output_order) {
  dims <- as.integer(dims)
  d <- length(dims)
  output_order <- as.integer(output_order)
  if (length(output_order) != 1L || is.na(output_order) || output_order < 0L ||
      output_order >= d) {
    stop("The number of output modes must lie between 0 and d - 1.",
         call. = FALSE)
  }
  out_idx <- if (output_order) seq_len(output_order) else integer()
  in_idx <- seq.int(output_order + 1L, d)
  list(
    dims = dims,
    d = d,
    output_order = output_order,
    out_idx = out_idx,
    in_idx = in_idx,
    out_dims = dims[out_idx],
    in_dims = dims[in_idx],
    q = as.integer(prod(dims[out_idx])),
    p_in = as.integer(prod(dims[in_idx]))
  )
}

.kron_factors <- function(component, split) {
  if (split$output_order) {
    out <- component$U[split$out_idx]
    u_out <- if (length(out) == 1L) out[[1L]] else
      Reduce(kronecker, rev(out))
  } else {
    u_out <- 1
  }
  inn <- component$U[split$in_idx]
  u_in <- if (length(inn) == 1L) inn[[1L]] else
    Reduce(kronecker, rev(inn))
  list(u_out = as.numeric(u_out), u_in = as.numeric(u_in))
}

.apply_rank1_operator <- function(component, X, split) {
  k <- .kron_factors(component, split)
  score <- as.numeric(X %*% k$u_in)
  prediction <- if (split$q == 1L) {
    matrix(component$lambda * score, ncol = 1L)
  } else {
    tcrossprod(component$lambda * score, k$u_out)
  }
  list(prediction = prediction, score = score,
       u_out = k$u_out, u_in = k$u_in)
}

.tensor_tangent_factors <- function(X, U, dims, output_order) {
  split <- .split_internal_dims(dims, output_order)
  if (!split$output_order) {
    stop("Tensor-response tangent factors require output modes.")
  }
  B_out <- .build_tangent_basis(U[split$out_idx])
  A_in <- .compute_tangent_design(X, U[split$in_idx])
  list(B_out = B_out, A_in = A_in, split = split)
}

.tensor_tangent_rhs <- function(target, B_out, A_in) {
  df_out <- ncol(B_out)
  df_in <- ncol(A_in)
  rhs <- numeric(df_out + df_in - 1L)
  u_out <- B_out[, 1L]
  score <- A_in[, 1L]
  rhs_in <- as.numeric(crossprod(A_in, target %*% u_out))
  rhs[[1L]] <- rhs_in[[1L]]
  if (df_out > 1L) {
    rhs[2L:df_out] <- as.numeric(crossprod(
      B_out[, -1L, drop = FALSE],
      as.numeric(crossprod(target, score))
    ))
  }
  if (df_in > 1L) {
    rhs[(df_out + 1L):(df_out + df_in - 1L)] <- rhs_in[-1L]
  }
  rhs
}

.tensor_tangent_statistics <- function(X, target, U, dims, output_order,
                                       include_gram = FALSE,
                                       full_gram = TRUE) {
  factors <- .tensor_tangent_factors(X, U, dims, output_order)
  rhs <- .tensor_tangent_rhs(target, factors$B_out, factors$A_in)
  if (!include_gram) return(list(rhs = rhs))
  gram_in <- crossprod(factors$A_in)
  df_out <- ncol(factors$B_out)
  df_in <- ncol(factors$A_in)
  df <- df_out + df_in - 1L
  gram <- NULL
  if (full_gram) {
    gram <- matrix(0, df, df)
    gram[1L, 1L] <- gram_in[1L, 1L]
    if (df_out > 1L) {
      oi <- 2L:df_out
      gram[oi, oi] <- diag(gram_in[1L, 1L], length(oi))
    }
    if (df_in > 1L) {
      ii <- (df_out + 1L):df
      gram[1L, ii] <- gram_in[1L, -1L]
      gram[ii, 1L] <- gram_in[-1L, 1L]
      gram[ii, ii] <- gram_in[-1L, -1L, drop = FALSE]
    }
  }
  list(rhs = rhs, gram = gram, gram_in = gram_in,
       df_out = df_out, df_in = df_in)
}

.solve_tensor_tangent <- function(statistics, ridge = 0) {
  df_out <- statistics$df_out
  df_in <- statistics$df_in
  df <- df_out + df_in - 1L
  result <- numeric(df)
  joint <- c(1L, if (df_in > 1L) (df_out + 1L):df else integer())
  gram <- statistics$gram_in
  if (ridge > 0) gram <- gram + ridge * diag(df_in)
  result[joint] <- tryCatch(
    solve(gram, statistics$rhs[joint]),
    error = function(e) stop(
      "The restricted tangent normal operator is singular; use a better ",
      "initialization or set a positive `ridge` in the practical profile.",
      call. = FALSE
    )
  )
  if (df_out > 1L) {
    denom <- statistics$gram_in[1L, 1L] + ridge
    if (!is.finite(denom) || denom <= .Machine$double.eps) {
      stop("The output tangent curvature is singular.", call. = FALSE)
    }
    result[2L:df_out] <- statistics$rhs[2L:df_out] / denom
  }
  result
}
