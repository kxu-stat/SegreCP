# Tucker RGN for tensor-on-tensor regression.
#
# Port of RGN_tenreg() with retrac = 'hosvd' from Luo and Zhang's public code
# (github.com/yuetianluo/RGN-for-Tensor-Estimation, RGN_Funs.R), extended from
# order-3 scalar regression to an order-d coefficient with input modes followed
# by output modes.  Each output entry y_{m,j} is a scalar measurement with
# sensing tensor x_m (outer) e_j.  The iteration is the same as in the original:
#   1. tangent coordinates: core B (r_1 x ... x r_d) and arms D_k ((p_k - r_k) x r_k),
#      with sketching directions V_k from the current core (their V1t, V2t, V3t);
#   2. least squares for (B, D_1, ..., D_d)  (their sketch_method);
#   3. HOSVD retraction of B x_k U_k + sum_k fold_k(U_k_perp D_k V_k' (kron U)').
# Only the computation of the least-squares normal equations differs: we form
# them from Kronecker structure instead of building the n*q-row design matrix.
#
# Data convention (as in regression_adapters.R): X is n x prod(input
# dims) with rows vec(x_m); Y is n x prod(output dims) with rows vec(y_m).

unfold <- function(x, k) {
  d <- length(dim(x))
  matrix(aperm(x, c(k, setdiff(seq_len(d), k))), nrow = dim(x)[k])
}

fold <- function(M, k, dims) {
  d <- length(dims)
  perm <- c(k, setdiff(seq_len(d), k))
  aperm(array(M, dims[perm]), order(perm))
}

mode_product <- function(x, M, k) {        # x times_k M
  dims <- dim(x)
  dims[k] <- nrow(M)
  fold(M %*% unfold(x, k), k, dims)
}

multi_mode_product <- function(x, mats) {  # x times_k mats[[k]] for all k
  for (k in seq_along(mats)) x <- mode_product(x, mats[[k]], k)
  x
}

# Kronecker product of a list, the first matrix varying fastest in the rows
# and columns (matches column-major vec of a tensor).
kron_fastest_first <- function(mats) {
  out <- mats[[1]]
  for (M in mats[-1]) out <- kronecker(M, out)
  out
}

hosvd_truncate <- function(x, ranks) {
  U <- lapply(seq_along(ranks), function(k) {
    svd(unfold(x, k), nu = ranks[k], nv = 0)$u
  })
  list(U = U, est = multi_mode_product(x, lapply(U, tcrossprod)))
}

complement <- function(U) {
  if (ncol(U) == nrow(U)) return(matrix(0, nrow(U), 0))   # full rank: empty
  qr.Q(qr(U), complete = TRUE)[, (ncol(U) + 1):nrow(U), drop = FALSE]
}

# Contract every x_m with one matrix per input mode: returns n x prod(cols).
input_contraction <- function(X, input_dims, mats) {
  n <- nrow(X)
  arr <- array(X, c(n, input_dims))
  for (k in seq_along(input_dims)) {
    arr <- mode_product(arr, t(mats[[k]]), k + 1)
  }
  matrix(arr, nrow = n)
}

# Rows e_j contracted with one matrix per output mode: q x prod(cols).
output_contraction <- function(mats) kron_fastest_first(mats)

# Mixing matrix for arm k: maps (e, c) to the vec index of a tensor whose
# mode-k size is p_k - r_k and other modes have size r_l.
arm_mixing <- function(k, sizes, V) {
  d <- length(sizes)
  others <- setdiff(seq_len(d), k)
  n_other <- prod(sizes[others])
  Q <- matrix(0, prod(sizes), sizes[k] * ncol(V))
  # vec index of (i_1, ..., i_d)
  strides <- cumprod(c(1, sizes))[seq_len(d)]
  other_grid <- as.matrix(expand.grid(lapply(sizes[others], seq_len)))
  base <- 1 + as.numeric((other_grid - 1) %*% strides[others])
  for (e in seq_len(sizes[k])) {
    rows <- base + (e - 1) * strides[k]
    for (c in seq_len(ncol(V))) {
      Q[rows, e + (c - 1) * sizes[k]] <- V[, c]
    }
  }
  Q
}

tucker_rgn_tot <- function(X, Y, input_dims, output_dims, ranks,
                           iterations, truth = NULL, keep_history = TRUE) {
  dims <- c(input_dims, output_dims)
  d <- length(dims)
  d_in <- length(input_dims)
  in_modes <- seq_len(d_in)
  out_modes <- d_in + seq_along(output_dims)

  # Initialization as in the original example: HOSVD of A*(y).
  W <- array(crossprod(X, Y), dims)
  init <- hosvd_truncate(W, ranks)
  U <- init$U
  Xt <- init$est
  history <- if (keep_history) list(Xt) else NULL

  for (iteration in seq_len(iterations)) {
    U_perp <- lapply(U, complement)
    Z <- multi_mode_product(Xt, lapply(U, t))
    V <- lapply(seq_len(d), function(k) qr.Q(qr(t(unfold(Z, k)))))

    # Blocks: 0 = core, k = arm k.  Block b uses U_perp in mode b, U elsewhere.
    # A mode with rank equal to its dimension has no arm.
    mats_for <- function(b) lapply(seq_len(d), function(k) {
      if (k == b) U_perp[[k]] else U[[k]]
    })
    blocks <- c(0, which(dims > ranks))
    Pin <- lapply(blocks, function(b) {
      input_contraction(X, input_dims, mats_for(b)[in_modes])
    })
    Out <- lapply(blocks, function(b) {
      output_contraction(lapply(mats_for(b)[out_modes], function(M) M))
    })
    Q <- lapply(blocks, function(b) {
      if (b == 0) return(diag(prod(ranks)))
      sizes <- ranks
      sizes[b] <- dims[b] - ranks[b]
      arm_mixing(b, sizes, V[[b]])
    })

    # Normal equations: F_b = kron(Out_b, Pin_b) %*% Q_b, rows (m, j).
    sizes_param <- vapply(Q, ncol, numeric(1))
    offsets <- c(0, cumsum(sizes_param))
    G <- matrix(0, sum(sizes_param), sum(sizes_param))
    rhs <- numeric(sum(sizes_param))
    for (a in seq_along(blocks)) {
      ia <- offsets[a] + seq_len(sizes_param[a])
      rhs[ia] <- crossprod(Q[[a]],
                           as.numeric(crossprod(Pin[[a]], Y %*% Out[[a]])))
      for (b in seq_along(blocks)) {
        if (b < a) next
        ib <- offsets[b] + seq_len(sizes_param[b])
        K <- kronecker(crossprod(Out[[a]], Out[[b]]),
                       crossprod(Pin[[a]], Pin[[b]]))
        G[ia, ib] <- crossprod(Q[[a]], K %*% Q[[b]])
        if (b > a) G[ib, ia] <- t(G[ia, ib])
      }
    }
    theta <- solve(G, rhs)

    # Tangent-space point and HOSVD retraction (as in the original code).
    B <- array(theta[offsets[1] + seq_len(sizes_param[1])], ranks)
    tilde <- multi_mode_product(B, U)
    for (a in seq_along(blocks)[-1]) {
      k <- blocks[a]
      D <- matrix(theta[offsets[a] + seq_len(sizes_param[a])],
                  nrow = dims[k] - ranks[k])
      others <- setdiff(seq_len(d), k)
      kron_others <- kron_fastest_first(U[others])
      tilde <- tilde + fold(U_perp[[k]] %*% D %*% t(kron_others %*% V[[k]]),
                            k, dims)
    }
    retracted <- hosvd_truncate(tilde, ranks)
    U <- retracted$U
    Xt <- retracted$est
    if (keep_history) history[[iteration + 1]] <- Xt
  }
  list(coefficient = Xt, coefficient_history = history)
}
