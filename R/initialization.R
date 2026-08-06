.mode_unfold <- function(x, mode) {
  dx <- dim(x)
  perm <- c(mode, setdiff(seq_along(dx), mode))
  matrix(aperm(x, perm), nrow = dx[[mode]])
}

.group_unfold <- function(x, modes) {
  dx <- dim(x)
  rest <- setdiff(seq_along(dx), modes)
  matrix(aperm(x, c(modes, rest)), nrow = prod(dx[modes]))
}

.leading_svd <- function(x, rank) {
  if (rank > min(dim(x))) {
    stop(sprintf("Requested rank %d exceeds the unfolding rank bound %d.",
                 rank, min(dim(x))), call. = FALSE)
  }
  base::svd(x, nu = rank, nv = rank)
}

.normalize_vector <- function(x, label = "vector") {
  x <- as.numeric(x)
  nrm <- sqrt(sum(x * x))
  if (!is.finite(nrm) || nrm <= .Machine$double.eps) {
    stop(sprintf("Cannot normalize a zero %s.", label), call. = FALSE)
  }
  x / nrm
}

.auto_cpca_split <- function(dims) {
  d <- length(dims)
  if (d < 2L) stop("CPCA requires a tensor with at least two modes.")
  if (d > 20L) {
    stop("Automatic CPCA split enumeration is limited to 20 modes; supply `split`.",
         call. = FALSE)
  }
  # Algorithm CPCA in the manuscript minimizes the *raw product-scale*
  # imbalance |prod(p[S]) - sqrt(prod(p))|.  This is deliberately not replaced
  # by a log-scale imbalance, because the two order candidate subsets
  # differently when mode dimensions are unequal.
  target <- sqrt(prod(as.numeric(dims)))
  candidates <- unlist(lapply(seq_len(d - 1L), function(k) {
    z <- utils::combn(seq_len(d), k, simplify = FALSE)
    z
  }), recursive = FALSE)
  scores <- vapply(candidates, function(S) {
    abs(prod(as.numeric(dims[S])) - target)
  }, numeric(1L))
  candidates[[which.min(scores)]]
}

.validate_cpca_split <- function(split, d) {
  if (!is.numeric(split) || any(!is.finite(split)) ||
      any(split != floor(split))) {
    stop("`split` must contain whole-number mode indices.", call. = FALSE)
  }
  split <- sort(unique(as.integer(split)))
  if (!length(split) || length(split) >= d ||
      any(split < 1L | split > d)) {
    stop("`split` must be a non-empty proper subset of the tensor modes.",
         call. = FALSE)
  }
  split
}

.refold_rank1_factors <- function(x, dims) {
  x <- as.numeric(x)
  if (length(dims) == 1L) return(list(.normalize_vector(x)))
  arr <- array(x, dim = dims)
  U <- lapply(seq_along(dims), function(k) {
    .normalize_vector(base::svd(.mode_unfold(arr, k), nu = 1L, nv = 0L)$u[, 1L])
  })
  product <- .rank1_vector(1, U)
  if (sum(product * x) < 0) U[[1L]] <- -U[[1L]]
  U
}

.cpca_initialize <- function(x, rank, split = NULL) {
  dims <- as.integer(dim(x))
  d <- length(dims)
  if (is.null(split)) split <- .auto_cpca_split(dims)
  split <- .validate_cpca_split(split, d)
  complement <- setdiff(seq_len(d), split)
  unfolded <- .group_unfold(x, split)
  s <- .leading_svd(unfolded, rank)

  factors <- vector("list", d)
  left <- lapply(seq_len(rank), function(i) {
    .refold_rank1_factors(s$u[, i], dims[split])
  })
  right <- lapply(seq_len(rank), function(i) {
    .refold_rank1_factors(s$v[, i], dims[complement])
  })
  for (pos in seq_along(split)) {
    factors[[split[[pos]]]] <- do.call(cbind, lapply(left, `[[`, pos))
  }
  for (pos in seq_along(complement)) {
    factors[[complement[[pos]]]] <- do.call(cbind, lapply(right, `[[`, pos))
  }

  xvec <- as.numeric(x)
  comps <- lapply(seq_len(rank), function(i) {
    Ui <- lapply(factors, function(A) A[, i])
    lambda <- sum(xvec * .rank1_vector(1, Ui))
    .canonicalize_component(list(lambda = lambda, U = Ui))
  })
  list(comps = .validate_components(comps, dims, rank), split = split,
       singular_values = s$d[seq_len(rank)])
}

.hosvd_initialize <- function(x, rank) {
  dims <- as.integer(dim(x))
  U <- lapply(seq_along(dims), function(k) {
    .leading_svd(.mode_unfold(x, k), rank)$u[, seq_len(rank), drop = FALSE]
  })
  xvec <- as.numeric(x)
  comps <- lapply(seq_len(rank), function(i) {
    Ui <- lapply(U, function(A) A[, i])
    lambda <- sum(xvec * .rank1_vector(1, Ui))
    .canonicalize_component(list(lambda = lambda, U = Ui))
  })
  list(comps = .validate_components(comps, dims, rank))
}

.with_local_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  seed <- .as_whole_number(seed, "seed", allow_zero = TRUE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.random_initialize <- function(x, rank, seed = NULL, project_weights = TRUE) {
  dims <- as.integer(dim(x))
  comps <- .with_local_seed(seed, lapply(seq_len(rank), function(i) {
    Ui <- lapply(dims, function(p) .normalize_vector(stats::rnorm(p)))
    lambda <- if (project_weights) {
      sum(as.numeric(x) * .rank1_vector(1, Ui))
    } else 1
    .canonicalize_component(list(lambda = lambda, U = Ui))
  }))
  list(comps = .validate_components(comps, dims, rank))
}

#' Initialize a CP tensor
#'
#' The CPCA implementation follows the manuscript algorithm: it uses a
#' non-empty proper mode split and, when the split is omitted, minimizes
#' `abs(prod(dims[split]) - sqrt(prod(dims)))`. It deterministically resolves
#' ties, aligns the two refolded singular-vector sides, re-estimates each
#' component weight by projection, and absorbs a negative weight into the first
#' factor. The HOSVD initializer uses the diagonal of a multilinear HOSVD core;
#' it is an initializer, not a general CP decomposition algorithm.
#'
#' @param x A numeric array/matrix or `rTensor::Tensor` to initialize.
#' @param rank Positive target CP rank.
#' @param method One of `"cpca"`, `"hosvd"`, or `"random"`. A function or an
#'   already constructed component list is also accepted. A function is called
#'   as `method(x, rank, ...)` and must return components or
#'   `list(comps = components)`.
#' @param split Optional CPCA mode subset, using one-based R indices.
#' @param seed Optional local seed for random initialization. The caller's RNG
#'   state is restored.
#' @param project_random_weights If `TRUE`, estimate random-initialization
#'   weights by projection onto `x`; otherwise use unit weights.
#' @param ... Additional arguments passed to a custom initialization function.
#' @return A list with `comps`, `method`, and method-specific details.
#' @examples
#' set.seed(1)
#' x <- array(rnorm(4 * 5 * 6), c(4, 5, 6))
#' init_cpca <- cp_initialize(x, rank = 2, method = "cpca")
#' init_hosvd <- cp_initialize(x, rank = 2, method = "hosvd")
#' @export
cp_initialize <- function(x, rank, method = c("cpca", "hosvd", "random"),
                          split = NULL, seed = NULL,
                          project_random_weights = TRUE, ...) {
  x <- as_segre_array(x)
  rank <- .as_whole_number(rank, "rank")
  dims <- as.integer(dim(x))

  if (is.function(method)) {
    ans <- method(x, rank, ...)
    comps <- .validate_components(ans, dims, rank)
    return(list(comps = comps, method = "custom_function"))
  }
  if (is.list(method)) {
    comps <- .validate_components(method, dims, rank)
    return(list(comps = comps, method = "custom_components"))
  }
  method <- match.arg(method)
  ans <- switch(
    method,
    cpca = .cpca_initialize(x, rank, split),
    hosvd = .hosvd_initialize(x, rank),
    random = .random_initialize(x, rank, seed, project_random_weights)
  )
  ans$method <- method
  ans
}

.resolve_initialization <- function(moment, rank, init, split = NULL,
                                    seed = NULL, init_args = list()) {
  if (!is.list(init_args)) stop("`init_args` must be a list.", call. = FALSE)
  args <- c(list(x = moment, rank = rank, method = init, split = split,
                 seed = seed), init_args)
  do.call(cp_initialize, args)
}
