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
  if (d < 2L) {
    stop("CPCA requires a tensor with at least two modes.", call. = FALSE)
  }
  if (d > 20L) {
    stop("Automatic CPCA split enumeration is limited to 20 modes; supply `split`.",
         call. = FALSE)
  }
  # The criterion of the CPCA algorithm in the paper, |prod(p[S]) -
  # sqrt(prod(p))|, on the original scale. (On the log scale, a different
  # subset can be chosen when the dimensions differ.) Ties go to the first
  # subset in the order of utils::combn().
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

#' Initial CP components
#'
#' Computes `rank` rank-one components from a tensor, as a start for
#' [segre_decompose()] (the observed tensor) or [segre_regress()] (the moment
#' of the design and the response).
#'
#' * `"cpca"`, composite PCA: the tensor is unfolded into a matrix whose rows
#'   are indexed by the modes in `split`. Each of the `rank` leading pairs of
#'   singular vectors is refolded into two tensors, and each of these is
#'   replaced by the rank-one tensor formed from the leading left singular
#'   vectors of its unfoldings. By default `split` is the set of modes whose
#'   dimensions have the product closest to `sqrt(prod(dim(x)))`.
#' * `"hosvd"`: the factors of component `i` are the `i`-th leading left
#'   singular vectors of the mode unfoldings of `x`.
#' * `"random"`: factors with independent standard normal entries, scaled to
#'   unit length.
#'
#' The weight of each component is the inner product of `x` with the rank-one
#' tensor of its unit factors (or 1 for random factors with
#' `project_random_weights = FALSE`). A negative weight is made positive by
#' changing the sign of the first factor.
#'
#' @param x A numeric array or an `rTensor::Tensor`.
#' @param rank Number of components.
#' @param method `"cpca"`, `"hosvd"` or `"random"`; or a function called as
#'   `method(x, rank, ...)` that returns a component list; or a component list
#'   (see [cp_components()]), which is checked and returned.
#' @param split For `"cpca"`, the modes that index the rows of the unfolding,
#'   a non-empty proper subset of `seq_along(dim(x))`.
#' @param seed Seed for `"random"`. The random number stream of the session
#'   is left unchanged.
#' @param project_random_weights For `"random"`, whether the weights are
#'   computed from `x` (`TRUE`) or set to 1.
#' @param ... Further arguments for a function `method`.
#' @return A list with `comps`, the list of components `list(lambda, U)` with
#'   unit factors, and `method`. For `"cpca"` also `split` and the leading
#'   `singular_values` of the unfolding.
#' @examples
#' set.seed(1)
#' x <- array(rnorm(4 * 5 * 6), c(4, 5, 6))
#' init <- cp_initialize(x, rank = 2, method = "cpca")
#' init$split
#' vapply(init$comps, `[[`, numeric(1), "lambda")
#' cp_initialize(x, rank = 2, method = "random", seed = 1)$comps[[1]]$lambda
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

# Start of a fit from a supplied component list.
.component_start <- function(init, dims, rank, init_args) {
  if (!is.list(init_args)) stop("`init_args` must be a list.", call. = FALSE)
  if (length(init_args)) {
    stop("`init_args` is used only when `init` is a method name or a ",
         "function.", call. = FALSE)
  }
  list(comps = .validate_components(init, dims, rank),
       method = "custom_components")
}

# Start of a fit from a method name or a function, computed from `moment`
# (the observed tensor, or the moment of the design and the response).
.resolve_initialization <- function(moment, rank, init, split = NULL,
                                    seed = NULL, init_args = list()) {
  if (!is.list(init_args)) stop("`init_args` must be a list.", call. = FALSE)
  args <- c(list(x = moment, rank = rank, method = init, split = split,
                 seed = seed), init_args)
  do.call(cp_initialize, args)
}
