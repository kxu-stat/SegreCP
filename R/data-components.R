# Input checks, conversion of tensors, designs and responses, and CP
# components.

# TRUE for an rTensor::Tensor, without attaching rTensor.
.is_rtensor <- function(x) {
  inherits(x, "Tensor") && methods::is(x, "Tensor")
}

.as_whole_number <- function(x, name, allow_zero = FALSE) {
  lower <- if (allow_zero) 0 else 1
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x != floor(x) || x > .Machine$integer.max) {
    qualifier <- if (allow_zero) "nonnegative" else "positive"
    stop(sprintf("`%s` must be one %s integer.", name, qualifier),
         call. = FALSE)
  }
  as.integer(x)
}

.as_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.as_mode_dims <- function(x, name, allow_empty = FALSE) {
  if (is.null(x)) return(NULL)
  if (!is.numeric(x) || any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < 1) || any(x > .Machine$integer.max) ||
      (!allow_empty && !length(x))) {
    stop(sprintf("`%s` must contain positive integers%s.", name,
                 if (allow_empty) " (or be empty for a scalar)" else ""),
         call. = FALSE)
  }
  as.integer(x)
}

#' Convert a tensor to a base R array
#'
#' `SegreCP` computes with base R arrays. The fitting functions call
#' `as_segre_array()` on their input, so an `rTensor::Tensor` can be passed
#' wherever an array is expected.
#'
#' @param x A numeric array or matrix, or an `rTensor::Tensor`.
#' @return A numeric (double) array with the dimensions and entries of `x`.
#' @seealso [as_rtensor()] for the conversion in the other direction.
#' @examples
#' x <- rTensor::as.tensor(array(1:24, c(2, 3, 4)))
#' dim(as_segre_array(x))
#' @export
as_segre_array <- function(x) {
  if (.is_rtensor(x)) {
    x <- methods::slot(x, "data")
  }
  if (!is.numeric(x) || is.null(dim(x))) {
    stop("`x` must be a numeric base R array/matrix or an rTensor::Tensor.",
         call. = FALSE)
  }
  if (!length(dim(x)) || any(dim(x) < 1L)) {
    stop("Tensor dimensions must all be positive.", call. = FALSE)
  }
  out <- array(as.numeric(x), dim = as.integer(dim(x)))
  if (any(!is.finite(out))) {
    stop("Tensor entries must all be finite.", call. = FALSE)
  }
  out
}

.format_tensor <- function(x, format = c("array", "rTensor")) {
  format <- match.arg(format)
  if (format == "array") return(x)
  rTensor::as.tensor(x)
}

.rank1_vector <- function(lambda, U) {
  as.numeric(make_rank1_vec_cpp(as.numeric(lambda), lapply(U, as.numeric)))
}

.canonicalize_component <- function(component) {
  component$lambda <- as.numeric(component$lambda)
  component$U <- lapply(component$U, as.numeric)
  if (component$lambda < 0) {
    component$lambda <- -component$lambda
    component$U[[1L]] <- -component$U[[1L]]
  }
  component
}

# Checks a component list (or list(comps = ...)) and returns it with unit
# factors and positive weights.
.validate_components <- function(components, dims = NULL, rank = NULL,
                                 normalize = TRUE) {
  if (is.list(components) && !is.null(components$comps)) {
    components <- components$comps
  }
  if (!is.list(components) || !length(components)) {
    stop("CP components must be a non-empty list of list(lambda, U).",
         call. = FALSE)
  }
  if (!is.null(rank) && length(components) != rank) {
    stop(sprintf("Expected %d CP components (the rank), got %d.",
                 rank, length(components)), call. = FALSE)
  }

  inferred <- NULL
  out <- vector("list", length(components))
  for (i in seq_along(components)) {
    comp <- components[[i]]
    if (!is.list(comp) || length(comp$lambda) != 1L || !is.list(comp$U) ||
        !length(comp$U)) {
      stop(sprintf(
        "Component %d must be list(lambda = scalar, U = list of vectors).", i
      ), call. = FALSE)
    }
    Ui <- lapply(comp$U, function(u) {
      u <- as.numeric(u)
      if (!length(u) || any(!is.finite(u))) {
        stop("Every factor must be a non-empty finite numeric vector.",
             call. = FALSE)
      }
      u
    })
    dims_i <- vapply(Ui, length, integer(1L))
    if (is.null(inferred)) inferred <- dims_i
    if (!identical(dims_i, inferred)) {
      stop("All CP components must have the same mode dimensions.",
           call. = FALSE)
    }
    lambda <- as.numeric(comp$lambda)
    if (!is.finite(lambda)) {
      stop("All component weights must be finite.", call. = FALSE)
    }
    if (normalize) {
      for (k in seq_along(Ui)) {
        nrm <- sqrt(sum(Ui[[k]] * Ui[[k]]))
        if (!is.finite(nrm) || nrm <= .Machine$double.eps) {
          stop(sprintf("Component %d, mode %d has zero norm.", i, k),
               call. = FALSE)
        }
        Ui[[k]] <- Ui[[k]] / nrm
        lambda <- lambda * nrm
      }
    }
    out[[i]] <- .canonicalize_component(list(lambda = lambda, U = Ui))
    if (out[[i]]$lambda == 0) {
      stop(sprintf("Component %d has zero weight and is not a Segre point.", i),
           call. = FALSE)
    }
  }
  if (!is.null(dims) && !identical(as.integer(dims), as.integer(inferred))) {
    stop(sprintf("CP components have dimensions (%s); expected (%s).",
                 paste(inferred, collapse = ", "),
                 paste(dims, collapse = ", ")), call. = FALSE)
  }
  out
}

#' CP components from weights and factor matrices
#'
#' `SegreCP` represents a CP tensor
#' \eqn{\sum_i \lambda_i\, u_{i1} \otimes \cdots \otimes u_{id}}{sum_i lambda_i u_i1 o ... o u_id}
#' as a list of components, component `i` being
#' `list(lambda = lambda_i, U = list(u_i1, ..., u_id))`. `cp_components()`
#' builds this list from factor matrices. The factors are scaled to unit length
#' (the scales are absorbed into the weights) and a negative weight is made
#' positive by changing the sign of the first factor.
#'
#' @param lambda Weights, one per component.
#' @param factors A list with one factor matrix per mode. Column `i` of
#'   `factors[[k]]` is the mode-`k` factor of component `i`.
#' @return A list of components `list(lambda, U)`.
#' @examples
#' A <- matrix(c(1, 0, 0, 1, 1, 0), 3)  # 3 x 2
#' B <- diag(2)
#' comps <- cp_components(c(2, 1), list(A, B))
#' comps[[2]]
#' reconstruct_cp(comps)
#' @export
cp_components <- function(lambda, factors) {
  lambda <- as.numeric(lambda)
  if (!length(lambda) || any(!is.finite(lambda))) {
    stop("`lambda` must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (!is.list(factors) || !length(factors)) {
    stop("`factors` must be a non-empty list of factor matrices.",
         call. = FALSE)
  }
  rank <- length(lambda)
  factors <- lapply(seq_along(factors), function(k) {
    A <- as.matrix(factors[[k]])
    if (!is.numeric(A) || ncol(A) != rank || any(!is.finite(A))) {
      stop(sprintf("Factor matrix %d must be finite and have %d columns.",
                   k, rank), call. = FALSE)
    }
    A
  })
  comps <- lapply(seq_len(rank), function(i) {
    list(lambda = lambda[[i]], U = lapply(factors, function(A) A[, i]))
  })
  .validate_components(comps, rank = rank)
}

#' Dense tensor from CP components
#'
#' @param components A list of components `list(lambda, U)` (see
#'   [cp_components()]) or a fit from [segre_decompose()] or [segre_regress()].
#' @param format `"array"` or `"rTensor"`.
#' @return The tensor \eqn{\sum_i \lambda_i\, u_{i1} \otimes \cdots \otimes u_{id}}{sum_i lambda_i u_i1 o ... o u_id}.
#' @examples
#' comps <- cp_components(2, list(c(1, 0), c(0, 1, 0)))
#' reconstruct_cp(comps)
#' @export
reconstruct_cp <- function(components, format = c("array", "rTensor")) {
  if (inherits(components, "segre_cp_fit")) components <- components$components
  components <- .validate_components(components)
  dims <- vapply(components[[1L]]$U, length, integer(1L))
  value <- Reduce(`+`, lapply(components, function(comp) {
    .rank1_vector(comp$lambda, comp$U)
  }))
  .format_tensor(array(value, dim = dims), match.arg(format))
}

# <sum_i x_i, sum_j y_j> for two component lists, from the factor Gram
# matrices.
.cp_inner_product <- function(x, y) {
  lx <- vapply(x, `[[`, numeric(1L), "lambda")
  ly <- vapply(y, `[[`, numeric(1L), "lambda")
  d <- length(x[[1L]]$U)
  gram <- matrix(1, length(x), length(y))
  for (k in seq_len(d)) {
    Ax <- do.call(cbind, lapply(x, function(z) z$U[[k]]))
    Ay <- do.call(cbind, lapply(y, function(z) z$U[[k]]))
    gram <- gram * crossprod(Ax, Ay)
  }
  as.numeric(crossprod(lx, gram %*% ly))
}

#' Relative error between two CP tensors
#'
#' Computes \eqn{\|\hat T - T\|_F / \|T\|_F}{||estimate - truth||_F / ||truth||_F}
#' from the inner products of the factors, without forming the tensors.
#'
#' @param estimate,truth Lists of CP components (see [cp_components()]) or
#'   fits from [segre_decompose()] or [segre_regress()], with the same mode
#'   dimensions.
#' @return The relative Frobenius error.
#' @examples
#' truth <- cp_components(c(2, 1), list(diag(3)[, 1:2], diag(4)[, 1:2]))
#' estimate <- cp_components(c(2.1, 0.9), list(diag(3)[, 1:2], diag(4)[, 1:2]))
#' cp_relative_error(estimate, truth)
#' @export
cp_relative_error <- function(estimate, truth) {
  if (inherits(estimate, "segre_cp_fit")) estimate <- estimate$components
  if (inherits(truth, "segre_cp_fit")) truth <- truth$components
  estimate <- .validate_components(estimate)
  dims <- vapply(estimate[[1L]]$U, length, integer(1L))
  truth <- .validate_components(truth, dims = dims)
  ee <- .cp_inner_product(estimate, estimate)
  tt <- .cp_inner_product(truth, truth)
  et <- .cp_inner_product(estimate, truth)
  if (!is.finite(tt) || tt <= 0) {
    stop("`truth` has zero Frobenius norm.", call. = FALSE)
  }
  sqrt(max(0, ee + tt - 2 * et)) / sqrt(tt)
}

# Design as an n x prod(input_dims) matrix with rows vec(X_m).
.as_design_matrix <- function(X, input_dims = NULL) {
  if (.is_rtensor(X)) X <- methods::slot(X, "data")

  if (is.list(X) && !is.data.frame(X)) {
    if (!length(X)) {
      stop("`X` must contain at least one observation.", call. = FALSE)
    }
    first <- as_segre_array(X[[1L]])
    inferred <- as.integer(dim(first))
    if (is.null(input_dims)) input_dims <- inferred
    input_dims <- .as_mode_dims(input_dims, "input_dims")
    if (!identical(input_dims, inferred)) {
      stop("`input_dims` does not match the tensors in `X`.", call. = FALSE)
    }
    rows <- lapply(seq_along(X), function(i) {
      Xi <- as_segre_array(X[[i]])
      if (!identical(as.integer(dim(Xi)), as.integer(input_dims))) {
        stop(sprintf("X[[%d]] has inconsistent dimensions.", i),
             call. = FALSE)
      }
      as.numeric(Xi)
    })
    return(list(X = do.call(rbind, rows), input_dims = as.integer(input_dims),
                n = length(rows), source = "list"))
  }

  if (!is.numeric(X) || is.null(dim(X))) {
    stop("`X` must be a numeric matrix/array, rTensor Tensor, or list of tensors.",
         call. = FALSE)
  }
  dx <- as.integer(dim(X))
  if (!length(dx) || dx[[1L]] < 1L || any(dx[-1L] < 1L)) {
    stop("`X` must contain at least one observation and positive tensor modes.",
         call. = FALSE)
  }
  if (length(dx) == 2L) {
    Xmat <- unname(as.matrix(X))
    if (is.null(input_dims)) input_dims <- ncol(Xmat)
    input_dims <- .as_mode_dims(input_dims, "input_dims")
    if (prod(input_dims) != ncol(Xmat)) {
      stop("`prod(input_dims)` must equal `ncol(X)`.", call. = FALSE)
    }
    if (any(!is.finite(Xmat))) {
      stop("Design entries must all be finite.", call. = FALSE)
    }
    return(list(X = Xmat, input_dims = input_dims,
                n = nrow(Xmat), source = "matrix"))
  }

  inferred <- dx[-1L]
  if (is.null(input_dims)) input_dims <- inferred
  input_dims <- .as_mode_dims(input_dims, "input_dims")
  if (!identical(input_dims, inferred)) {
    stop("For an array `X`, dimensions must be `c(n, input_dims)`.",
         call. = FALSE)
  }
  n <- dx[[1L]]
  Xmat <- matrix(as.numeric(X), nrow = n)
  if (any(!is.finite(Xmat))) {
    stop("Design entries must all be finite.", call. = FALSE)
  }
  list(X = Xmat, input_dims = input_dims, n = n, source = "array")
}

# Response as an n x prod(output_dims) matrix with rows vec(Y_m).
.as_response_block <- function(Y, n, output_dims = NULL) {
  if (.is_rtensor(Y)) Y <- methods::slot(Y, "data")

  if (is.list(Y) && !is.data.frame(Y)) {
    if (length(Y) != n) {
      stop("`Y` and `X` contain different sample sizes.", call. = FALSE)
    }
    first <- as_segre_array(Y[[1L]])
    inferred <- as.integer(dim(first))
    if (is.null(output_dims)) output_dims <- inferred
    output_dims <- .as_mode_dims(output_dims, "output_dims")
    q <- prod(output_dims)
    rows <- lapply(seq_along(Y), function(i) {
      Yi <- as_segre_array(Y[[i]])
      if (!identical(as.integer(dim(Yi)), as.integer(output_dims))) {
        stop(sprintf("Y[[%d]] has inconsistent dimensions.", i),
             call. = FALSE)
      }
      as.numeric(Yi)
    })
    return(list(Y = do.call(rbind, rows), output_dims = as.integer(output_dims),
                q = q, source = "list"))
  }

  if (is.numeric(Y) && is.null(dim(Y))) {
    if (length(Y) != n) {
      stop("A vector `Y` must have length `nrow(X)`.", call. = FALSE)
    }
    output_dims <- .as_mode_dims(output_dims, "output_dims", allow_empty = TRUE)
    if (length(output_dims) && prod(output_dims) != 1L) {
      stop("A vector `Y` represents a scalar response.", call. = FALSE)
    }
    if (any(!is.finite(Y))) {
      stop("Response entries must all be finite.", call. = FALSE)
    }
    return(list(Y = matrix(as.numeric(Y), ncol = 1L),
                output_dims = integer(), q = 1L, source = "vector"))
  }

  if (!is.numeric(Y) || is.null(dim(Y))) {
    stop("`Y` must be a vector, matrix/array, rTensor Tensor, or list.",
         call. = FALSE)
  }
  dy <- as.integer(dim(Y))
  if (!length(dy) || any(dy < 1L)) {
    stop("`Y` must contain at least one observation and positive tensor modes.",
         call. = FALSE)
  }
  if (dy[[1L]] != n) {
    stop("The first dimension of `Y` must equal `nrow(X)`.", call. = FALSE)
  }
  if (length(dy) == 2L) {
    Ymat <- unname(as.matrix(Y))
    if (ncol(Ymat) == 1L && (is.null(output_dims) || !length(output_dims))) {
      output_dims <- integer()
    } else if (is.null(output_dims)) {
      output_dims <- ncol(Ymat)
    }
    output_dims <- .as_mode_dims(output_dims, "output_dims", allow_empty = TRUE)
  } else {
    inferred <- dy[-1L]
    if (is.null(output_dims)) output_dims <- inferred
    output_dims <- .as_mode_dims(output_dims, "output_dims")
    if (!identical(output_dims, inferred)) {
      stop("For an array `Y`, dimensions must be `c(n, output_dims)`.",
           call. = FALSE)
    }
    Ymat <- matrix(as.numeric(Y), nrow = n)
  }
  if (prod(output_dims) != ncol(Ymat)) {
    stop("`prod(output_dims)` must equal the response size per observation.",
         call. = FALSE)
  }
  if (any(!is.finite(Ymat))) {
    stop("Response entries must all be finite.", call. = FALSE)
  }
  list(Y = Ymat, output_dims = output_dims, q = ncol(Ymat),
       source = if (length(dy) == 2L) "matrix" else "array")
}
