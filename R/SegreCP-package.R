#' SegreCP: Riemannian CP tensor estimation on the Segre manifold
#'
#' Riemannian gradient descent (RGD) and Riemannian Gauss--Newton (RGN) for
#' noisy low-CP-rank tensor estimation, the algorithms of *Guaranteed Noisy CP
#' Tensor Recovery via Riemannian Optimization on the Segre Manifold* (NeurIPS
#' 2026, <https://arxiv.org/abs/2510.00569>). The CP components are updated
#' one at a time on the Segre manifold of rank-one tensors.
#' `citation("SegreCP")` gives the reference.
#'
#' @section Main functions:
#' * [segre_decompose()]: CP decomposition of a noisy tensor.
#' * [segre_regress()]: scalar-on-tensor and tensor-on-tensor regression with
#'   a low-CP-rank coefficient.
#' * [segre_control()]: number of iterations, stopping rule and step size.
#' * [cp_initialize()]: CPCA, HOSVD and random starts.
#' * [cp_components()], [reconstruct_cp()], [cp_relative_error()]: CP
#'   components, dense tensors and errors.
#' * [convergence_data()]: objective and error at each iteration.
#'
#' @section Data layout:
#' Tensors are base R arrays; `rTensor::Tensor` objects are also accepted. A
#' CP tensor is a list of components `list(lambda, U)`, each with a weight
#' `lambda` and a list `U` of one factor vector per mode. In regression, the
#' rows of the design matrix are the vectorized covariate tensors.
#'
#' @keywords internal
#' @useDynLib SegreCP, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
