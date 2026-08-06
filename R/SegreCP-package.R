#' SegreCP: Riemannian CP tensor estimation
#'
#' `SegreCP` implements the RGD and RGN algorithms described in
#' *Guaranteed Noisy CP Tensor Recovery via Riemannian Optimization on the
#' Segre Manifold*. Computation is performed on base R arrays and dense numeric
#' matrices; `rTensor::Tensor` inputs are converted once at the API boundary.
#'
#' @keywords internal
#' @useDynLib SegreCP, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
