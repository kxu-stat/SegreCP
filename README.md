# SegreCP

`SegreCP` is the R package for *Guaranteed Noisy CP Tensor Recovery via
Riemannian Optimization on the Segre Manifold* (NeurIPS 2026,
[arXiv:2510.00569](https://arxiv.org/abs/2510.00569)). It implements Riemannian
gradient descent (RGD) and Riemannian Gauss–Newton (RGN) on the Segre manifold
for CP tensor decomposition, scalar-on-tensor regression and tensor-on-tensor
regression, with CPCA and HOSVD initialization.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kxu-stat/SegreCP", build_vignettes = TRUE)
```

The package requires R 4.2 or later and a C++ compiler.

## Example

```r
library(SegreCP)
set.seed(1)
dims <- c(5, 6, 4)
factors <- lapply(dims, function(p) {
  z <- rnorm(p)
  z / sqrt(sum(z^2))
})
signal <- array(3 * Reduce(kronecker, rev(factors)), dim = dims)
observed <- signal + array(rnorm(prod(dims), sd = 0.05), dim = dims)

fit <- segre_decompose(observed, rank = 1, method = "rgn", init = "cpca",
                       control = segre_control(max_iter = 10))
summary(fit)
convergence_data(fit)
```

| Function | Purpose |
|---|---|
| `segre_decompose()` | CP decomposition by RGD or RGN |
| `segre_regress()` | scalar-on-tensor and tensor-on-tensor regression by RGD or RGN |
| `cp_initialize()` | CPCA, HOSVD or random initialization |
| `segre_control()` | number of iterations, step size, stopping rule |
| `reconstruct_cp()`, `coef()`, `predict()`, `convergence_data()` | fitted tensors, coefficients, predictions and objective paths |

By default (`profile = "manuscript"`), the components are updated sequentially
in a fixed order, every update uses the rank-(1, ..., 1) T-HOSVD retraction,
RGD uses the fixed step size of the paper's convergence results (0.5 for
decomposition, 1/σ<sub>X</sub><sup>2</sup> for regression), RGN is
unregularized, and all requested iterations are run. The vignettes
`getting-started`, `regression-data-layout` and `algorithm-fidelity` describe
the interface, the data layouts (base R arrays and `rTensor::Tensor` objects)
and the options.

## Reproducing the paper

The [`reproducibility/`](reproducibility/README.md) directory contains the
code for all simulations, figures and tables and for the trade-data analysis,
with instructions. It is part of the GitHub repository, not of the installed
package.

## Citation

Run `citation("SegreCP")` in R for the reference to the paper.

## License

MIT; see [`LICENSE.md`](LICENSE.md).
