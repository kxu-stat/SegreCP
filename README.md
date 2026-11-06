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

## Examples

CP decomposition of a noisy rank-two tensor:

```r
library(SegreCP)
set.seed(1)
dims <- c(8, 7, 6)
factors <- lapply(dims, function(p) qr.Q(qr(matrix(rnorm(2 * p), p))))
truth <- cp_components(c(5, 3), factors)
observed <- reconstruct_cp(truth) + array(rnorm(prod(dims), sd = 0.05), dims)

fit <- segre_decompose(observed, rank = 2, method = "rgn",
                       control = segre_control(max_iter = 10))
fit
cp_relative_error(fit, truth)
```

Regression of a scalar response on 3 x 4 matrix covariates:

```r
n <- 200
X <- matrix(rnorm(n * 12), n)     # row m is as.vector(X_m)
beta <- 2 * outer(c(1, 0, 0), c(0, 1, 0, 0))
y <- as.numeric(X %*% as.vector(beta)) + rnorm(n, sd = 0.1)

fit <- segre_regress(X, y, rank = 1, input_dims = c(3, 4), method = "rgn")
coef(fit)
predict(fit, X[1:5, ])
```

| Function | Purpose |
|---|---|
| `segre_decompose()` | CP decomposition by RGD or RGN |
| `segre_regress()` | scalar-on-tensor and tensor-on-tensor regression by RGD or RGN |
| `cp_initialize()` | CPCA, HOSVD or random initialization |
| `segre_control()` | number of iterations, step size, stopping rule |
| `cp_components()`, `reconstruct_cp()`, `cp_relative_error()` | CP components, dense tensors, errors |
| `coef()`, `predict()`, `summary()`, `convergence_data()` | estimates, predictions, objective and error paths |

By default (`profile = "manuscript"`), each iteration updates the components
one at a time in a fixed order, every update uses the rank-(1, ..., 1) T-HOSVD
retraction, RGD uses the fixed step size of the paper's convergence results
(0.5 for decomposition, 1/σ<sub>X</sub><sup>2</sup> for regression), RGN is
unregularized, and all requested iterations are run. The vignettes describe
the interface (`vignette("getting-started", package = "SegreCP")`), the data
layouts for regression (`"regression-data-layout"`) and the algorithms,
options and starts (`"algorithms"`).

## Reproducing the paper

The [`reproducibility/`](reproducibility/README.md) directory contains the
code for all simulations, figures and tables and for the trade-data analysis,
with instructions. It is part of the GitHub repository, not of the installed
package.

## Citation

Run `citation("SegreCP")` in R for the reference to the paper.

## License

MIT; see [`LICENSE.md`](LICENSE.md).
