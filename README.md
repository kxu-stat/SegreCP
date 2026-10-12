# SegreCP

`SegreCP` is the R package accompanying [*Guaranteed Noisy CP Tensor Recovery
via Riemannian Optimization on the Segre Manifold*](https://arxiv.org/abs/2510.00569)
by Ke Xu and Yuefeng Han. It implements the paper's Riemannian gradient descent
(RGD) and Riemannian Gauss–Newton (RGN) methods for CP tensor decomposition,
scalar-on-tensor regression, and tensor-on-tensor regression. CPCA and HOSVD
initializations are included.

## Installation

```r
install.packages("remotes")
remotes::install_github(
  "kxu-stat/SegreCP", build_vignettes = TRUE, dependencies = TRUE
)
library(SegreCP)
```

The package requires R 4.2 or later and a C++ build toolchain. R package
dependencies are listed in [`DESCRIPTION`](DESCRIPTION).

## Example

Fit a rank-one CP decomposition with CPCA initialization and RGD:

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

fit <- segre_decompose(
  observed, rank = 1, method = "rgd", init = "cpca",
  control = segre_control(max_iter = 10)
)
summary(fit)
convergence_data(fit)
```

Choose `method = "rgn"` for RGN or `init = "hosvd"` for HOSVD. The
[getting-started vignette](vignettes/getting-started.Rmd) gives complete examples
for decomposition and both regression tasks, including a convergence plot.
Base R arrays and `rTensor::Tensor` inputs are supported; see the
[data-layout vignette](vignettes/regression-data-layout.Rmd).

## Numerical experiments

The repository separates new algorithm runs from figures and tables computed
from saved experiment records:

- [`reproducibility/simulations/`](reproducibility/simulations/README.md)
  generates the documented decomposition and tensor-on-tensor regression
  designs and reruns the authors' RGD/RGN methods.
- The [initialization and convergence
  vignette](vignettes/rebuttal-initialization-convergence.Rmd) and
  [tensor-regression results
  vignette](vignettes/rebuttal-regression-reproducibility.Rmd) regenerate plots
  and summaries from the bundled replication-level records, including the
  reported comparison methods.
- [`reproducibility/real_data/`](reproducibility/real_data/README.md) provides
  official data sources, preprocessing code, and instructions for rerunning
  the authors' methods on the trade-data analysis. Third-party raw and
  processed data are not distributed here.

The [experiment guide](reproducibility/README.md) gives the run commands and
explains which figures are rebuilt from recorded results.

## Citation and contact

Please cite [Xu and Han, arXiv:2510.00569](https://arxiv.org/abs/2510.00569).
The package also provides `citation("SegreCP")` in R. For questions, contact
[Ke Xu](mailto:kxu6@nd.edu). The code is released under the
[MIT License](LICENSE.md).
