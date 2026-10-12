# SegreCP

`SegreCP` implements the Riemannian gradient descent (RGD) and Riemannian
Gauss–Newton (RGN) methods from *Guaranteed Noisy CP Tensor Recovery via
Riemannian Optimization on the Segre Manifold*. It fits low-rank CP tensor
decompositions, scalar-on-tensor regressions, and tensor-on-tensor regressions.
CPCA and HOSVD initializations are included; random and user-supplied starts
are also supported. The installed estimators are **RGD and RGN only**.

## Installation

From the `SegreCP` source directory (the directory containing `DESCRIPTION`):

```sh
R CMD INSTALL .
```

You need R >= 4.2, a C++ build toolchain, and the R packages `Rcpp` and
`rTensor`. The repository URL for a GitHub installation command has not yet
been set. Before public release, the copyright holders must also replace the
current all-rights-reserved research-review license with a chosen public
license; the present source should not be described as open source.

## Quick start

The examples below use small simulated data and run without downloading the
paper's datasets. `rank` is the number of CP components; `method` selects RGD
or RGN, and `init` selects CPCA or HOSVD. The default `segre_control()` profile
follows the displayed paper algorithms. For applications, see
`?segre_control` for explicitly recorded numerical safeguards in the
`"practical"` profile.

### CP tensor decomposition

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

fit_decomp <- segre_decompose(
  observed, rank = 1, method = "rgd", init = "cpca",
  control = segre_control(max_iter = 10)
)
dim(reconstruct_cp(fit_decomp))
summary(fit_decomp)
```

To inspect the observed objective after each sweep, use
`convergence_data(fit_decomp)`. For simulated data, set `keep_history = TRUE`
when fitting and pass the known CP truth to obtain reconstruction error as
well. For decomposition with an identity observation operator, an RGN
component update equals a unit-step RGD update.

### Scalar-on-tensor regression

Each row of `X` below is a vectorized `3 x 4` predictor tensor. A numeric
array with dimensions `c(n, 3, 4)` is accepted too.

```r
set.seed(2)
n <- 80
input_dims <- c(3, 4)
X <- matrix(rnorm(n * prod(input_dims)), nrow = n)
beta <- 2 * outer(c(1, 0, 0), c(1, 0, 0, 0))
y <- as.numeric(X %*% as.vector(beta) + rnorm(n, sd = 0.05))

fit_scalar <- segre_regress(
  X, y, rank = 1, input_dims = input_dims,
  method = "rgn", init = "hosvd",
  control = segre_control(max_iter = 5)
)
dim(coef(fit_scalar))
predict(fit_scalar, X[1:3, , drop = FALSE])
```

### Tensor-on-tensor regression

The coefficient modes are ordered as input modes followed by output modes;
here the response has dimensions `c(n, 2, 2)`.

```r
set.seed(3)
n <- 60
input_dims <- c(2, 3)
output_dims <- c(2, 2)
X <- matrix(rnorm(n * prod(input_dims)), nrow = n)
beta_factors <- list(c(1, 0), c(1, 0, 0), c(1, 0), c(1, 0))
beta <- array(2 * Reduce(kronecker, rev(beta_factors)),
              dim = c(input_dims, output_dims))
Y <- array(as.numeric(X %*% matrix(beta, nrow = prod(input_dims))) +
             rnorm(n * prod(output_dims), sd = 0.05),
           dim = c(n, output_dims))

fit_tensor <- segre_regress(
  X, Y, rank = 1, input_dims = input_dims,
  output_dims = output_dims, method = "rgd", init = "cpca",
  control = segre_control(max_iter = 5)
)
dim(coef(fit_tensor))
dim(predict(fit_tensor, X[1:3, , drop = FALSE]))
```

The same calls can use `method = "rgd"` or `"rgn"` and `init = "cpca"` or
`"hosvd"`. `cp_initialize()` also accepts random starts, a list of CP
components, or a custom initialization function. See `?segre_decompose`,
`?segre_regress`, and `?cp_initialize` for the full arguments and return
values. The [getting-started vignette](vignettes/getting-started.Rmd) runs
these three tasks and plots an observed-objective convergence trace.

## Tensor inputs

The computational core uses base R numeric arrays. For regression, a numeric
matrix with one vectorized predictor tensor per row is the most efficient
representation. The package also accepts arrays with dimensions
`c(n, input_dims)`, lists of individual tensors, and `rTensor::Tensor`
objects, converting them at the API boundary. Responses can be a numeric
vector, matrix, array, `rTensor::Tensor`, or list of tensors as appropriate.
See the [data-layout vignette](vignettes/regression-data-layout.Rmd) for a
worked matrix/array/`rTensor` comparison.

## Reproduce the paper's numerical results

The [tensor-regression simulation
vignette](vignettes/rebuttal-regression-reproducibility.Rmd) reconstructs
tables and plots from bundled replication-level records. The
[initialization and convergence
vignette](vignettes/rebuttal-initialization-convergence.Rmd) visualizes the
CPCA/HOSVD comparison and optimization trajectories. These vignettes
recompute summaries from recorded runs; building them does **not** rerun the
full, expensive multi-method simulations. RGD/RGN and initialization are the
only supported estimators. Comparison code belongs to the source checkout's
[`reproducibility/` directory](reproducibility/README.md), not to the
installed package.

For a new simulation rather than a redraw of stored results, the source
checkout also includes [fresh RGD/RGN decomposition and tensor-regression
runners](reproducibility/simulations/README.md). They generate data under the
rebuttal-stage designs and save fit metrics, convergence traces, random seeds,
and data fingerprints outside the package directory. They do **not** rerun the
other five estimators or reproduce hardware-dependent timing values.

The [real-trade vignette](vignettes/real-trade-reproducibility.Rmd) and
[`reproducibility/real_data/` guide](reproducibility/real_data/README.md)
explain how to obtain official data into a **private directory**, rebuild and
check the 2018–2019 analysis inputs, and rerun rank-three RGD/RGN through the
public package API. Neither raw nor processed third-party trade data are
distributed. The verified offline source-to-input rebuild and conditional
RGD/RGN refit are documented in [`inst/VALIDATION.md`](inst/VALIDATION.md).
A fresh download from every provider and a new fit of the other five methods
have not been verified. Reproduction scripts require a source checkout;
`R CMD INSTALL` does not install `reproducibility/`.

## Citation and tests

```r
citation("SegreCP")
```

With `testthat` installed, run the package tests from the `SegreCP` source
directory using `testthat::test_local(".")`. The detailed algorithm-to-paper
mapping and historical compatibility notes are in
[`inst/ALGORITHM_FIDELITY.md`](inst/ALGORITHM_FIDELITY.md) and
[`inst/RELEASE_PROVENANCE.md`](inst/RELEASE_PROVENANCE.md).
