# SegreCP

`SegreCP` is an R package for the sequential Riemannian gradient descent
(RGD) and Riemannian Gauss--Newton (RGN) algorithms in *Guaranteed Noisy CP
Tensor Recovery via Riemannian Optimization on the Segre Manifold*. It covers

- noisy CP tensor decomposition (identity observation operator),
- scalar-on-tensor regression, and
- tensor-on-tensor regression.

Every component iterate remains a nonzero rank-one tensor. Updates use the
exact rank-`(1, ..., 1)` T-HOSVD retraction specialized to a rank-one-plus-
tangent proposal.

## Supported software scope

The supported estimators are **RGD and RGN only**. The package also provides the
initialization layer required by those estimators: CPCA, HOSVD, random starts,
and user-supplied initializations. It does not install or export CP-ALS,
Orth-ALS, ICO, CP-RRR, Tucker-RGD, TT-IHT, HOPLS, or another competing
optimizer.

For manuscript experiments, the single authors-method entry point is installed
at `inst/scripts/ours_methods.R`. The repository-only comparison entry point is
`reproducibility/comparison_methods.R`; `.Rbuildignore` excludes that entire
directory from the built package. The table validator under `inst/reproduce/`
only reads frozen numerical records and contains no comparator implementation.

## Installation

When this folder is published as its own GitHub repository:

```r
# install.packages("remotes")
remotes::install_github("OWNER/SegreCP")
```

When it remains under `software/SegreCP` in the manuscript repository:

```r
remotes::install_github(
  "OWNER/REPOSITORY",
  subdir = "software/SegreCP"
)
```

For a local checkout, run from a shell at the manuscript-repository root:

```sh
R CMD INSTALL software/SegreCP
```

The package currently carries a research-review, all-rights-reserved license
placeholder. The copyright holders should select an open-source license before
advertising the package as open-source software.

## Data formats and efficiency

The computational core uses base R numeric arrays and matrices. They have
contiguous column-major storage and can be passed directly to the compiled C++
kernels. An `rTensor::Tensor` is fully supported and converted once at the API
boundary, preserving its dimensions and entry order.

For decomposition, use either a base array with dimensions
`c(p1, ..., pd)` or an `rTensor::Tensor` of the same size.

For regression, the most efficient design representation is an
`n`-by-`prod(input_dims)` numeric matrix. The following are equivalent:

- a matrix whose rows are vectorized covariate tensors;
- an array or `rTensor::Tensor` with dimensions `c(n, input_dims)`; or
- a length-`n` list of arrays or `rTensor::Tensor` objects.

Responses may be a length-`n` vector, an `n`-by-`prod(output_dims)` matrix, an
array or `rTensor::Tensor` with dimensions `c(n, output_dims)`, or a list of
response tensors. Public coefficient modes are always ordered as
**input modes followed by output modes**.

## CP decomposition

```r
library(SegreCP)

set.seed(1)
dims <- c(8, 9, 7)
u <- lapply(dims, function(p) {
  z <- rnorm(p)
  z / sqrt(sum(z^2))
})
signal <- array(3 * Reduce(kronecker, rev(u)), dims)
observed <- signal + array(rnorm(prod(dims), sd = 0.05), dims)

fit_rgd <- segre_decompose(
  observed,
  rank = 1,
  method = "rgd",
  init = "cpca",
  control = segre_control(profile = "manuscript", max_iter = 20)
)

fit_rgn <- segre_decompose(
  observed,
  rank = 1,
  method = "rgn",
  init = "hosvd",
  control = segre_control(profile = "manuscript", max_iter = 10)
)

reconstructed <- reconstruct_cp(fit_rgn)
summary(fit_rgn)
```

To retain an auditable optimization trajectory, set `keep_history = TRUE` and
use `convergence_data()`:

```r
fit_trace <- segre_decompose(
  observed, rank = 1, method = "rgn", init = "hosvd",
  control = segre_control(profile = "manuscript", max_iter = 10),
  keep_history = TRUE
)
trace <- convergence_data(fit_trace)
```

If simulated truth is supplied, `convergence_data(fit_trace, truth)` also
returns whole-tensor relative Frobenius error at initialization and after every
completed sweep.

For decomposition with the identity observation operator, RGN is exactly the
unit-step RGD component update.

## Tensor regression

```r
set.seed(2)
n <- 100
input_dims <- c(4, 5)
output_dims <- c(3, 2)
coefficient_dims <- c(input_dims, output_dims)

factors <- lapply(coefficient_dims, function(p) {
  z <- rnorm(p)
  z / sqrt(sum(z^2))
})
beta <- 2 * Reduce(kronecker, rev(factors))

X <- matrix(rnorm(n * prod(input_dims)), nrow = n)
beta_matrix <- matrix(beta, nrow = prod(input_dims),
                      ncol = prod(output_dims))
Y <- X %*% beta_matrix +
  matrix(rnorm(n * prod(output_dims), sd = 0.05), nrow = n)

fit <- segre_regress(
  X, Y,
  rank = 1,
  input_dims = input_dims,
  output_dims = output_dims,
  method = "rgn",
  init = "cpca",
  control = segre_control(profile = "manuscript", max_iter = 10)
)

dim(coef(fit))
dim(predict(fit, X[1:4, , drop = FALSE]))
```

## Initialization

Both required spectral initializers are available:

```r
cp_initialize(observed, rank = 1, method = "cpca")
cp_initialize(observed, rank = 1, method = "hosvd")
```

`"random"`, a supplied list of CP components, and a custom function are also
supported. A custom function is called as `function(x, rank, ...)` and returns
either a component list or `list(comps = component_list)`.

## Manuscript and practical profiles

`segre_control(profile = "manuscript")` uses the current displayed algorithms:
sequential residual refresh, fixed component order, fixed RGD step size,
exact structured T-HOSVD, a fixed sweep cap, and zero RGN ridge. Regression
inputs are normalized by `sqrt(n)` unless they are already on that scale.

`segre_control(profile = "practical")` makes empirical safeguards explicit:
decreasing-weight component order, an observable stopping tolerance,
backtracking for regression RGD, and a small recorded RGN ridge. These settings
are useful in applications but are not identical to the displayed manuscript
algorithm.

See `vignette("algorithm-fidelity", package = "SegreCP")` and
`inst/ALGORITHM_FIDELITY.md` in the source tree for the exact implementation
boundary.

## Migrating the existing experiment scripts

- Existing regression objects that already store both `X` and `Y` divided by
  `sqrt(n)` must use `normalize_regression = "none"`; raw observations should
  use the default `"sqrt_n"` normalization.
- The active legacy tensor-response implementation stores factors output-first.
  The package API is input-first, so legacy custom components and explicit CPCA
  splits must be permuted before they are supplied.
- Supply `init_split` when exact replay of a legacy CPCA tie choice is required.
  The package's automatic tie rule is deterministic.
- Historical Jacobi/simultaneous-update benchmark variants are outside version
  0.1.1. The package covers the current manuscript and main sequential paths.

## Rebuttal-stage numerical reproducibility

The package documentation reproduces numerical experiments reported as tables
during the rebuttal. The main seven-method regression vignette does not use the
former manuscript Figure 1 data. The initializer vignette displays only the
new rebuttal CPCA-versus-HOSVD ablation; it does not read or publish the old
figure's image or plotting output.

```r
vignette("rebuttal-regression-reproducibility", package = "SegreCP")
vignette("rebuttal-initialization-convergence", package = "SegreCP")
```

The main vignette rebuilds the seven-method tensor-on-tensor regression table
from 1,680 paired replicate-level records and visualizes the recorded
checkpoints, high-coherence stress cells, and accuracy--time relationship. The
second vignette shows CPCA/HOSVD sensitivity, all ten decomposition sweeps, and
observed-objective convergence for CP-ALS, RGD, and RGN. CP-ALS is read from a
frozen experiment artifact; it is not a callable package estimator.

To rebuild the rebuttal tables without rerunning the expensive benchmark:

```r
driver <- system.file(
  "reproduce", "validate_rebuttal_tables.R", package = "SegreCP"
)
system2(
  file.path(R.home("bin"), "Rscript"),
  c(driver, "--output", "rebuttal_table_rebuild")
)
```

To use only the authors' methods in a repository-level experiment:

```r
ours <- system.file("scripts", "ours_methods.R", package = "SegreCP")
source(ours)
# fit_ours_decomposition(..., method = "RGN", initialization = "cpca")
# fit_ours_regression(..., method = "RGD", initialization = "hosvd")
```

The frozen initialization artifact was generated before the code-scope split
and includes CP-ALS records. Its exact mixed runner remains recoverable from the
repository-only archived 0.1.0 source tarball and is not installed by version
0.1.1. Runtime validation uses the well-separated tier order
`{RGD, RGN} < CP-ALS < HOPLS < TT-IHT < {RRR, Tucker-RGD}`. Absolute seconds
and the order within the two near-tied pairs are not portability claims. The
regression initializer rows are reconstructed from frozen rebuttal records; a
full public-API regression initializer rerun is not claimed here.

## Tests

```r
testthat::test_local("software/SegreCP")
```

The tests include the installed-script scope boundary, direct
tangent-projection and dense normal-equation oracles,
RGD/RGN equivalence for decomposition, exact structured T-HOSVD comparisons,
CPCA/HOSVD initialization, data-format parity, input validation, frozen
rebuttal-table reconstruction, runtime-tier checks, and convergence-history
contracts.

See [`inst/VALIDATION.md`](inst/VALIDATION.md) for the recorded package-check and
active-source parity results and
[`inst/RELEASE_PROVENANCE.md`](inst/RELEASE_PROVENANCE.md) for the exact
historical-source and public-release boundaries.
