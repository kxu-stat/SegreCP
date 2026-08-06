# Algorithm fidelity and scope

This file records the implementation contract for `SegreCP` 0.1.1. It is
intended to prevent a practical software option from being mistaken for a
displayed manuscript algorithm.

## Live references

- Current manuscript: `neurips_2026/neurips2026.tex`.
- Active decomposition implementation used by the current experiment pipeline:
  `rgn_benchmark_decomp/tensor-CP-decomposition.R`.
- Active regression implementation used by the current experiment pipeline:
  `rgn_benchmark/tensor-regression.R`.

Archived files under `supplementary_material_submission/` are historical
snapshots and are not package sources.

## Shared algorithmic invariants

Both `segre_decompose()` and `segre_regress()` implement complete sequential
(Gauss--Seidel) sweeps. After component `i` is updated, its cache, prediction,
and residual are refreshed immediately before component `i + 1` is processed.
The default manuscript component order is `1, ..., r`.

Every component is represented by a positive weight and unit-norm mode factors.
Each RGD or RGN tangent proposal is retracted with the exact
rank-`(1, ..., 1)` T-HOSVD. Because the proposal is a rank-one tensor plus a
tangent tensor, every leading mode singular vector is recovered from a `2 x 2`
symmetric Gram eigenproblem; no ambient unfolding and no one-step power-method
approximation is used.

## RGD

For component `i`, the package forms the current partial least-squares residual,
projects the partial Euclidean gradient onto the current Segre tangent space,
takes the requested tangent step, and applies the exact retraction. Under the
manuscript profile, component order and step sizes are fixed. The default RGD
step size is `0.5`, matching the manuscript's displayed fixed-step setting;
experiment-specific scripts may override the step rule or step size.

The practical regression profile may instead use monotone backtracking. This is
an explicit empirical safeguard, not a claim about the fixed-step theorem.

## RGN

For CP decomposition, the observation operator is the identity. Its restricted
tangent normal operator is therefore the identity, so the RGN component update
is the unit-step RGD component update.

For scalar-on-tensor and tensor-on-tensor regression, the package solves the
restricted tangent least-squares normal equation. Tensor-response sufficient
statistics are contracted modewise without materializing an ambient tangent
basis. The manuscript profile sets `ridge = 0`. A positive `ridge` is available
only as an explicitly recorded regularized practical variant.

## Initialization

`cp_initialize(method = "cpca")` implements Composite PCA with a nonempty proper
mode split. The automatic split is deterministic; refolded singular-vector
signs are aligned; weights are re-estimated by projection; and a negative sign
is absorbed into the first factor.

`cp_initialize(method = "hosvd")` takes leading mode singular vectors and the
diagonal entries of the multilinear HOSVD core. This is a warm-start method,
not a general CP-decomposition solver. Random, supplied-component, and custom
function initializers use the same validated component interface.

## Data and normalization

Base R arrays and dense matrices are the internal representation. An
`rTensor::Tensor` is converted once at the package boundary. For regression,
coefficient modes are public in input-first/output-second order; a fixed
output-first permutation is used only inside the structured contractions and is
reversed before returning the fit.

With `normalize_regression = "sqrt_n"`, both design and response are divided by
`sqrt(n)`, matching the manuscript's normalized observation operator. Use
`"none"` only when data have already been normalized or an unnormalized
least-squares convention is intended.

## Deliberate scope boundary

Version 0.1.1 exposes the two observation structures used by the current
software: identity-operator decomposition and tensor regression. It does not
yet expose an arbitrary user-defined linear operator. It also does not expose
historical Jacobi/simultaneous-update benchmark variants or any competing
estimator. Consequently, package documentation should claim fidelity to the
current manuscript algorithms and the active sequential implementation, not to
every archived experimental branch.

## Numerical acceptance criteria

The package tests compare:

- tangent projections with explicit tangent bases;
- one-step RGD with a direct projection plus full rank-one T-HOSVD;
- decomposition RGN with unit-step RGD;
- scalar and tensor-response RGN with dense tangent normal-equation oracles;
- CPCA and HOSVD on an exactly recoverable orthogonal CP tensor; and
- base-array, matrix, list, and `rTensor` input routes.

Build success and numerical tests are reported separately from any theorem or
proof audit.
