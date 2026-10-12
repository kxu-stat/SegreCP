# Validation record

Date: 2026-08-04

Version-0.1.2 package check update: 2026-10-03.
Version-0.1.3 source-replay update: 2026-10-03.
Version-0.1.4 reader-facing package and fresh-simulation update: 2026-10-03.

This record separates package mechanics, executable numerical oracles, frozen
rebuttal-artifact reconstruction, fresh package reruns, and migration
comparisons with the active manuscript code. It is not a theorem or proof
audit.

## Installed code scope

The exported namespace and the installed experiment script implement only RGD,
RGN, and their CPCA/HOSVD/random/user-supplied initialization layer. Competing
method names occur in frozen data, reconstruction tests, and vignettes, but no
competing optimizer is installed or exported. The sole current comparator
entry point is repository-only `reproducibility/comparison_methods.R`, which is
excluded by `.Rbuildignore`.

## Historical version-0.1.1 package check

The source package, including all four vignettes, was built and checked on macOS
arm64 with R 4.4.3:

```sh
cd software
R CMD build SegreCP
R CMD check --no-manual SegreCP_0.1.1.tar.gz
```

Final status: `OK` (0 errors, 0 warnings, 0 notes). All examples, vignettes,
compiled-code checks, Rd checks, and 156 testthat expectations passed; there
were no test warnings or skips.

The version-0.1.2 source tarball was built and checked on macOS arm64 with
R 4.4.3 using `R CMD check --no-manual`. Final status was `OK`: zero errors,
warnings, or notes; all five vignettes built and all 156 testthat expectations
passed. The built tarball excludes `reproducibility/` and contains no real
trade tensors, design rows, or downloaded source files.

The version-0.1.3 source tarball was built and checked on the same platform
with `R CMD check --no-manual`: final status `OK`, zero package errors,
warnings, or notes, all five vignettes rebuilt, and all 156 testthat
expectations passed without skips. A separate source-checkout test parsed the
12 real-data R scripts, validated all 82 source and four derived manifest rows,
checked both command-line help modes and package write guards. The detached
checkout comparator test passed its routing contract and a two-sweep synthetic
CP-ALS smoke fit. These two repository-only tests are intentionally outside
the installed package's testthat suite. The source tarball was inspected to
confirm that no real trade data or `reproducibility/` scripts were included.

The version-0.1.4 source tarball built and passed `R CMD check --no-manual`
on the same macOS arm64/R 4.4.3 platform: `Status: OK`, zero errors, warnings,
or notes. All six vignettes rebuilt and all 156 testthat expectations passed
without skips. In the source checkout, `reproducibility/simulations/tests/smoke.R`
also passed for both a fresh decomposition fit and a fresh tensor-regression
fit using the checked 0.1.4 package. These source-only scripts are excluded
from the tarball; they are not installed estimators.
The real-data contract check also parsed 13 R scripts and validated all 82
source and four derived manifest entries; the detached-checkout comparator
contract and synthetic CP-ALS smoke fit passed.

The new getting-started vignette demonstrates the public decomposition,
scalar-response regression, tensor-response regression, and convergence APIs
on self-generated data. The real-trade vignette now accepts a private data
directory through a source-only reader, rather than asking users to type an
internal historical dataset label.

## Version-0.1.2 real-trade authors-method replay

The repository-only `reproducibility/real_data/fit_authors_methods.R` was run
with the installed version-0.1.2 package on the four frozen v3 processed
inputs. Their recorded SHA-256 hashes matched 4/4. The runner used the frozen
five folds, 30 complete sweeps, CPCA, and five repetitions per fold, for 50
RGD/RGN fits. A row-wise comparison with the archived v3 seven-method raw
record matched all 50 fold-by-repetition-by-method keys. Maximum absolute
differences were `1.02e-14` for training RMSE, `9.77e-15` for test RMSE,
`2.00e-15` for training R-squared, and `1.30e-14` for test R-squared. Mean
test RMSE was 2.700299739447382 for RGD and 2.668389603947836 for RGN,
matching the archived table to printed precision. Hardware-dependent fit
seconds were not used as parity criteria. Against the frozen comparator rows,
RGN remains second and RGD seventh by mean test RMSE; the five comparators
were not rerun through this package.

This was a fresh package-API replay of the **authors' two methods conditional
on the frozen processed inputs**, not a fresh seven-method refit. At that
0.1.2 checkpoint the local source tree lacked the CEPII Gravity ZIP and RDS;
both were recovered and the upstream replay was validated separately for
0.1.3 below. The historical sample screen uses 2019 source-availability
information, so it is not strictly selected using only 2018 information.

## Version-0.1.3 official-source-to-input replay

The missing `Gravity_rds_V202211.zip` was downloaded from the official CEPII
route, and both its archive and extracted `Gravity_V202211.rds` matched the
frozen byte sizes and SHA-256 hashes. The RDS is readable (4,699,296 rows,
87 columns) and has the required 2018 and 2019 observations. It remains in
the private real-data directory, not in SegreCP source or its tarball.

In an isolated private directory, the repository-only scripts reconstructed
five reporter-query metadata files from official raw inputs; each was bytewise
identical to its historical counterpart. A second isolated directory initially
contained only `raw_official/` and `source_archives/`, with no prebuilt
metadata or processed outputs. A single `run.R --prepare --offline` command
checked all 82 frozen official-source content hashes, generated its own query
metadata, parsed both years of WITS/TRAINS records, and regenerated the
covariate panel, tensors, selected six-feature design, common folds, and
registry. All four selected analysis inputs matched the historical files byte
for byte and by SHA-256. The reconstructed tensors
had `X = 80 x 25 x 21 x 6` and `Y = 80 x 25 x 21 x 8`; the design had 80 rows
and six feature columns, and the five folds had 16 exporters each.

This validates a local **frozen-official-source-to-input replay**. Together
with the separate 50-fit RGD/RGN package-API parity on the exact same input
bytes, it supports the authors' two-method result chain. It does not show
that a fresh download from every upstream API still succeeds, nor rerun the
other five estimators. A fresh full-download test would have to revalidate
current WITS API availability and response content. The repository-level,
data-free contract test additionally checks all 82 source and four derived
manifest rows, script syntax, CLI help, and safeguards against writing private
data into the package source.

## Executable package oracles

Tests under `tests/testthat/` include:

- explicit tangent-basis projection checks;
- a direct full rank-one T-HOSVD oracle for a one-step decomposition update;
- rank-two sequential residual-refresh and component-order oracles;
- scalar and tensor-response dense tangent normal-equation RGN oracles;
- an exact-quadratic RGD line-search oracle;
- normalized-versus-pre-scaled regression parity;
- monotone backtracking with an exact no-op fallback;
- CPCA's displayed product-scale split rule and sign/weight conventions;
- CPCA and HOSVD recovery on an orthogonal CP tensor; and
- base array, matrix, list, and `rTensor::Tensor` format parity;
- component-history and `convergence_data()` alignment; and
- frozen rebuttal means, runtime tiers, provenance gates, and convergence
  contracts.

## Rebuttal-stage regression table reconstruction

The main reproducibility vignette and
`inst/reproduce/validate_rebuttal_tables.R` read replicate-level artifacts,
not the formatted response table. The deterministic reconstruction verifies:

- all 1,680 rows in the `4 x 3 x 20 x 7` paired grid;
- one shared seed and all seven methods in every setting/replication pair;
- all 360 six-method checkpoint records at iterations 3, 5, and 10;
- the SHA-256 hash of the derived checkpoint archive;
- the rebuttal Panel A, B, and replication-balanced Panel C means; and
- the stable runtime tier
  `{RGD, RGN} < CP-ALS < HOPLS < TT-IHT < {RRR, Tucker-RGD}`.

Exact seconds and the within-tier order of RGD/RGN or RRR/Tucker-RGD are not
portable pass/fail conditions. The main seven-method vignette does not use the
old manuscript Figure 1 data. The initializer vignette uses only the numerical
CPCA-versus-HOSVD ablation reported in the rebuttal and does not read the old
figure's image or plotting output.

The final grid's source manifest records the historical solver with SHA-256
`50ee798a...`. The later active repository solver has changed. Accordingly,
the package labels this stage as deterministic reconstruction of frozen raw
records, not as a bitwise fresh rerun from a silently substituted active file.

## Fresh authors-only simulations with the public API

The repository-only `reproducibility/simulations/run_decomposition.R` was run
for the complete 20-replication, 10-sweep rebuttal decomposition design. Its
generated-data fingerprints matched the frozen package rerun in all 20
replications. Its mean objective and error paths matched that package rerun to
floating-point precision (maximum mean-path discrepancy about `5.6e-17`;
maximum objective discrepancy about `5.1e-11`). Comparison with the earlier
historical reference retained the documented HOSVD stability limits: mean
error-path discrepancy below `1e-3`, iteration-3/5/10 discrepancy below
`5e-4`, and preserved RGD-versus-RGN checkpoint ranking. The HOSVD branch is
not claimed to be iterate-by-iterate identical to the historical solver.

The separate `run_tensor_regression.R` generates new data under the frozen
12-cell tensor-on-tensor design and fits only CPCA-initialized RGD and RGN.
Default and selected 10-sweep smoke runs succeeded. In checked cells 1 and 7,
the complete generated-data fingerprints matched the historical records. For
noisy cell 7, the public CPCA initial tensor differed from the historical one
by only `1.31e-14` in relative Frobenius norm. The current active manuscript
solver and the package then produced the same 10-sweep method metrics; the
frozen seven-method table did not. Its recorded solver hash differs from the
current source hash. This is evidence of a solver-version difference after
the historical run, not a data-generation or API mode-order mismatch. The
exact historical code change is not available here, so the new script is a
fresh current-method study, **not** pointwise reproduction of the frozen
regression table or a fresh seven-method benchmark. Compare historical rows
by coherence, noise, replication, and method; its `cell_id` is not globally
unique across coherence subruns.

The complete 0.1.4 tensor-regression rerun subsequently finished all 12 cells
and 20 replications per cell: 480 finite rank-three RGD/RGN fits, each with 10
sweeps, and 5,280 objective rows. All 240 data and base-data fingerprints
matched the historical generated-data archive. Relative to the frozen
RGD/RGN rows, the maximum absolute coefficient-error difference was
`0.01573823`, and the maximum absolute test-RMSE difference was `0.0002788839`;
neither is reported as historical parity. The lower mean test-RMSE author
method agreed in 11 of 12 cells; the sole reversal was a near tie at
coherence `0.75`, noise standard deviation `0.5`. Three fixed-step RGN paths
had small observed-objective increases (largest relative to initial loss
`2.93e-6`), so the fresh regression study does not establish monotonicity of
every path.

## Initialization sensitivity and empirical convergence

The bundled legacy decomposition study retains its original 8-of-8 report. The
regression reference archives the original 18-of-18 report; one historical
check concerns four candidate plot files that are intentionally not bundled or
used by this package. Package-subset manifests separately hash every numerical
file actually shipped. The archived version-0.1.0 formal package run used the
rebuttal decomposition design, 20 paired replications, 10 sweeps, and the
recorded seed schedule. It passed 8 of 8 checks, including non-increasing
observed loss in all `20 x 2 x 3 = 120` CP-ALS/RGD/RGN paths. Version 0.1.1
does not install that mixed runner; the exact pre-split source tarball is
retained outside the built package, and the frozen numerical bundle and its
manifests are unchanged.

CPCA paths match the legacy artifact within `2.74e-14`, and both CPCA/HOSVD
CP-ALS paths match within `3.11e-14`. The historical HOSVD RGD/RGN code used a
signed-weight component order, while the public package uses canonical
nonnegative weights. That branch is therefore validated as a stability rerun:
the maximum mean-path discrepancy is `8.68e-4`, the maximum iteration-3/5/10
discrepancy is `3.08e-4`, and all recorded checkpoint method ranks and
CPCA-versus-HOSVD contrasts are preserved.

Observed loss is the optimization target. Oracle error to the simulated truth
can increase slightly in individual paths near the noise floor even when loss
continues to decrease. The vignettes therefore report empirical convergence
diagnostics, not a proof of global convergence or per-replication monotonicity
of oracle error.

## Active-source migration parity

A separate repository-level migration harness matched initialization, mode
order, normalization, component order, and step rules between this package and
the active source files:

- `rgn_benchmark_decomp/tensor-CP-decomposition.R`;
- `rgn_benchmark/tensor-regression.R`.

The maximum entrywise discrepancies after matched one-step updates were:

| Route | RGD | RGN |
|---|---:|---:|
| CP decomposition | 4.44e-16 | 2.78e-16 |
| Scalar-on-tensor regression | 4.00e-14 | 2.78e-16 |
| Tensor-on-tensor regression | 9.24e-14 | 1.11e-16 |

Matched regression step-rule discrepancies were 2.22e-16 for
`exact_quadratic` and 1.05e-15 for backtracking. Four-sweep rank-two checks were
also at floating-point precision: 3.33e-16 for decomposition RGD, at most
2.22e-16 for scalar regression RGD, and 8.88e-16 for scalar regression RGN.

These comparisons cover the current sequential/Gauss--Seidel routes. They do
not claim parity with historical Jacobi benchmark variants or archived frozen
supplementary snapshots.
