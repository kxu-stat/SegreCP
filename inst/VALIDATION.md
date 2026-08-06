# Validation record

Date: 2026-08-04

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

## Package check

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
