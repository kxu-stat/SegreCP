# Repository-only reproducibility code

This directory is intentionally excluded from the built `SegreCP` source
package. The installed package exposes only the authors' RGD and RGN methods,
their initialization layer, and supporting data/diagnostic utilities.

`comparison_methods.R` is the single current experiment-facing entry point for
competing estimators. It contains the CP-ALS path formerly mixed into the
initialization study and records the task-specific workers used for the other
comparators. Those routines are experiment adapters, not supported `SegreCP`
estimators.

## Historical mixed runner

The frozen `decomp_pkg` artifact was created before this scope split. Its source
manifest records the mixed runner at
`inst/reproduce/rebuttal_initialization_sensitivity.R` with SHA-256
`fa2d3ef0bdc6e4827f50a9caad7a7f00cf9ade0d3a6e5844b311d91276079de2`.
That exact file remains recoverable as a member of
`archive/SegreCP_0.1.0_pre_scope_split.tar.gz`; `LEGACY_SOURCE_MAP.csv` records
the member path and hash. The frozen bundle and each of its manifests are left
unchanged.

The current authors-only entry point is installed at
`inst/scripts/ours_methods.R`. New runs should use that script for CPCA, HOSVD,
RGD, and RGN, and use `comparison_methods.R` separately when a comparator is
required.

## Fresh synthetic experiments

`simulations/README.md` gives two source-checkout commands for each task:
a quick local run and an explicit full-design run. `run_decomposition.R`
regenerates the rebuttal decomposition tensors, then fits CPCA/HOSVD-initialized
RGD and RGN. `run_tensor_regression.R` regenerates the paired
tensor-on-tensor regression grid and fits RGD/RGN with shared CPCA. Use new
output paths outside the package; neither script changes the installed methods.

The full decomposition rerun agrees with the frozen package study to numerical
precision, with the historical HOSVD comparison assessed by prespecified
stability tolerances. The full tensor-regression rerun also finished all
12 cells and 20 replications per cell: all 240 generated-data fingerprints
matched the historical archive, and 480 RGD/RGN fits produced finite metrics.
The noisy results do **not** reproduce the old method metrics pointwise because
the archived benchmark used a different solver version. This is an independent
current-method rerun, not a replay of the seven-method table. The bundled
frozen tables can still be inspected and redrawn, but redrawing them must not
be described as fresh method fitting.

## Real trade data preparation

See `real_data/README.md` for a source-pinned workflow that downloads official
trade inputs into a private directory and rebuilds the frozen 2018–2019 trade
tensor, design matrix, and common folds. This folder contains code, URLs, and
hashes only; it does not distribute raw or processed third-party data. Its strict
`run.R --check-inputs` mode checks both frozen sources and derived artifacts.
The data-preparation workflow alone does not recreate the historical
seven-method numerical table. `real_data/fit_authors_methods.R` separately
checks the four frozen processed inputs and refits RGD/RGN with the public
package; see its README for smoke and full five-fold commands.
