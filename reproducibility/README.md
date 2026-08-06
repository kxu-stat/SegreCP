# Repository-only comparison code

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

