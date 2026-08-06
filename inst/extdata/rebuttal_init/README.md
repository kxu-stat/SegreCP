# Rebuttal-stage initializer sensitivity

This bundle supports the initializer-sensitivity and convergence vignette. It
shows only the CPCA-versus-HOSVD ablation reported during the rebuttal. The
historical regression run was produced in the same simulation workflow as an
earlier figure, but this package does not read or publish that figure's image or
plotting output.

## Frozen studies

- `decomp_ref/` contains the 20-replication, 10-sweep decomposition
  study with dimensions `60 x 30 x 30`, rank 3, condition number 30,
  coherence 0.5, noise standard deviation 1, and signal multiplier 3. Its
  frozen validation report passes 8 of 8 checks.
- `reg_ref/` contains the four rebuttal regression cells comparing
  CPCA and HOSVD. Its frozen validation report passes 18 of 18 checks.
- `decomp_pkg/` is the formal version-0.1.0 rerun through the public `SegreCP`
  API with the same decomposition design and seed schedule. It stores all
  10-sweep oracle-error and observed-objective paths and passes 8 of 8 package
  validation checks.

Each reference directory includes a `package_subset_manifest.csv` that hashes
the numerical files actually shipped. `historical_artifact_manifest.csv`
retains the broader source-run inventory, which can mention historical plot or
binary artifacts that are intentionally not bundled. `source_manifest.csv`
records the frozen runner and solver hashes.

The archived `decomp_pkg/` directory also has an artifact manifest and a package
subset manifest. Its `source_manifest.csv` records both portable relative
labels and the exact installed R database/shared-library hashes used for that
run. Those installed hashes are a macOS build fingerprint, not a
cross-platform source identity; portable R/C++ source hashes are recorded in
the package-level `inst/PACKAGE_SOURCE_MANIFEST.csv`, and the release tarball
has a separate adjacent SHA-256 checksum.

The three methods in a decomposition replication share the realized CPCA or
HOSVD initializer. Observed least-squares loss is non-increasing in all 120
stored package paths. Oracle error is not the optimized objective and can rise
slightly in individual runs near the noise floor.

## Historical HOSVD ordering

The historical rebuttal code ordered components using signed HOSVD diagonal
core entries. `SegreCP` canonicalizes every component to a nonnegative weight,
so that signed ordering is not a gauge-invariant public convention. CPCA
trajectories and CP-ALS trajectories match the legacy artifact to numerical
precision. For HOSVD RGD/RGN, the package artifact is a canonical stability
rerun: it passes the recorded mean-path and checkpoint tolerances and preserves
the checkpoint method ranks and CPCA-versus-HOSVD contrasts, but it is not
claimed to be an exact iterate-by-iterate legacy replay.

## Source boundary after the scope split

The version-0.1.0 runner combined CP-ALS with RGD/RGN and is therefore no
longer installed by version 0.1.1. In the repository checkout, its exact source
remains recoverable from
`reproducibility/archive/SegreCP_0.1.0_pre_scope_split.tar.gz`; the tar member
and original SHA-256 are recorded in `reproducibility/LEGACY_SOURCE_MAP.csv`.
The `decomp_pkg/` files and their manifests have not been rewritten.

For new work, source the installed `inst/scripts/ours_methods.R` entry point for
CPCA, HOSVD, RGD, and RGN. The single repository-only
`reproducibility/comparison_methods.R` entry point contains comparison code and
is excluded from the built package.
