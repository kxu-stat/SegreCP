# SegreCP 0.1.4

- Reorganized the public README around installation and executable examples
  for CP decomposition, scalar-on-tensor regression, and tensor-on-tensor
  regression; added a getting-started vignette with an objective trace.
- Kept the historical real-trade design identifier inside the source-only
  reproducibility layer. A reader-facing helper now loads and verifies the
  prepared design, response, and folds without exposing archival filenames.
- Gave the simulation vignettes task-focused titles and separated frozen-record
  table reconstruction from fresh authors-method simulation reruns.

# SegreCP 0.1.3

- Restored the pinned CEPII Gravity V202211 raw archive and RDS in the private
  research data directory after matching their frozen sizes and SHA-256 hashes;
  no third-party numerical data are included in the package or source checkout.
- Made the real-trade vignette and source-checkout instructions standalone,
  clarified the trade-design archival identifier, and repaired the convergence
  plotting example to use the recorded repetition column.
- Verified a one-command offline rebuild from 82 hash-matched official source
  files to four byte-identical analysis inputs; offline mode now reconstructs
  reporter-query metadata from raw sources automatically.
- Added fast, data-free contract checks for the repository-only data workflow
  and help text for the authors-method fitting command. The separate CP-ALS
  comparison script now works from a standalone source checkout.

# SegreCP 0.1.2

- Added a repository-only, data-free real-trade workflow with official-source
  download, preprocessing, frozen-input checks, authors-method fitting, and a
  Markdown vignette.
- Documented the historical 2019 source-availability sample screen and the
  boundary between the 50-row RGD/RGN numerical parity check and a full
  seven-method, source-to-result reproduction.

# SegreCP 0.1.1

- Restricted the installed experiment-facing code to RGD, RGN, and the
  initialization layer; no competing optimizer is installed or exported.
- Added `inst/scripts/ours_methods.R` as the single documented entry point for
  authors' methods and moved comparison code to the repository-only,
  package-excluded `reproducibility/comparison_methods.R` script.
- Preserved the exact pre-split 0.1.0 source tarball and its mixed-runner hash
  without rewriting the frozen rebuttal artifacts or their manifests.

# SegreCP 0.1.0

- Added sequential RGD and RGN for CP decomposition, scalar-on-tensor
  regression, and tensor-on-tensor regression.
- Added exact structured rank-`(1, ..., 1)` T-HOSVD retraction kernels.
- Added CPCA, HOSVD, random, supplied-component, and custom-function
  initialization routes.
- Added base R array/matrix, list-of-tensors, and `rTensor::Tensor`
  compatibility.
- Added manuscript and practical control profiles, documentation, vignettes,
  and direct numerical-oracle tests.
- Added optional component-history recording and `convergence_data()` for
  observed-objective and simulated-truth trajectories.
- Added frozen replicate-level rebuttal artifacts, a deterministic seven-method
  table validator, and a formal CPCA/HOSVD initialization rerun with SHA-256
  provenance and legacy-parity gates.
- Added vignettes that visualize the rebuttal-stage regression tables,
  initialization sensitivity, complete ten-sweep trajectories, and
  objective convergence. Runtime validation uses separated tiers rather than
  hardware-specific seconds or within-near-tie ordering.
