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
