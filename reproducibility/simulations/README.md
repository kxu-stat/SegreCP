# Fresh simulations with RGD and RGN

These repository-only scripts generate new simulated data and fit **RGD and
RGN only** through the public `SegreCP` API. They do not read stored outcomes
as substitutes for new fits, implement comparison methods, or include generated
data and results in the source tree.

Install `SegreCP` from its source directory and install `digest` for
generated-data fingerprints. Start the commands below in the `SegreCP` source
directory and choose new output paths outside that directory.

## CP decomposition

```sh
Rscript reproducibility/simulations/run_decomposition.R \
  --output-dir=/path/to/new-smoke-results
Rscript reproducibility/simulations/run_decomposition.R \
  --full --output-dir=/path/to/new-full-results
```

The first command runs one replication for three sweeps; the explicit `--full`
option runs the frozen 20-replication, 10-sweep design. The script refuses to
overwrite an existing output path. The full design uses a `60 x 30 x 30`
tensor, rank 3, condition number 30, coherence 0.5, Gaussian noise standard
deviation 1, signal multiplier 3, and seed `20260829 + replication`. Each
generated tensor receives CPCA and HOSVD initializers; RGD and RGN share each
realized initializer. RGD uses step size 0.5. All fits use sequential
Gauss--Seidel updates, exact rank-one retraction, and no early stopping.

`trajectory.csv` contains observed least-squares objective and whole-tensor
relative Frobenius error against the generated truth at initialization and
every completed sweep. `terminal_accuracy.csv` contains terminal values.
`seed_schedule_and_data_fingerprints.csv`, `settings.csv`,
`validation_checks.csv`, and `session_info.txt` record provenance and checks.
The objective should not increase; oracle error can increase near the noise
floor because it is not the optimized quantity.

If the installed package contains the frozen rebuttal tables,
`reference_comparison.csv` reports numeric differences *after* fitting. For a
full run, the frozen package rerun should match to numerical precision under
the same R and platform environment. Against the earlier historical solver,
compare CPCA
pointwise (tolerance `1e-10`), but assess HOSVD using the 20-replication mean
path (`1e-3`) and iterations 3, 5, and 10 (`5e-4`). The one-replication smoke
comparison is diagnostic only and cannot validate a mean-path claim. Different
R/BLAS versions may change fingerprints or last-bit numeric agreement. The
historical reference used a different serialized data representation, so its
data fingerprints are not treated as directly comparable with this script's;
the frozen package-rerun fingerprints are directly comparable. The comparison
also records whether the RGD-versus-RGN error ranking at the checkpoints is
preserved; it says nothing about the excluded comparison methods.

The historical HOSVD code used signed core entries to order its first sweep.
This script matches that first-sweep order, then uses SegreCP's canonical
nonnegative-weight order. Thus CPCA is eligible for near-exact comparison with
the archived reference, while HOSVD should be assessed by mean-path and
checkpoint stability, not iterate-by-iterate equality. This is a fresh
authors-method rerun, **not** a fresh rerun of CP-ALS or the seven-method
regression table. The frozen reference tables and historical comparison code
are documented separately under `inst/extdata/rebuttal_init/`.

## Tensor-on-tensor regression

```sh
Rscript reproducibility/simulations/run_tensor_regression.R \
  --output-dir=/absolute/private/new-smoke-results
Rscript reproducibility/simulations/run_tensor_regression.R \
  --full --output-dir=/absolute/private/new-full-results
```

The default runs one cell, one replication, and two sweeps; `--help` describes
optional smoke selectors. `--full` requests the rebuttal grid: coefficient
dimensions `10 x 10 x 10 x 10 x 10`, rank 3, condition number 30, 1,000
training and 1,000 test observations, coherence `0, 0.2, 0.5, 0.75`, noise
standard deviation `0, 0.5, 1`, 20 paired replications per cell, and 10
sweeps. The two methods share a CPCA initializer in each generated task. The
public coefficient has predictor modes first and response modes last; the
script transposes the historical output-first generator accordingly.

`authors_regression_metrics.csv` records coefficient error and train/test
prediction errors for both methods. `authors_regression_objective.csv` records
the optimization paths. `generated_data_fingerprints.csv`,
`run_settings.csv`, and `sessionInfo.txt` record the data seeds, design, and
software environment. Outputs go only to the user-selected new directory.

This is an independent current-package rerun, **not** a reconstruction of the
frozen seven-method table from its CSVs and **not** a fresh seven-method
benchmark. Generated-data fingerprints matched the historical inputs for
replication 1 in checked cells 1 and 7. However, for noisy cell 7
(`coherence = 0.5`, `noise_sd = 0.5`), the current RGD/RGN metrics did not
agree pointwise with the historical table. They did agree with a rerun of the
current manuscript solver on the same generated data; the archived table
records a different solver-file hash. The specific historical code change is
not available in this checkout. The complete `--full` run also finished with
480 finite rank-three method fits and 5,280 objective records. All 240 paired
generated-data fingerprints matched the historical archive. In 11 of 12
cells, the lower mean test-RMSE method among RGD/RGN stayed the same; the one
changed cell (`coherence = 0.75`, `noise_sd = 0.5`) is a near tie. Three RGN
paths had a small observed-objective increase under the fixed-step paper
profile, so the fresh run is not evidence of monotonic loss in every path.
The script supports a fresh test of the authors' algorithms under the rebuttal
design, but does not claim exact historical numerical parity or portable
wall-clock timing. Match archived rows by coherence, noise, replication, and
method: the archived `cell_id` restarts for each coherence subrun.

For an automated smoke check of both simulation entrances after installation,
run:

```sh
Rscript reproducibility/simulations/tests/smoke.R
```
