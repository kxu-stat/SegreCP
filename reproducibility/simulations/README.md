# Numerical experiments

These scripts generate new data under the rebuttal-stage simulation designs
and fit RGD and RGN with the installed `SegreCP` package. Install the package,
run the commands from the repository root, and choose new output directories
outside the repository.

## CP decomposition

```sh
Rscript reproducibility/simulations/run_decomposition.R \
  --output-dir=/path/to/new-smoke-results
Rscript reproducibility/simulations/run_decomposition.R \
  --full --output-dir=/path/to/new-full-results
```

The default runs one replication for three sweeps. `--full` runs 20
replications for 10 sweeps. The design uses a `60 x 30 x 30` tensor, rank 3,
condition number 30, coherence 0.5, Gaussian noise standard deviation 1,
signal multiplier 3, and seed `20260829 + replication`. Each generated tensor
has CPCA and HOSVD starts, shared between RGD and RGN. RGD uses step size 0.5.
The HOSVD start is ordered by signed core weights for the first sweep.

`trajectory.csv` records the objective and relative Frobenius error at each
sweep. `terminal_accuracy.csv` gives the final values; `settings.csv`,
`validation_checks.csv`, and `session_info.txt` describe the run. When bundled
reference trajectories are available, `reference_comparison.csv` summarizes
their numerical differences from the fresh fits. The objective is the
optimized quantity; reconstruction error against the known tensor may rise
near the noise floor.

## Tensor-on-tensor regression

```sh
Rscript reproducibility/simulations/run_tensor_regression.R \
  --output-dir=/absolute/path/to/new-smoke-results
Rscript reproducibility/simulations/run_tensor_regression.R \
  --full --output-dir=/absolute/path/to/new-full-results
```

The default runs one design cell, one replication, and two sweeps; `--help`
lists selectors for a different quick run. `--full` covers 12 cells, 20 paired
replications per cell, and 10 sweeps. The coefficient has dimensions
`10 x 10 x 10 x 10 x 10` and rank 3, with 1,000 training and 1,000 test
observations. Coherence is `0, 0.2, 0.5, 0.75` and noise standard deviation is
`0, 0.5, 1`. The same CPCA start is used for RGD and RGN in each task. The
script generates the coefficient in response-first order and converts it to
the package's predictor-first format for estimation.

`authors_regression_metrics.csv` records coefficient and prediction errors
and locally measured runtimes. `authors_regression_objective.csv` gives the
convergence paths. `run_settings.csv` and `sessionInfo.txt` record the design
and software environment. Results are written to the selected output directory.

To check both entry points with small runs after installing the package:

```sh
Rscript reproducibility/simulations/tests/smoke.R
```
