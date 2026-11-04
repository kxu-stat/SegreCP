# Current numerical experiments

These entry points generate fresh data for the fixed `signal4_formal20`
study and call the accelerated public API in installed `SegreCP >= 0.1.5`.
They run RGD and RGN without ridge regularization. No manuscript workspace
or stored experimental outcomes are needed for these authors-method runs.

```sh
Rscript reproducibility/simulations/run_decomposition.R \
  --output-dir=/existing/parent/new-decomposition-smoke
Rscript reproducibility/simulations/run_tensor_regression.R \
  --output-dir=/existing/parent/new-regression-smoke
```

The default uses the first dimension setting, replication 1, both noise levels,
and three decomposition or two regression iterations. `--p-index`,
`--replication`, and `--iterations` select another smoke subset. Add `--full`
without these selectors for all dimensions, 20 independent paired replications,
and 10 iterations. Choose a new output directory outside the package checkout;
its parent directory must exist.

| Task | Dimensions | Full RGD/RGN fits |
|---|---|---:|
| Decomposition | `P x 30 x 30 x 30`, `P = 30, 60, 90, 120` | 320 |
| Tensor-on-tensor regression | input `P x 10 x 10`, output `10 x 10`, `P = 40, 60, 100` | 240 |

Both tasks use rank 3, condition number 30, coherence 0.5, signal multiplier 4,
and noise levels 0 and 1. Regression uses 1,000 training and 1,000 test samples.
The Mersenne-Twister seed is
`741103000 + task_offset + 1000 * p_index + replication`, with task offsets
0 for decomposition and 100000 for regression. The two noise levels share
factor directions and the same realized noise. The coefficient and noise
scaling follow the fixed study; regression uses a normalized design and
`normalize_regression = "none"`.

Within every dataset, RGD and RGN share the same realized CPCA start.
CPCA groups modes `1,2` for decomposition and output modes `4,5` for regression.
Both use decreasing-weight component order. Decomposition RGD uses step 0.5,
and RGN uses step 1. Regression RGD uses its exact-quadratic step rule; RGN
uses the unregularized tangent least-squares update.

`trajectories.csv` records relative coefficient recovery error and the observed
least-squares objective at every iteration, including initialization.
Regression also records train, test, and clean-test RMSE.
`fits.csv` records actual iterations, ridge, and the public solver entry point.
`settings.csv`, `settings.rds`, and `sessionInfo.txt` record the configuration.
`inputs.csv` contains input and initialization hashes; `cases/` saves compact
truth and initialization objects, and `fits/` saves every component state,
accepted step, control, and solver objective. Recovery error is reported as
computed, without smoothing or imposing monotonicity. These runs do not
measure solver speed.

The separate `run_comparisons.R` entry point runs comparison methods through
the current manuscript adapters; it requires the complete manuscript workspace.
The authors-only entry points above remain standalone. Bundled `inst/extdata`
rebuttal results are historical records, separate from this fresh study.

To check both authors-method entry points:

```sh
Rscript reproducibility/simulations/tests/smoke.R
```
