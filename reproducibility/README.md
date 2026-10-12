# Reproduce the experiments

Install `SegreCP` first, then run these scripts from the repository root.
They use the package's RGD and RGN methods with CPCA or HOSVD initialization.

## Synthetic experiments

[`simulations/README.md`](simulations/README.md) gives commands, seeds, and
settings for fresh noisy CP decomposition and tensor-on-tensor regression runs.
The [initialization](../vignettes/rebuttal-initialization-convergence.Rmd) and
[tensor-regression](../vignettes/rebuttal-regression-reproducibility.Rmd)
vignettes recompute the rebuttal-stage tables and plots from per-replication
results.

The comparison script is [`comparison_methods.R`](comparison_methods.R). Its
CP-ALS decomposition route runs here; the regression comparison workers use
the larger manuscript experiment workspace.

## Real-data analysis

[`real_data/README.md`](real_data/README.md) gives official data sources,
preprocessing commands, and the five-fold RGD/RGN analysis workflow. Downloaded
and processed trade data stay outside this repository.
