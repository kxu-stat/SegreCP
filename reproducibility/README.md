# Reproduce the experiments

Install `SegreCP >= 0.1.5` first, then run these scripts from the repository
root. Fresh authors-method fits use the accelerated public RGD/RGN API
with zero ridge regularization.

## Synthetic experiments

[`simulations/README.md`](simulations/README.md) gives commands, seeds, and
settings for the current `signal4_formal20` decomposition and tensor-on-tensor
regression study: 20 paired replications, both noise levels, and shared CPCA
starts. The default commands run small subsets before a full run.
The [initialization](../vignettes/rebuttal-initialization-convergence.Rmd) and
[tensor-regression](../vignettes/rebuttal-regression-reproducibility.Rmd)
vignettes recompute historical rebuttal-stage tables and plots from their
bundled per-replication records; these archives are separate from fresh runs.

The comparison script is [`comparison_methods.R`](comparison_methods.R). Its
CP-ALS decomposition route runs here. The separate
[`simulations/run_comparisons.R`](simulations/run_comparisons.R) entry point
uses the current zero-ridge regression adapters from the manuscript workspace.

## Real-data analysis

[`real_data/README.md`](real_data/README.md) gives official data sources,
preprocessing commands, and the five-fold RGD/RGN analysis workflow. Downloaded
and processed trade data stay outside this repository.
