# Tensor-on-Tensor Regression: Final Seven-Method Experiment

This directory contains the submission-facing source tables for the paired
synthetic tensor-on-tensor regression experiment.  Every replication uses the
same generated coefficient tensor, training sample, test sample, and noise
realization for all seven methods:

`CP-ALS`, `RRR`, `Tucker-RGD`, `HOPLS`, `TT-IHT`, `RGD`, and `RGN`.

The experiment uses coefficient dimensions `10 x 10 x 10 x 10 x 10`, two
response modes and three predictor modes, CP rank 3, condition number 30,
1,000 training and 1,000 test observations, 20 paired replications, coherence
`0, 0.2, 0.5, 0.75`, noise standard deviation `0, 0.5, 1`, and 10 updates.
CP-ALS, RGD, and RGN share the computed CPCA initialization.  The other
methods use their native parameterizations and initializations; HOPLS extracts
three latent blocks.

This is the experiment reported as tables during the rebuttal. It is not the
experiment behind the manuscript's former Figure 1.

## Frozen source tables

- `final_seven_method_grid_raw.csv`: all 1,680 method-cell-replication rows.
- `final_seven_method_replication_balanced.csv`: each metric first averaged
  equally over the 12 settings within a replication.
- `final_seven_method_grid_summary.csv`: the response-table means and standard
  errors over the 20 paired replication averages.
- `final_seven_method_noisy_grid_summary.csv`: the corresponding eight noisy
  settings only.
- `final_seven_method_cell_summary.csv`: method summaries for every coherence
  and noise cell.
- `final_seven_method_serial_timing_raw.csv` and
  `final_seven_method_serial_timing_summary.csv`: 20 serial runs at
  `(coherence, noise_sd) = (0.5, 0.5)`, including method-specific
  initialization.
- `final_settings.csv`, `final_validation_checks.csv`, and
  `final_source_artifact_sha256.csv`: frozen design, pairing checks, and hashes
  of the exact source code and result files.
- `checkpoint_metrics_raw.csv`: 360 replicate-level records at iterations 3,
  5, and 10 for the six methods with comparable outer-update checkpoints.
  HOPLS is absent because its latent blocks are not commensurate with those
  checkpoints.
- `checkpoint_source_manifest.csv`: hashes of the four frozen inputs used to
  assemble the checkpoint archive and of the derived archive itself.
- `package_subset_manifest.csv`: hashes of every data and documentation file
  actually shipped in this installed bundle.

The files are regenerated and validated from the completed run directories by

```bash
Rscript rgn_benchmark/freeze_final_tensor_on_tensor_seven_method.R
```

The underlying runner is
`rgn_benchmark/run_tensor_on_tensor_baseline_grid.R`; it records the complete
settings, seed schedule, generated-data fingerprints, validation checks, and R
session information in every run directory.  The exact completed source-run
paths are listed in `final_source_artifact_sha256.csv`.

The package vignette and `inst/reproduce/validate_rebuttal_tables.R` rebuild
Panels A--C deterministically from these replicate-level files. They do not
silently substitute a later active solver file for the frozen source named in
the manifest. The archived version-0.1.0 initialization rerun and the current
package's numerical-oracle tests separately validate the public RGD/RGN API and
convergence recording; comparator implementations are not installed.

Absolute wall-clock seconds are reference-machine measurements. Portable
runtime validation uses only the separated tier order
`{RGD, RGN} < CP-ALS < HOPLS < TT-IHT < {RRR, Tucker-RGD}`; neither exact
seconds nor order within the two near-tied pairs is a pass/fail condition.
