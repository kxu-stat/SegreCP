# Tensor-on-tensor regression simulation records

This bundle contains replication-level results for the paired seven-method
comparison: CP-ALS, RRR, Tucker-RGD, HOPLS, TT-IHT, RGD, and RGN. Each method
in a cell and replication uses the same generated coefficient, training and
test samples, and noise realization. The coefficient has five modes of size
10, with three predictor and two response modes; the CP rank is 3. The design
uses condition number 30, 1,000 training and 1,000 test observations, four
coherence levels (`0, 0.2, 0.5, 0.75`), three noise levels (`0, 0.5, 1`),
20 replications per cell, and ten updates.

`final_seven_method_grid_raw.csv` holds the 1,680 method-cell-replication
rows. The other `final_seven_method_*summary.csv` files contain cell and
replication-balanced summaries. `checkpoint_metrics_raw.csv` records
iterations 3, 5, and 10 for the six methods with comparable outer-update
checkpoints. The timing files contain a separate serial run; exact seconds
depend on the machine. `final_settings.csv` records the experiment design.

The package's tensor-regression results vignette and
`inst/reproduce/validate_rebuttal_tables.R` recompute tables and plots from
these records. To generate new data and refit the authors' RGD/RGN methods,
use the source-checkout script in `reproducibility/simulations/`.
