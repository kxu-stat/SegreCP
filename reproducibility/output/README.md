# Stored results

Results of the runs reported in the paper. See `../README.md` for the commands
that write them.

| File | Content |
|---|---|
| `simulations/results_<study>.csv.gz` | one row per case (coherence, snr, replication), method and iteration: relative Frobenius error and objective; for `timing_*`, seconds of one fit |
| `figures/regression_cp.pdf`, `regression_other_formats.pdf` | Figures 1 and 2 |
| `figures/decomposition.pdf`, `decomposition_gauss_newton.pdf` | Figures 3 and 4 |
| `figures/profile_regression.pdf`, `profile_decomposition.pdf`, `tables/profile_*.csv` | Figure 5 |
| `tables/table_timing.tex`, `timing_summary_*.csv` | Table 3 (median seconds per fit and per iteration) |
| `tables/table_time_regression.tex`, `table_time_decomposition.tex`, `time_to_level_*.csv` | Tables 4 and 5 |
| `tables/table_errors_regression.tex`, `table_errors_decomposition.tex` | Tables 6 and 7 |
| `tables/table_initialization_regression.tex`, `table_initialization_decomposition.tex`, `errors_initialization.csv` | Tables 8 and 9 |
| `tables/rmse_curves_*.csv` | the curves of Figures 1–4 |
| `real_data/accuracy.csv` | test RMSE, test R² and training RMSE per fold, method and checkpoint (Tables 1 and 10) |
| `real_data/session.txt` | CPU and R session of the trade-data run |
