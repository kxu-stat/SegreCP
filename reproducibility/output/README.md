# Stored results

Results of the runs reported in the paper. See `../README.md` for the commands
that write them.

| File | Content |
|---|---|
| `simulations/results_<study>.csv.gz` | one row per case (coherence, snr, replication), method and iteration: relative Frobenius error and objective; for `timing_*`, seconds of one fit |
| `figures/*.pdf` | Figures 1–4 |
| `tables/table_errors_regression.tex`, `table_errors_decomposition.tex` | Tables 3 and 4 |
| `tables/timing_summary_*.csv` | Table 5 (median seconds per iteration) |
| `tables/table_time_to_level.tex`, `time_to_level_*.csv` | Table 6 |
| `tables/table_initialization.tex`, `errors_initialization.csv` | Table 7 |
| `tables/rmse_curves_*.csv` | the curves of Figures 1–4 |
| `real_data/accuracy.csv` | test RMSE, test R² and training RMSE per fold, method and checkpoint (Tables 1 and 8) |
| `real_data/session.txt` | CPU and R session of the trade-data run |
