# Reproducing the results of the paper

This folder contains the code that produces every figure and table of
*Guaranteed Noisy CP Tensor Recovery via Riemannian Optimization on the Segre
Manifold* (NeurIPS 2026): the simulations of Section 5 and Appendix C and the
trade-data analysis of Section 5 and Appendix D. RGD and RGN are run with the
`SegreCP` package in this repository; the competing methods are implemented in
`code/`.

## Contents

| Path | Content |
|---|---|
| `run_all.sh` | runs every step below in order |
| `code/simulations.R` | data generation for the regression simulations |
| `code/regression_adapters.R` | regression methods: CP-ALS, RGD, RGN, RRR, Tucker-RGD, TT-IHT, HOPLS, and the shared start of CP-ALS, RGD and RGN |
| `code/tucker_rgn.R` | Tucker-RGN for tensor-on-tensor regression |
| `code/decomposition_comparators.R` | decomposition methods: Orth-ALS, ICO, NLS-GN, GN-CG, CP-ALS+ELS, hybrid AMDM |
| `simulations/run_simulations.R` | simulation studies (one case per coherence, noise level and replication) |
| `simulations/plot_figures.R` | Figures 1–5 and Tables 3–9 from the simulation results |
| `real_data/` | trade data: download and preparation (`1_download_sources.R`, `2_download_tariffs.R`, `3_build_data.R`), fits (`run_real_data.R`) and Tables 1 and 10 (`make_tables.R`); see `real_data/README.md` |
| `output/` | results of the runs reported in the paper (see below) |

## Software and hardware

R 4.3.3 with `SegreCP` (this repository), `rTensor` 1.5.0,
`MultiwayRegression` 1.2, `Rcpp` 1.0.12, `ggplot2` 3.4.4 and `scales` 1.3.0;
the real-data preparation also uses `data.table`, `jsonlite` and `xml2`. All
reported numbers were computed on one machine with an Intel Xeon 2.1 GHz CPU
and 8 GB of memory, one fit per core and one BLAS thread.

```sh
Rscript -e 'install.packages(c("rTensor", "MultiwayRegression", "ggplot2", "scales"))'
R CMD INSTALL .        # SegreCP, from the repository root
```

## Running

All commands are run from the repository root; all paths are relative to it.
`bash reproducibility/run_all.sh` runs the steps in this order:

| Step | Command | Paper | Time on one core |
|---|---|---|---|
| 1 | `run_simulations.R regression <out> 1:20` | Figures 1, 2, 5; Tables 4, 6 | 2 h |
| 2 | `run_simulations.R decomposition <out> 1:20` | Figures 3, 4, 5; Tables 5, 7 | 2 h |
| 3 | `run_simulations.R decomposition_long <out> 1:20` | Figures 3–5, Table 5 (ρ = 0.75) | 2.5 h |
| 4 | `run_simulations.R initialization <out> 1:20` | Tables 8, 9 | 1.5 h |
| 5 | `run_simulations.R timing_regression <out> 1:5` | Tables 3, 4; Figure 5 | 0.5 h |
| 6 | `run_simulations.R timing_decomposition <out> 1:5` | Tables 3, 5; Figure 5 | 0.5 h |
| 7 | `plot_figures.R <out> <figures> <tables> paper` | Figures 1–5; Tables 3–9 | 1 min |
| 8 | `real_data/1_download_sources.R <data>` and `2_download_tariffs.R <data>` | trade data sources | download |
| 9 | `real_data/3_build_data.R <data>` | prepared trade data | 1 min |
| 10 | `real_data/run_real_data.R <data>/prepared <out>` | Tables 1, 10 | 10 min |
| 11 | `real_data/make_tables.R <out>` | Tables 1, 10 | seconds |

Steps 1–4 can run cases in parallel (`SEGRECP_CORES=<cores>`). The timing
studies (steps 5–6) must run on one core with one BLAS thread
(`SEGRECP_CORES=1 OPENBLAS_NUM_THREADS=1`) and nothing else running. A run
that stops can be restarted with the same command: finished cases are kept in
`<out>/cases_<study>/` and skipped.

To rebuild the figures and tables of the simulations without rerunning them,
use the stored results:

```sh
Rscript reproducibility/simulations/plot_figures.R reproducibility/output/simulations \
  reproducibility/output/figures reproducibility/output/tables paper
```

## Stored results (`output/`)

- `simulations/results_<study>.csv.gz`: one row per case, method and
  iteration (relative Frobenius error and objective; for the timing studies,
  seconds per fit), as written by `run_simulations.R` for the paper.
- `figures/`, `tables/`: Figures 1–5 and Tables 3–9 as made by
  `plot_figures.R` from these results, with the CSV summaries behind them.
- `real_data/accuracy.csv`, `session.txt`: the trade-data run behind Tables 1
  and 10. The prepared trade data are in the supplementary material of the
  paper (see `real_data/README.md`).

Rerunning a case reproduces the stored errors up to floating-point rounding,
which is largest for CP-ALS+ELS (relative differences below 1e-7).

## Simulation design

- Settings: coherence ρ ∈ {0, 0.25, 0.5, 0.75} times signal-to-noise ratio of
  the weakest component s ∈ {∞, 30, 10, 3}; 20 replications per setting; rank
  3; condition number κ = 30; 20 iterations (100 in `decomposition_long`).
  With df = 1 + Σ(p_ℓ − 1), s = λ_r / (σ √df) in decomposition and
  s = λ_r / (σ √(df/n)) in regression.
- Factors: in each mode, the Cholesky factor of the Gram matrix with entries
  ρ^|i−j|, embedded in R^p and rotated by a random orthogonal matrix.
- Decomposition: 60 × 30 × 30 × 30, λ_i = √r κ^(1/2−(i−1)/(r−1)) · 60, seed
  20261104 + 4000 + 100 k + replication (k = index of ρ).
- Regression: X_m ∈ R^(10×10×10), Y_m ∈ R^(10×10), n = 1000, design rows
  vec(X_m)/√n with standard normal entries (σ_X = 1), λ_i =
  κ^(1/2−(i−1)/(r−1)), seed 20261204 + 100 k + replication.
- Starts: CP-ALS, RGD and RGN share the CPCA start (`SegreCP::cp_initialize`,
  in regression applied to X'Y with the input and output modes as the two
  groups of the unfolding). The other methods use their own starts; RRR starts
  at random.
- RGD and RGN: `SegreCP` with ridge 0 and no early stopping. Decomposition:
  RGD step 0.5, components updated in a fixed order. Regression: RGD step
  α = 1 = σ_X^(−2), components updated in decreasing order of the current
  weights. RGN for decomposition is RGD with step 1.
- Initialization study: the settings with noise, starts CPCA, T-HOSVD (the
  i-th component from the i-th leading singular vectors of the mode unfoldings
  of Y, or of X'Y in regression) and random unit factors
  (`cp_initialize(method = "random")`, weights fitted to Y or X'Y, seed = data
  seed).

## Summaries in the figures and tables

- Curves: root mean square over the 20 replications of the relative Frobenius
  error ‖T̂ − T*‖_F / ‖T*‖_F, on a log10 scale; values below 1e-15 are drawn
  at 1e-15 and values above 1 (the random start of RRR) are not drawn.
- Tables 6 and 7: mean (standard deviation) over the replications of the
  relative error at iteration 20, times 10³, in the settings with noise; bold
  is the lowest mean of a setting at the printed precision. Tables 8 and 9:
  mean of the relative error, times 10³, at iterations 0, 3, 5, 10 and 20
  from each start; iteration 0 is the error of the start.
- RRR stops before iteration 20 in 17 of the 320 regression paths, when its
  update no longer changes the fit. Such a path keeps its last estimate for the
  remaining iterations (`fill_paths` in `plot_figures.R`).
- Table 3: median over all settings of the time of 20 iterations divided by
  20, with the error evaluation switched off (the estimates are unchanged), one
  fit at a time, one warm-up fit per method and the method order rotated. The
  shared CPCA start is timed separately.
- Tables 4 and 5: the error level of a setting is the lowest root mean square
  error among CP-ALS, RGD and RGN at its last iteration (20, or 100 at
  ρ = 0.75 in decomposition). The entry is the first iteration at which a
  method comes within 1% of the level, times its median time per iteration in
  that setting.
- Figure 5: performance profiles (Dolan and Moré, 2002) over the 12 settings
  with noise. For each method and factor τ, the fraction of settings in which
  its time to the error level (left) or its root mean square error at the last
  iteration (right) is at most τ times the best value of the setting; a method
  that does not reach the level of a setting is not within any factor there.

## Competing methods

- CP-ALS (decomposition): the update of the `cp` function of `rTensor` (same
  unfolding, Khatri–Rao product and normal equations), without its convergence
  check. CP-ALS (regression): alternating least squares over the input and
  output factors.
- Orth-ALS: orthogonalization of the factors before each ALS step, following
  the MATLAB code of Sharan and Valiant (2017).
- ICO: iterative concurrent orthogonalization (Han and Zhang, 2023).
- NLS-GN and GN-CG: Gauss–Newton with at most 20 conjugate-gradient inner
  iterations (Sorber et al., 2013; Singh et al., 2021).
- CP-ALS+ELS: ALS with enhanced line search (Rajih et al., 2008).
- Hybrid AMDM: all r singular values are inverted in iterations 1–10, one
  fewer in each following block of 10 iterations (Singh and Solomonik, 2023),
  so that it equals CP-ALS from iteration 31.
- RRR: `MultiwayRegression::rrr` (Lock, 2018) with its own random start. Its
  outer loop is given an iteration limit and records the coefficient after
  every outer iteration (`.make_rrr_fixed` in `regression_adapters.R`).
- Tucker-RGD: Riemannian gradient step with exact line search and backtracking,
  and T-HOSVD retraction (Luo and Zhang, 2024). Tucker-RGN: port of
  `RGN_tenreg` (HOSVD retraction) from the code of Luo and Zhang (2023),
  extended to tensor responses, with the normal equations formed from the
  Kronecker structure. TT-IHT: gradient step with step 1 = σ_X^(−2) (real
  data: 1/λ_max(X'X)) followed by TT-SVD (Qin et al., 2024). HOPLS: three
  latent blocks (Zhao et al., 2013). Multilinear and TT ranks are 3.

## Implementation notes on SegreCP

- Components are updated one at a time (Gauss–Seidel); each update uses the
  current values of the components already updated in the same iteration.
- The rank-(1, …, 1) T-HOSVD retraction of a rank-one tensor plus a tangent
  vector is computed from one 2 × 2 eigenproblem per mode, without forming the
  tensor.
- Regression stores the output modes first internally and permutes the
  factors back on output. The tangent design of the input modes is formed by
  modewise contractions of the design; on the output side an orthonormal basis
  of the output tangent space (p_out × df_out) is used. Neither the ambient
  tangent basis (p × df) nor the full design is formed.
