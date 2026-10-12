# Initialization and convergence records

These replication-level records support the package's initialization and
convergence vignette.

- `decomp_ref/`: 20 paired replications of a `60 x 30 x 30`, rank-three noisy
  CP decomposition, with CPCA and HOSVD starts and ten complete sweeps.
- `reg_ref/`: four tensor-regression settings comparing CPCA and HOSVD starts.
- `decomp_objective/`: observed-objective trajectories for the decomposition
  convergence plot, with their validation checks.

The convergence vignette shows the observed least-squares objective and error
against the simulated truth. The latter may fluctuate near the noise floor.

For new fits, use `segre_decompose()` or `segre_regress()` with
`init = "cpca"` or `"hosvd"`.
