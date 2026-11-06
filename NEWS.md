# SegreCP (development version)

- Reuse batched tensor contractions and refresh decomposition caches in one
  pass, preserving sequential component updates and exact T-HOSVD retractions.
- Avoid unused cross-moments and full tangent Gram matrices in regression,
  and reuse component predictions when the normalization permits.
- Set the default ridge to zero in both control profiles. Explicit positive
  ridge values remain available for regularized RGN.
- The default regression RGD step size is now `1 / design_variance`
  (σ_X^{-2}), the fixed step of the paper's regression convergence result; the
  decomposition default remains 0.5.
- Extend tests for general tensor orders, CP ranks, unequal mode dimensions,
  initializer routes, and regression normalization conventions.
- `reproducibility/` now contains the code for the simulations, figures,
  tables and trade-data analysis of the NeurIPS 2026 paper. The archived
  results of earlier experiments, their vignettes and tests were removed.
- Rewritten help pages: each exported function describes its model,
  arguments and returned object and has runnable examples. `print()` and
  `summary()` of a fit show the dimensions, the stopping reason and the
  component weights; `trace = TRUE` reports per iteration.
- The tangent basis of the output modes in tensor-on-tensor regression is
  built without a loop over its columns, with identical results.
- `init_args` is now an error when `init` is a list of components (it was
  ignored).
- The vignette `algorithm-fidelity` is now `algorithms`.
