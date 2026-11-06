# SegreCP 0.1.5

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
