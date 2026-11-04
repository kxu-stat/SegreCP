# SegreCP 0.1.5

- Reuse batched tensor contractions and refresh decomposition caches in one
  pass, preserving sequential component updates and exact T-HOSVD retractions.
- Avoid unused cross-moments and full tangent Gram matrices in regression,
  and reuse component predictions when the normalization permits.
- Set the default ridge to zero in both control profiles. Explicit positive
  ridge values remain available for regularized RGN.
- Extend tests for general tensor orders, CP ranks, unequal mode dimensions,
  initializer routes, and regression normalization conventions.
- Align fresh simulation and trade-data runs with the current zero-ridge
  experiments, including the normalized trade-data CPCA start. Preserve the
  bundled rebuttal results as historical records.
