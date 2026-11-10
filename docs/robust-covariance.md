# Robust Sandwich Covariance Primitives

Status: initial internal HC, cluster-robust, and Newey-West covariance primitives added; validation is ongoing.

## Generic sandwich

Given a caller-supplied bread matrix B and row-major observation score vectors sᵢ, the base estimator computes B (Σ sᵢsᵢᵀ) Bᵀ. HC0 is uncorrected; HC1 applies n/(n−p) and requires n>p. The caller is responsible for constructing scores appropriate to the estimating equations; no regression-model-specific API is introduced.

## Cluster robust

Cluster scores are summed within dense caller-provided cluster IDs, then outer products of cluster sums form the meat. CR0 is uncorrected; CR1 applies G/(G−1). Cluster IDs must lie in [0,G).

## HAC

Newey-West HAC uses Bartlett lag weights and includes both lagged cross-products and their transposes. Bandwidth is clamped to n−1.

## Numerical and ownership contracts

Inputs must be finite; dimensions and allocations are checked. Kahan compensated summation is used for meat accumulation. The final covariance is symmetrized to remove floating-point asymmetry. Output is written only after the complete result is finite. Dense workspaces are O(p²), with O(Gp) additional storage for cluster aggregation.

## Validation

Native tests cover hand-computable HC0/HC1, CR0/CR1, and one-lag Newey-West results, invalid cluster labels, HC1 degrees-of-freedom requirements, non-finite scores, and output preservation.
