# Generic GMM Objective Primitives

Status: internal GMM moment-objective evaluation added; validation is ongoing.

Given a moment vector g(θ), symmetric weighting matrix W, and caller-supplied moment Jacobian J, the evaluator computes Q=gᵀWg, ∇Q=2JᵀWg, and optional Gauss-Newton curvature 2JᵀWJ. This supports one-step/two-step GMM assembly without embedding a model-specific estimator.

The weighting matrix must be symmetric. Positive semidefiniteness is expected from the caller; materially negative objective values are rejected. Outputs are committed only after all requested quantities have been computed and checked for finiteness. Accumulations use compensated summation and dimensions are checked before allocation.

The layer also provides the Hansen J overidentification log-p-value using chi-square degrees of freedom equal to moment count minus parameter count. Native tests cover hand-computable objective, gradient, curvature, multivariate weighting, symmetry validation, overidentification inference, and output preservation.
