# Delta-Method and Wald Inference Primitives

Status: initial internal generic inference primitives added; validation is ongoing.

## Delta method

The delta-method variance primitive computes gᵀVg from a caller-supplied gradient and covariance matrix. A standard-error wrapper returns its square root. Accumulation uses compensated summation; materially negative variances are rejected rather than silently clipped.

## Wald tests

The generic Wald test evaluates Rθ = r using W = dᵀ(RVRᵀ)⁻¹d, d=Rθ−r. It forms the restricted covariance and solves the positive-definite system with Cholesky forward/back substitution; it does not explicitly invert matrices. The p-value is returned in log form using the chi-square log-survival primitive.

The caller supplies the estimate, covariance, restriction matrix, and null vector. No model-specific fitting API is introduced. Singular restricted covariance is reported distinctly.

## Validation

Native tests cover scalar and joint restrictions, a hand-computable delta-method variance, chi-square log-tail integration, singular covariance, negative variance rejection, and output preservation.
