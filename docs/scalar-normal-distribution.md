# Scalar Normal Distribution

Status: internal normal-distribution primitives added; validation is ongoing.

## Contract

The internal scalar normal family provides log-PDF, PDF, CDF, survival function, log-CDF, log-survival, and quantile for finite location and strictly positive finite scale.

- Inputs are scalar and no memory is allocated.
- The result pointer is required and remains unchanged on failure.
- Quantile probabilities include endpoints 0 and 1, which map to negative and positive infinity.
- Density underflow is represented by zero; log-density may be negative infinity for extreme finite inputs.
- The normal CDF and survival tails are evaluated through the regularized upper incomplete gamma function. Log tails are formed directly from a log-domain incomplete-gamma upper tail and log1p, not by taking log(1-CDF); tails remain representable in log form even after the probability itself underflows.

## Validation

Native tests cover standard-normal reference values, location/scale transformation, symmetry, log-tail behavior at 10 and 40 standard deviations, quantiles, invalid parameters, non-finite values, and extreme finite location differences.
