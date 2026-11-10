# Scalar Distributions: Normal, Student-t, and Chi-square

Status: internal normal, Student-t, and chi-square scalar distribution primitives added; validation is ongoing.

## Contract

The internal scalar normal family provides log-PDF, PDF, CDF, survival function, log-CDF, log-survival, and quantile for finite location and strictly positive finite scale.

- Inputs are scalar and no memory is allocated.
- The result pointer is required and remains unchanged on failure.
- Quantile probabilities include endpoints 0 and 1, which map to negative and positive infinity.
- Density underflow is represented by zero; log-density may be negative infinity for extreme finite inputs.
- The normal CDF and survival tails are evaluated through the regularized upper incomplete gamma function. Log tails are formed directly from a log-domain incomplete-gamma upper tail and log1p, not by taking log(1-CDF); tails remain representable in log form even after the probability itself underflows.

## Student-t and chi-square

Student-t provides log-PDF/PDF, CDF/survival and their log tails, and quantiles. Its tails are evaluated through the log-lower regularized beta primitive, and quantiles use the inverse regularized beta. Chi-square provides the same scalar family using regularized gamma P/Q and the log-Q primitive. Support endpoints are explicit: chi-square has zero density outside its nonnegative support, and quantiles at probability endpoints return their mathematical limits.

## Discrete Poisson and Binomial

Poisson and Binomial provide PMF/log-PMF, CDF/survival, log tails, and integer quantiles. Poisson tails use log-domain regularized gamma P/Q; Binomial tails use log-domain regularized beta lower/upper tails. Count parameters are limited to at most 2^53−1 so k+1 and factorial arguments remain exactly representable in double-precision special-function parameters. The Poisson quantile at probability one saturates to UINT64_MAX because its support is unbounded.

## F distribution

F distribution uses the standard numerator/denominator degrees-of-freedom parameterization. Its density is evaluated in the log domain; CDF and survival tails use the regularized beta lower/upper logs, and quantiles use the inverse regularized beta with a log-domain odds transform.

## Gamma and beta

Gamma distribution uses shape/scale parameterization and regularized gamma P/Q tails; beta uses alpha/beta parameters and regularized incomplete beta log tails. Both expose density/log-density, CDF/survival, log tails, and quantiles. Support boundaries are explicit, and CDF/survival use complementary special functions rather than subtracting a rounded CDF.

## Validation

Native tests cover standard-normal reference values, location/scale transformation, symmetry, log-tail behavior at 10 and 40 standard deviations, quantiles, invalid parameters, non-finite values, and extreme finite location differences. Student-t and chi-square tests cover closed-form reference cases, tails, quantiles, and parameter validation. Gamma and beta tests cover reference values, quantiles, tail logs, support boundaries, and invalid parameters. F-distribution tests cover closed-form reference cases, upper-tail logs, quantiles, and validation. Poisson and Binomial tests cover PMF/CDF/survival references, quantiles, degenerate endpoints, and invalid counts/probabilities.
