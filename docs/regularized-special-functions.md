# Regularized Incomplete Special Functions

Status: initial internal C implementation added; numerical validation is ongoing.

## Scope

The internal special-function layer currently provides:
- Regularized lower incomplete gamma, P(a,x).
- Regularized upper incomplete gamma, Q(a,x).
- Regularized incomplete beta, I_x(a,b).
- Inverses of P(a,x), Q(a,x), and I_x(a,b), solved with monotone bracketing and bisection.
- Digamma and trigamma for positive finite real arguments, using recurrence and asymptotic Bernoulli expansions.

These are internal numerical primitives. They are not yet part of the installed public C ABI or language-binding API.

## Numerical methods

### Incomplete gamma

For x < a + 1, P(a,x) is evaluated with its convergent power series. For x >= a + 1, Q(a,x) is evaluated with a modified continued fraction using the modified Lentz method. The opposite tail is obtained by complementing the directly evaluated tail.

This branch selection computes the small tail directly rather than subtracting it from one. At x = 0, the exact endpoint values are P(a,0)=0 and Q(a,0)=1.

### Incomplete beta

The regularized incomplete beta uses a modified Lentz continued fraction and the standard symmetry transformation around (a+1)/(a+b+2). The exponential prefactor is assembled in the log domain using lgamma, log, and log1p; it is not formed as a product of potentially underflowing intermediate factors.

At x=0 and x=1, exact endpoint values are returned.

## Digamma and trigamma

The positive-argument digamma and trigamma functions shift their arguments by recurrence to x >= 8, then evaluate asymptotic Bernoulli expansions. Restricting the current contract to x > 0 avoids reflection near non-positive poles. Non-finite intermediate results are reported as numerical failures and do not modify the output.

## Inverse functions

The inverse gamma routines bracket the nonnegative root by expanding the upper endpoint, then bisect while evaluating the requested tail directly (the Q inverse does not convert its target to 1-p). The inverse beta routine bisects on [0,1] and reflects upper-tail probabilities through I_x(a,b)=1-I_{1-x}(b,a). The algorithms stop when the floating-point midpoint can no longer refine the bracket. These are robust baseline inverses; safeguarded Newton acceleration and high-precision tail sweeps remain future work.

## Contracts

- All shape parameters must be finite and strictly positive.
- Gamma x must be finite and nonnegative.
- Beta x must be finite and lie in [0,1].
- A required output pointer must be non-NULL.
- Outputs are written only after successful completion.
- Invalid arguments, non-finite inputs, domain errors, and numerical non-convergence are distinct statuses.
- Results are constrained to [0,1]; only tiny floating-point excursions attributable to rounding are clamped.
- Iterative algorithms have a finite iteration cap and report numerical failure if convergence is not achieved.

## Validation

Independent high-precision reference vectors are embedded as constants in the native tests (generated offline with arbitrary-precision arithmetic; no runtime dependency is introduced). Coverage includes small shape parameters, central regions, large shapes, a gamma upper tail near 1e-20, monotonicity, and the digamma derivative identity psi\u2032(x)=trigamma(x).

Native tests cover:
- Closed-form references for gamma shape 1 and shape 1/2.
- A polynomial reference value for I_0.2(2,3).
- Gamma lower/upper tail complementarity.
- Small-x and far-upper-tail gamma values.
- Beta symmetry and exact endpoints.
- Round trips for inverse gamma P/Q and inverse beta over lower, central, and upper probabilities.
- Inverse probability-domain validation and output preservation.
- Digamma/trigamma closed-form references at 1, 1/2, and 8, recurrence identities, and invalid-input handling.
- Invalid domains, non-finite inputs, and output preservation.

These tests are regression coverage, not a substitute for broad high-precision reference sweeps. Before exposing these functions publicly, add independently generated high-precision vectors across logarithmically spaced shape and argument ranges, monotonicity tests, derivative identities, and stress tests for extreme finite parameters.

## Deferred work

- Digamma and trigamma.
- Inverse regularized gamma and beta functions.
- Error function inverses and normal quantile.
- High-precision reference-vector generation and broad parameter sweeps.
- Public API exposure after numerical contracts and compatibility review.
