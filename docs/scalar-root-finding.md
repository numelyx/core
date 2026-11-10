# Scalar Root-Finding Primitives

Status: internal bracketed bisection implementation added; validation is ongoing.

## Contract

`numelyx_core_bisect` solves a scalar continuous-root problem on a caller-supplied finite interval. The callback receives an x value, immutable context, and output pointer. It must return true only when it writes a complete function value.

- Endpoints must be finite and strictly ordered.
- Endpoint roots are accepted without iteration.
- Otherwise, endpoint function values must have opposite signs; sign bits are compared without multiplying values.
- Every callback value must be finite.
- The solver allocates no memory and does not retain callback context.
- On callback/domain/bracketing failure, the result structure remains unchanged.
- On iteration-limit termination, the best current midpoint and diagnostics are returned with a distinct status.
- Absolute and relative interval tolerances are both supported; at least one must be positive.

## Numerical details

The midpoint calculation avoids overflow for intervals crossing zero. Convergence is tested using half-width arithmetic, which also avoids forming an overflowing full interval width. If floating-point rounding makes the midpoint equal an endpoint, the endpoint with the smaller residual is returned.

## Validation

Native tests cover sqrt(2), endpoint roots, brackets spanning extreme finite values, no-bracket behavior, callback failure, non-finite callback values, iteration-limit reporting, and output preservation.

## Safeguarded Newton

The safeguarded Newton solver maintains the same sign-changing bracket while attempting a Newton step from the current iterate. A zero derivative, non-finite candidate, or candidate outside the open bracket falls back to the safe midpoint. Function and derivative callbacks are counted separately. Callback failures and non-finite derivatives are reported distinctly and leave the result structure unchanged.
