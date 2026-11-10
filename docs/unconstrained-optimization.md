# Unconstrained Optimization

Status: internal BFGS and L-BFGS implementations with Armijo backtracking; validation is ongoing.

## Methods

### BFGS

The existing dense BFGS implementation maintains an inverse-Hessian approximation. Its workspace is O(n²) and is appropriate for moderate parameter dimensions.

### L-BFGS

The limited-memory implementation retains a configurable history of recent step/gradient-difference pairs and applies the standard two-loop recursion. Its workspace is O(nm), where n is the parameter count and m is the history size. No dense Hessian or inverse-Hessian matrix is allocated.

Both methods share the same objective/gradient callbacks, central finite-difference fallback, Armijo backtracking, convergence tolerances, result ownership, and termination diagnostics. L-BFGS rejects zero history size and checks all history-size multiplications before allocation. Curvature pairs with insufficient positive curvature are skipped rather than admitted into the recursion.

## Validation

Native tests cover analytic-gradient Rosenbrock convergence, a 4096-parameter quadratic with a small L-BFGS history, invalid history sizes, size-overflow rejection, and allocation-failure cleanup at the history workspace boundary. Existing BFGS tests remain unchanged.
