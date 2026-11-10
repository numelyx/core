# Box-Constrained Convex Quadratic Programming

Status: internal projected-FISTA implementation added; numerical validation is ongoing.

## Problem form

Minimize 0.5 xᵀHx + gᵀx subject to lower ≤ x ≤ upper, where H must be symmetric positive semidefinite. Infinite bounds are supported; NaN bounds are rejected.

## Algorithm

The implementation validates symmetry and positive semidefiniteness with a tolerance-aware Cholesky-style factorization. It uses the maximum absolute row sum as a safe Lipschitz upper bound for the gradient and applies projected FISTA with a monotone restart when extrapolation increases the objective. The projected-gradient mapping provides the stationarity criterion.

A zero Hessian is handled as a linear program over the box, including exact detection of unbounded coordinates. For nonzero semidefinite Hessians, the current solver reports iteration-limit or numerical failure if it cannot establish convergence; general recession-cone unboundedness detection is not yet implemented.

## Limits

This is a box-constrained QP primitive, not a general linearly constrained QP solver. General inequalities, equalities, sparse Hessians, active-set methods, and primal-dual interior-point QP remain future work.

## Validation

Native tests cover positive-definite and rank-deficient semidefinite Hessians, active bounds, linear bounded/unbounded cases, non-PSD rejection, and output preservation.
