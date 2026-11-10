# Linear Programming

Status: internal two-phase simplex implementation added; numerical validation is ongoing.

## Problem form

The solver maximizes cᵀx subject to Ax ≤ b and x ≥ 0. Inputs are dense row-major matrices. Equality constraints and finite upper bounds can be represented by pairs of inequalities. No public ABI exposure is made yet.

## Algorithm

A two-phase tableau simplex first finds a feasible basis when any right-hand side is negative, then optimizes the objective. Bland's variable-index selection is used for entering and ratio ties to reduce cycling risk. The implementation distinguishes optimal, infeasible, unbounded, iteration-limit, and numerical-failure outcomes.

## Contracts and limits

- Inputs and options must be finite; variable count must be positive.
- Zero constraints are supported.
- Solution/objective outputs are written only for an optimal result.
- Tableau and work-array sizes are checked before allocation.
- The current implementation uses a dense tableau and is intended as a correctness baseline for small/medium LPs. Sparse revised simplex, scaling, presolve, and basis warm-start are future work.

## Validation

Native tests cover a reference bounded LP, a negative-RHS Phase-I case, infeasibility, unboundedness, zero-constraint bounded/unbounded cases, a degenerate LP, iteration limit, non-finite input, invalid options, and output preservation.
