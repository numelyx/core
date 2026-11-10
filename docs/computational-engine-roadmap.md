# Numelyx Computational Engine Roadmap

## Mission

Numelyx Core is a language-independent, composable numerical and computational engine. It must provide reusable, numerically reliable primitives from which higher-level statistical, econometric, optimization, mathematical-programming, and fuzzy-computing systems can be assembled.

Core must not become a catalogue of model-specific fitting functions. Formula parsing, data-frame semantics, model orchestration, and language bindings belong in higher layers. Core owns numerical contracts, reusable data structures, solvers, algorithms, and diagnostics.

This roadmap is a capability and verification contract, not a claim that the listed capabilities already exist.

## Architectural rules

1. Keep the public C API, private C implementation, and language bindings separate.
2. Prefer opaque handles and explicit ownership; document borrowed versus owned inputs and outputs.
3. Make randomness explicit through caller-owned RNG streams. Define deterministic stream splitting before parallel stochastic algorithms are added.
4. Use checked size arithmetic, finite-input validation, stable algorithms, and cleanup-safe failure paths.
5. Prefer factorizations and linear solves over explicit matrix inverses.
6. Support dense, sparse, banded, triangular, symmetric, block, and structured operators where algorithms benefit; do not force every problem into dense storage.
7. Expose callbacks and composable primitives, not model-specific APIs such as fit_ols or fit_probit.
8. Keep third-party implementation details and types out of the public ABI. Review license, maintenance, portability, performance, and fallback behavior before adding dependencies.
9. Every public capability requires independent numerical reference tests, adversarial edge tests, ownership tests, and API documentation.
10. Do not declare a subsystem complete until its acceptance criteria and integration tests pass on the supported build matrix.

## Existing verified foundation

The current public surface includes:
- Opaque Matrix handles, composable creation, views, arithmetic, multiplication, materialization, and structured representations.
- LU, QR, reduced SVD, Cholesky, LDLT, symmetric eigen decomposition, least-squares solves, pseudoinverse, rank, norms, condition estimates, and scalar reductions.
- Explicit-state PCG32 RNG, uniform and normal variates, cloning, and index sampling with or without replacement.
- Callback-driven unconstrained BFGS with Armijo backtracking and optional central finite-difference gradients.
- Multivariate Gaussian log density and sampling.

The current public probability API is not a general distribution library. The current optimizer is not a general constrained or mathematical-programming solver. The capabilities below remain roadmap items until implemented and tested.

## Implementation status on this branch

The following reusable internal primitives have been added. They are not yet a complete public numerical API:
- Stable log-domain scalar operations, including log1mexp and log-difference-exp.
- Regularized incomplete gamma P/Q, their log tails, regularized incomplete beta and its log tails, plus bracketed inverse gamma/beta functions.
- Positive-argument digamma/trigamma and erf/erfc/inverse-erf/normal-quantile support.
- Scalar Normal, Student-t, chi-square, Gamma, Beta, F, Uniform, Logistic, Weibull, and Lognormal distributions with density/log-density, CDF/survival, log tails, and quantiles.
- Discrete Poisson and Binomial distributions with PMF/log-PMF, CDF/survival, log tails, and integer quantiles.
- Generic HC0/HC1, cluster CR0/CR1, and Bartlett Newey-West HAC sandwich covariance primitives from caller-supplied score arrays.
- Generic delta-method variance/standard-error and Cholesky-based Wald tests with log chi-square p-values.
- Generic GMM moment objective, Jacobian gradient, and Gauss-Newton curvature evaluation. Counts are capped at 2^53−1 where exact conversion to the double-precision special-function layer is required.
- Atomic affine Box-Muller sampling against a caller-owned RNG stream.
- Dense BFGS and O(nm)-memory L-BFGS, plus safeguarded scalar bisection/Newton root solvers.
- Dense two-phase simplex for max cᵀx subject to Ax ≤ b and x ≥ 0, with Phase-I feasibility and Bland anti-cycling selection.
- Box-constrained convex QP using PSD validation and projected FISTA. General linear constraints and recession-cone unboundedness detection are not yet implemented.
- Fuzzy membership functions, t-norms/s-norms, complement, alpha cuts, interval arithmetic, sampled-set defuzzification, and Mamdani/Sugeno aggregation.

Native regression tests and focused subsystem documentation have been added. Each subsystem remains internal pending broader independent reference sweeps, stress testing, performance work, and public API review. CI status must be checked against the latest branch commit; earlier successful runs do not validate later changes.

## Capability domains

### 1. Numerical foundations and special functions

Provide stable scalar and vector primitives:
- log-gamma, log-beta, erf, erfc, normal-tail and inverse-normal functions.
- Regularized incomplete gamma and beta functions and their complements.
- Stable log1p, expm1, log1pexp, logsumexp, logdiffexp, and log-domain probability helpers.
- Compensated summation, scaled sum-of-squares, stable norms, online mean/variance/covariance, and mergeable accumulators.
- Root finding: bracketed bisection, Brent-style methods, safeguarded Newton, and secant methods.
- Explicit tolerance, convergence, non-finite, domain, and overflow contracts.

### 2. Probability distributions

A consistent distribution contract should cover, when mathematically applicable:
- PDF/PMF and log-PDF/log-PMF.
- CDF, survival function, log-CDF, and log-survival.
- Quantile/inverse CDF.
- Sampling from an explicit RNG.
- Parameter validation, support boundaries, endpoint behavior, and status reporting.

Initial scalar families:
- Normal, Student-t, chi-square, F, uniform, exponential, gamma, beta, logistic, log-normal, Weibull.
- Bernoulli, binomial, Poisson, geometric, and negative-binomial.

Multivariate and structured families:
- Multivariate normal and multivariate Student-t.
- Categorical, multinomial, Dirichlet, and Wishart where stable reusable primitives justify them.

Tail probabilities must be computed directly or through complementary special functions; log(1 - cdf) is not an acceptable general implementation of log-survival.

### 3. Random sampling and resampling

- Weighted categorical sampling with validated, numerically stable weights.
- Stratified, grouped, cluster, block, moving-block, and permutation sampling.
- Bootstrap building blocks: percentile, basic, studentized, and BCa intervals.
- Delete-one and delete-cluster jackknife.
- Reservoir and streaming sampling.
- Reproducible independent RNG stream splitting for parallel workers.
- Sampling plans remain separate from the statistic or estimator callback.

### 4. Derivatives and general optimization

Derivative contracts:
- Objective, gradient, Hessian, residual, and Jacobian callbacks.
- Analytic callbacks plus scale-aware finite differences.
- Hessian symmetry/curvature diagnostics and derivative checking.
- Sparse and matrix-free operator callbacks for large problems.

Unconstrained algorithms:
- BFGS, L-BFGS, Newton, damped Newton, trust-region Newton, nonlinear conjugate gradient, and Nelder-Mead.
- Safeguarded line searches, trust-region acceptance, damping/regularization, scaling, and explicit convergence reasons.
- Levenberg-Marquardt and Gauss-Newton primitives for nonlinear least squares.

Constrained continuous optimization:
- Bounds and projected-gradient methods.
- Linear and nonlinear equality/inequality constraints.
- Active-set, augmented-Lagrangian, and interior-point building blocks.
- KKT residuals, feasibility diagnostics, constraint scaling, and multiplier reporting.

### 5. Mathematical programming

Provide generic problem representations and solvers, without model-specific wrappers:
- Linear programming: primal/dual forms, feasibility, presolve, revised simplex, and primal-dual interior-point methods.
- Quadratic programming: convex QP first, then explicitly classified non-convex cases where supported.
- Nonlinear programming: smooth constrained problems with KKT-based termination and robust globalization.
- Mixed-integer and combinatorial optimization: branch-and-bound foundations, node queues, incumbent management, bounds, cuts/interfaces, and pluggable branching/relaxation callbacks.
- Network optimization: graph representations, shortest paths, max-flow/min-cut, and minimum-cost flow.
- Multiobjective optimization: Pareto dominance, nondominated archives, and scalarization helpers.
- Scaling, presolve, infeasibility/unboundedness detection, certificates where feasible, and deterministic diagnostics.

Solver APIs must distinguish optimal, feasible-but-not-proven-optimal, infeasible, unbounded, iteration/node limit, numerical failure, and interrupted termination. Never label a local optimum as a global optimum without proof.

### 6. Fuzzy computing

A reusable fuzzy layer should include:
- Membership functions: triangular, trapezoidal, Gaussian, generalized bell, sigmoid, and user-defined callbacks.
- Fuzzy sets with explicit universe/domain semantics and validated membership range.
- T-norms, t-conorms, complements, implication operators, and aggregation.
- Fuzzy relations, composition, and inference primitives.
- Linguistic variables and rule evaluation primitives without embedding a particular application DSL.
- Fuzzy numbers and interval/alpha-cut representations, arithmetic, comparison, and defuzzification methods such as centroid.
- Type-1 fuzzy foundations first; type-2 systems only after their representation and computational costs are justified.
- Fuzzy optimization: alpha-cut reformulation, possibilistic constraints, fuzzy objectives/constraints, and explicit uncertainty semantics.

Probability distributions and fuzzy membership functions are not interchangeable. Their semantics, normalization, and inference rules differ.

### 7. Statistical inference and econometric primitives

- Residual, score, moment, and information accumulation.
- Classical covariance and sandwich estimators.
- HC0-HC5, HAC/Newey-West, and cluster-robust covariance building blocks.
- Delta method, Wald, score, and likelihood-ratio test primitives.
- Confidence intervals, p-value/tail helpers, quantile helpers, and multiple-testing adjustments.
- Weighted, grouped, masked, and missing-data-aware reductions.
- Projection, orthogonalization, rank diagnostics, and stable weighting-matrix solves for IV/GMM.
- Robust covariance and moment accumulation reusable by simultaneous-equation and GMM implementations.

### 8. Time-series, panel, and multivariate primitives

- Lag, lead, first/seasonal differences, rolling and expanding reductions.
- Autocovariance, autocorrelation, cross-correlation, convolution, filtering, and spectral primitives.
- Block, banded, sparse, and state-space-ready linear algebra.
- Panel entity/time indexing contracts, within/between transformations, grouped operations, and cluster accumulation.
- Multivariate covariance transformations, whitening, projections, decompositions, and matrix-variate building blocks.
- Formula parsing and data-frame policy remain outside Core.

## Delivery sequence

### Phase 0 — Contract and inventory
Reconcile public/private headers, implementation, tests, CMake, dependency policy, and documentation. Establish independent reference sources, tolerances, status semantics, ownership rules, and benchmark baselines. Correct documentation claims that exceed the implemented public API.

### Phase 1 — Scalar numerical foundation
Implement special functions, stable log-domain helpers, online/mergeable accumulators, and robust root finding. Test domains, endpoints, tails, monotonicity, identities, and extreme magnitudes.

### Phase 2 — Probability and resampling
Implement the distribution contract and initial continuous/discrete families, sampling plans, stream splitting, bootstrap, and jackknife. Add statistical reference vectors, distributional property tests, and deterministic RNG tests.

### Phase 3 — Derivatives and optimization
Implement derivative callbacks/checkers, Newton and damped Newton, L-BFGS, trust-region and nonlinear least-squares primitives, bounds, and constrained continuous optimization. Add comprehensive convergence and failure diagnostics.

### Phase 4 — Mathematical programming
Implement LP representations, presolve, simplex and primal-dual methods, then QP and nonlinear constrained solvers. Add mixed-integer branch-and-bound and network optimization foundations. Cross-check small problems against independent reference solvers.

### Phase 5 — Fuzzy computing
Implement membership functions, fuzzy sets/numbers, operators, alpha-cuts, inference primitives, defuzzification, and fuzzy optimization formulations. Test algebraic laws and independent reference examples.

### Phase 6 — Inference, econometrics, time series, and panel
Implement robust covariance and hypothesis-testing primitives, IV/GMM moment and weighting infrastructure, time-series and panel transformations, and structured algebra. Add composition tests for representative calculations without model-specific Core APIs.

### Phase 7 — Hardening and release readiness
Validate static/shared package consumers, ABI, allocation failures, sanitizers, thread safety, supported-platform CI, numerical stress, deterministic reproducibility, benchmarks, dependency/license audit, documentation, examples, and binding integration.

## Completion criteria

A capability is complete only when:
1. Public and private contracts are documented and consistent.
2. Algorithms have independent numerical reference tests and adversarial edge-case coverage.
3. Allocation failures, overflow, invalid dimensions, non-finite inputs, and cleanup paths are tested.
4. Ownership, concurrency, determinism, and callback behavior are specified and tested where relevant.
5. Static/shared builds and package-consumer tests pass on the supported CI matrix.
6. Benchmarks show acceptable complexity and no unjustified dense allocations or explicit inversions.
7. Composition tests demonstrate that primitives can be combined into representative statistical, econometric, optimization, and fuzzy workflows.

No subsystem or the overall engine should be described as complete merely because its API has been declared or a small set of examples passes.
