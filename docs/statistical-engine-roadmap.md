# Numelyx Statistical Computing Engine Roadmap

## Product objective

Numelyx is a general-purpose, language-independent numerical and statistical
computing engine. Its public API is a collection of stable, composable
primitives. Statistical models and estimators are composed above those
primitives by language bindings or higher-level packages.

The Core must not grow a separate model-fitting entry point for every estimator
(for example, dedicated `fit_tobit()`, `fit_panel()`, or `fit_ols()`
functions). Model-specific orchestration belongs in a higher layer unless a
small reusable numerical primitive is genuinely missing.

This roadmap is a capability target, not a claim that every item is currently
implemented.

## Current verified foundation

The current public C API provides:

- Opaque Matrix handles, structured storage, lazy views and arithmetic.
- LU, QR, reduced SVD, Cholesky, LDLT, symmetric eigen decomposition,
  triangular solves, ordinary and weighted least-squares primitives,
  pseudoinverse, numerical rank, norms, condition estimates, determinant and
  signed log-determinant, subspace bases and scalar reductions.
- Explicit-state PCG32 RNG, uniform and normal variates, cloning, and index
  sampling with and without replacement.
- Multivariate Gaussian log density and sampling.
- Unconstrained smooth-objective minimization using BFGS with Armijo
  backtracking, optional gradients, and finite-difference gradients.

These capabilities are a foundation. They do not yet amount to a complete
statistical distribution, resampling, inference, or optimization toolkit.

## Capability map by model family

| Model family | Required reusable primitives |
|---|---|
| Linear and weighted least squares | Matrix factorizations, rank diagnostics, triangular solves, stable residual and objective evaluation, covariance estimators |
| GLS / FGLS | Cholesky and triangular solves, structured covariance operations, log determinants, iterative reweighting primitives |
| Multivariate models | Multivariate distributions, covariance/correlation operations, symmetric eigensolvers, matrix-normal and related sampling foundations |
| Panel data | Grouped reductions, within/between transforms, entity/time indexing support, clustered covariance building blocks, weighted projections |
| Time series | Lag/lead and differencing views, autocovariance/cross-covariance, stable recursions, spectral primitives, state-space linear algebra and simulation primitives |
| Simultaneous equations / IV / GMM | QR/SVD-based projections, rank checks, moment accumulation, weighting-matrix solves, robust covariance building blocks |
| Probit and binary response likelihoods | Normal CDF, survival function, PDF, log-CDF/log-survival, quantile and stable likelihood/gradient primitives |
| Tobit and censored/truncated likelihoods | Distribution tails, log-CDF/log-survival, stable Mills-ratio-related scalar functions, derivative and optimizer primitives |
| ML and nonlinear statistical estimation | General objective/gradient/Hessian interfaces, automatic or numerical differentiation strategy, regularization primitives, scalable optimizers |
| Bootstrap and simulation inference | Resampling plans, sampling with/without replacement, stratified/cluster/block resampling, replicate execution support, reproducible independent streams |

The table describes the primitives needed to implement these methods. It does
not imply that each model is or should be implemented inside Core.

## Workstream A — probability distributions and special functions

### Distribution API contract

For each supported distribution, define only mathematically meaningful
operations and explicit parameter domains:

- PDF/PMF where representable.
- CDF and survival function.
- Log-PDF/log-PMF and log-CDF/log-survival for stable likelihood evaluation.
- Quantile/inverse-CDF where a robust implementation is available.
- Sampling using an explicit RNG handle.
- Parameter validation, endpoint behavior, non-finite behavior and status
  semantics.

Do not compute a likelihood as `log(pdf(x))` when a direct log-density formula
is available. Tail functions must avoid catastrophic cancellation; for example,
log-survival must not be implemented as `log(1 - cdf(x))` in extreme tails.

### Initial distribution families

1. Continuous univariate: Normal, Student t, chi-square, F, uniform,
   exponential, gamma, beta, logistic, Weibull, and log-normal.
2. Discrete univariate: Bernoulli, binomial, Poisson, geometric, and negative
   binomial.
3. Multivariate: multivariate normal first; then multivariate t, categorical
   and multinomial sampling where their numerical contracts are established.
4. Distribution utilities: random variates, truncation, and stable tail
   probability operations. Mixture distributions can follow once component
   composition and sampling contracts are clear.

### Special-function dependencies

Build and test the required scalar numerical foundation explicitly:

- log-gamma and gamma-related ratios.
- beta/log-beta.
- error function and complementary error function.
- regularized incomplete gamma and beta functions and their complements.
- Stable inverse functions needed by quantiles.
- Stable log-sum-exp, log-one-plus-exp, and tail/Mills-ratio helpers.

Select algorithms by domain and conditioning; include reference vectors,
monotonicity, normalization, symmetry, limiting cases and extreme-tail tests.
Use an established numerical dependency when it materially improves
correctness and maintenance, but keep the public API independent of a vendor's
types and ABI.

## Workstream B — sampling and resampling

Current index sampling is a useful base, but a statistical engine needs more
than a single index-array function.

Required capabilities:

- Sampling with and without replacement, with explicit reproducibility
  contracts.
- Weighted categorical sampling with validated nonnegative weights and
  documented zero-total behavior.
- Stratified and grouped sampling.
- Cluster resampling.
- Block and moving-block resampling for dependent/time-series observations.
- Permutation and randomization plans.
- Bootstrap replicate orchestration primitives, including percentile,
  basic, studentized and BCa interval building blocks.
- Jackknife delete-one and grouped/delete-cluster plans.
- Reservoir sampling or streaming-compatible sampling where appropriate.
- Explicit RNG stream management/splitting so parallel workers can reproduce
  independent, non-overlapping streams.

Resampling must separate the sampling plan from the statistic callback. Core
should provide deterministic plans and safe iteration primitives rather than
embedding model-specific estimators. Document whether each algorithm is
ordered, whether duplicates are possible, memory complexity, and exact RNG
consumption guarantees.

## Workstream C — optimization and derivatives

BFGS is the current optimizer. It is not sufficient as the only optimization
method for statistical estimation.

### Optimizer methods to add

Prioritize reusable unconstrained methods:

1. Newton-Raphson with Hessian solve and globalization.
2. Damped Newton / trust-region Newton for indefinite or poorly conditioned
   Hessians.
3. L-BFGS for large parameter vectors.
4. Nonlinear conjugate gradient.
5. Derivative-free Nelder-Mead as a fallback for low-dimensional objectives.

Then evaluate constrained methods:

- Box constraints and projected-gradient methods.
- General constraints only after a clear constraint/Jacobian contract exists.
- Levenberg-Marquardt for nonlinear least-squares residual problems.

Every optimizer must report convergence reason, iteration/evaluation counts,
gradient norm, step norm, objective value, and relevant diagnostics. A failed
line search or iteration limit must not be reported as convergence.

### Derivative infrastructure

Define reusable callback contracts for objective, gradient, Hessian and
residual/Jacobian functions. Provide scale-aware finite-difference fallbacks
with explicit step policies and tests. Consider automatic differentiation only
as a separately evaluated architectural decision; do not pretend finite
differences are equivalent to exact derivatives.

Use factorizations and linear solves, not explicit Hessian inversion. Add
curvature checks, regularization/damping, line-search safeguards and
ill-conditioning diagnostics.

## Workstream D — inference and covariance primitives

Model estimation needs reliable inference building blocks beyond coefficient
optimization:

- Residual and score accumulation.
- Observed and expected information matrices.
- Classical covariance from stable factorization/solve paths.
- Sandwich/robust covariance.
- Heteroskedasticity-consistent covariance estimators.
- HAC/Newey-West covariance with documented kernel and bandwidth choices.
- Cluster-robust covariance, including multiway support if its contract is
  justified.
- Delta-method transformations and covariance propagation.
- Wald, score and likelihood-ratio statistic primitives.
- Confidence-interval and p-value helpers built on distribution tails.
- Missing-data masks and weighted/grouped reductions where needed.

These should be low-level operations over matrices, vectors and callback
results, not model-specific fit APIs.

## Workstream E — structured data and time-series primitives

Some model families require reusable data transformations not supplied by
dense matrix algebra alone:

- Lag, lead, difference and seasonal-difference views or transforms.
- Stable rolling and grouped reductions.
- Autocovariance, autocorrelation and cross-correlation.
- Convolution/filtering and spectral density building blocks.
- Block/banded structured matrices and sparse operations suitable for
  state-space systems.
- Group/entity/time indexing contracts for panel transformations.
- Missing-observation and unbalanced-panel handling at the data-adapter layer.

Keep data-frame concerns, formula parsing and model syntax outside Core. Core
should accept validated numerical buffers, index/group descriptors and
well-defined layouts.

## Workstream F — external numerical dependencies

Audit candidate dependencies against actual algorithm needs before adoption.
For each dependency, record:

- Which algorithms it supplies and why the in-house implementation is not
  preferable.
- License compatibility, supported platforms, build modes and ABI exposure.
- Thread-safety, allocator ownership, error handling and reproducibility.
- Whether it is required or optional at build/runtime.
- Version pinning and security/update policy.

Do not add a dependency merely because it is popular. Do not expose third-party
types in the public ABI. Prefer a narrow internal adapter and retain a
consistent public Numelyx status/ownership contract.

## Delivery phases

### Phase 1 — statistical foundation

- Complete the distribution/special-function design and numerical contracts.
- Add core continuous and discrete distributions with log-tail operations.
- Expand sampling to weighted, stratified, cluster and block plans.
- Add bootstrap and jackknife plan primitives.
- Add Newton-Raphson, damped/trust-region Newton and L-BFGS.
- Add gradient/Hessian and convergence diagnostics contracts.

### Phase 2 — inference foundation

- Add robust, HC, HAC and cluster covariance primitives.
- Add score/information, Wald/score/LR and delta-method utilities.
- Add probability quantiles and confidence interval helpers.
- Add numerical reference and property-based tests for all new algorithms.

### Phase 3 — dependent and structured observations

- Add time-series lag/difference/autocovariance/spectral primitives.
- Add panel grouping and within/between transformations.
- Add block/bootstrap and cluster resampling integration tests.
- Add state-space-ready structured linear algebra as justified by workloads.

### Phase 4 — production hardening

- Benchmark representative workloads across matrix sizes and distributions.
- Test deterministic behavior across supported platforms where promised.
- Run sanitizer, race, allocation-failure and adversarial numerical tests.
- Validate public API ownership, ABI, package consumption and language binding
  integration.
- Publish explicit numerical accuracy and compatibility guarantees.

## Acceptance criteria

A capability is not complete merely because a function exists. Each addition
must include:

1. A written mathematical and API contract.
2. Stable algorithm selection and documented numerical limitations.
3. Input-domain, overflow, non-finite and failure-path behavior.
4. Native tests with trusted reference values and adversarial edge cases.
5. Property tests (for example CDF monotonicity, PDF non-negativity,
   normalization, inverse-CDF consistency, residual checks and sampling
   frequencies where appropriate).
6. Reproducibility tests when RNG is involved.
7. Allocation-failure and sanitizer coverage for ownership-bearing paths.
8. Public-header, package-consumer and cross-platform CI coverage.
9. Documentation and benchmark evidence for performance-sensitive paths.

## Explicit non-goals for Core

- Formula parsing, data frames, database access and file-format ingestion.
- PHP/Zend or any other language runtime dependency.
- Model-specific one-call fitting APIs for every statistical estimator.
- Cryptographic randomness.
- Silent changes to distribution algorithms or RNG sequences without a
  documented compatibility/version decision.
