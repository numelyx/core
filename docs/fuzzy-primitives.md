# Fuzzy Membership and Algebraic Primitives

Status: initial internal scalar fuzzy primitives added; validation is ongoing.

## Membership functions

The current layer provides triangular, trapezoidal, Gaussian, generalized-bell, and sigmoid membership functions. Parameter domains are checked and valid outputs are in [0,1]. Shoulder forms are supported by allowing adjacent equal breakpoints in triangular and trapezoidal functions while requiring a nonzero overall support interval.

## Fuzzy operators

The layer provides complement, minimum/product/Lukasiewicz t-norms, and maximum/probabilistic-sum/bounded-sum s-norms. All operator inputs must be finite memberships in [0,1].

## Rule aggregation

Mamdani aggregation accepts rule firing strengths and a row-major rule-by-sample consequent-membership matrix. It supports minimum/product implication and maximum/probabilistic-sum/bounded-sum aggregation. Sugeno aggregation computes the firing-strength-weighted average of crisp rule outputs.

## Defuzzification

Sampled membership arrays support centroid and mean-of-maxima defuzzification. Coordinates are scaled before weighted accumulation to reduce overflow risk; an all-zero membership set is rejected as having no defined crisp representative.

## Alpha cuts

Triangular and trapezoidal fuzzy-number alpha cuts are computed directly for alpha in [0,1]. The implementation allocates no memory and leaves outputs unchanged on validation failure.

## Deferred

General membership inversion, fuzzy relations, and fuzzy optimization remain future work. Sampled-set centroid and mean-of-maxima defuzzification are available, along with Mamdani implication/aggregation and weighted-average Sugeno rule aggregation. Interval addition, subtraction, multiplication, and division are now available as the arithmetic kernel for alpha-cut extension operations.

## Validation

Native tests cover membership reference values, operator identities, alpha-cut endpoints/interiors, interval arithmetic across sign combinations, zero-crossing division rejection, overflow handling, centroid/mean-of-maxima defuzzification, Mamdani/Sugeno aggregation, invalid parameters, non-finite inputs, and output preservation.
