# principalMoments

[![R-CMD-check](https://github.com/AngusMcLure/principalMoments/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AngusMcLure/principalMoments/actions/workflows/R-CMD-check.yaml)

Compute the two principal representations of probability distributions on a
known compact interval, given raw moments `m1, ..., mv` (`m0 = 1` is implicit).
The precision backend uses Rmpfr (and its GMP/MPFR libraries).

**Numerical status:** version 0.3.0 adds fail-safe rank handling and a
precision-refined MPFR backend. It recovers the 17-moment failure and has been
checked through 100 uniform moments with accurate inputs. This is not a universal
capacity guarantee. See [the progress report](ROBUSTNESS-PROGRESS.md) and the
[original 0.2.0 assessment](NUMERICAL-ROBUSTNESS.md) for evidence and limitations.

```r
# install.packages("pak")
pak::pak("AngusMcLure/principalMoments")

library(principalMoments)
reps <- principal_representations(c(1/2, 1/3, 1/4))
reps$lower  # two-node Gauss rule
reps$upper  # Simpson / Gauss-Lobatto rule
reps$diagnostics

evaluate_representations(reps, exp)
functional_bounds(reps, exp, derivative = "positive")
```

`lower` and `upper` name **representations**, not the ordering of arbitrary
functional values. For a continuous integrand whose `(v+1)`th derivative is
nonnegative on the interior, the lower representation minimises its expectation
and the upper maximises it (Markov-Krein). A nonpositive derivative reverses the
roles. `functional_bounds()` requires an explicit assertion of this condition;
it does not prove it or check it numerically. General oscillating integrands or
nonlinear distribution functionals need separate optimisation arguments.

For example, the clustered group-testing prevalence integrand satisfies the
positive derivative condition for every order when pool size `s > 1`:

```r
prevalence <- function(q, s) -expm1(log1p(-q) / s)
functional_bounds(reps, prevalence, s = 25, derivative = "positive")

# General distribution functionals can be evaluated without claiming bounds:
evaluate_representations(reps, function(nodes, weights) {
  sum(weights * nodes^2) - sum(weights * nodes)^2
}, type = "distribution")
```

Support defaults to `[0, 1]`; set `support = c(a, b)` for another finite interval.
The supplied moments must be those of the variable on that original interval.
Feasible boundary sequences return the same unique measure twice. Infeasible
moments are rejected; noisy moment estimates are not silently projected.

## Algorithms and numerical scope

| Number of moments | Lower | Upper |
|---|---|---|
| `v = 2n-1` | n-node Gauss | endpoints and n-1 interior nodes (Lobatto) |
| `v = 2n` | left endpoint and n interior nodes (Radau) | right endpoint and n interior nodes (Radau) |

Hausdorff localising matrices check feasibility. Interior quadrature uses
Cholesky orthonormalisation and a symmetric multiplication-matrix eigenproblem;
modified measures produce the endpoint-constrained rules. Boundary recovery
uses a null polynomial from a singular localising matrix. Every output must
reproduce **all** normalised supplied moments within `tol`.

The numeric-input fast path is intended for low to moderate moment orders;
use preserved-precision input for the MPFR backend. Raw moment matrices become ill-conditioned at high orders and near
boundary sequences. Diagnostic residuals are backward checks, not guarantees
of accurate nodes or bounds. Numerical rank depends on `tol`. Try nearby
tolerances, inspect warnings, and use higher precision or a better moment basis
for demanding problems. The package fails explicitly when it cannot recover
all constraints rather than dropping moments. See [DESIGN.md](DESIGN.md).

## Development

```r
install.packages(c("roxygen2", "testthat", "knitr", "rmarkdown", "rcmdcheck"))
roxygen2::roxygenise()
testthat::test_local()
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), error_on = "warning")
```

GitHub Actions builds the source package, builds the vignette, installs it,
and runs R CMD check (including testthat and examples) on Linux release/devel,
macOS release, and Windows release. The workflow uses standard GitHub-hosted
runners and read-only repository permissions. The repository is public.

## Background

The implementation follows classical moment-space and quadrature results:
Karlin & Studden (1966), *Tchebycheff Systems*, pp. 55-56;
Krein & Nudelman (1977), *The Markov Moment Problem and Extremal Problems*;
Gautschi (2004), *Orthogonal Polynomials*, Section 1.4.2.

[QuadMomPy](https://github.com/puetzmi/QuadMomPy) provides related moment-inversion
and endpoint quadrature algorithms. It was inspected as a reference; this is an
independent R implementation of classical mathematical algorithms, without
copied Python code or a Python dependency. The motivating application is
[clustered group testing](https://github.com/AngusMcLure/MX-PVS-sample-size/blob/clustered-group-testing-theory/clustered-group-testing-theory.Rmd).

## Selected powers and mixed pool sizes

`principal_representations()` interprets its input as consecutive raw moments.
For selected positive real orders on `[0,1]`, use the separate API:

```r
# One size-2 pool and one size-3 pool identify E(X^2), E(X^3), E(X^5),
# where X = 1-P. Here the example moments come from Uniform(0,1).
r <- generalised_principal_representations(1 / c(3, 4, 6), c(2, 3, 5))
r$lower
r$upper
b <- power_moment_bounds(r, exponent = 1)
c(minimum = 1 - unname(b[2]), maximum = 1 - unname(b[1]))
evaluate_representations(r, function(x) 1 - x)
```

The orientation for a target power depends on how many supplied exponents
exceed it. `power_moment_bounds()` handles that orientation automatically;
the ordinary derivative-sign `functional_bounds()` deliberately rejects these
objects. Arbitrary integrands can still be evaluated, without claiming bounds.

Equally spaced exponents reduce to ordinary moments by a power transformation.
Other exponents use MPFR Newton corrections with ordered positive support,
continuation, and precision-refined checks. Legacy double fits supply warm
starts only.
Both output measures must reproduce every constraint. Recognisable point-mass
and endpoint-mixture boundaries are recovered directly. Nonconvergence is
reported explicitly; it does not prove the input infeasible. This API covers
selected monomials, not arbitrary constraint functions or unbounded support.
It is intended for low to moderate dimensions; near boundaries, high or closely
spaced powers can be ill-conditioned. Check residuals and tolerance sensitivity.


## Preserving precision for longer sequences

```r
# Construct the moments accurately, rather than rounding first:
m <- 1 / Rmpfr::mpfr(2:51, 512)
r <- principal_representations(m)
r$diagnostics
functional_bounds(r, function(q) -expm1(log1p(-q) / 25),
                  derivative = "positive")

# Decimal strings are also accepted without an intermediate double conversion.
r <- generalised_principal_representations(
  c("0.9321971323", "0.8546368326", "0.8403270527"), c(10, 40, 50)
)
```

The short decimal example is interpreted as those supplied values, not the
unrounded analytical beta moments. Use enough digits for the desired problem.
Plain numeric input retains its 53-bit information content. Asking for higher
working precision cannot recreate digits lost in input rounding. Highly
ill-conditioned double input is rejected instead of silently losing moments.

`precision` selects an initial bit count, `max_precision` limits refinement,
`accuracy` controls node/weight and functional stability checks, and `tol` controls
absolute reproduction error. The MPFR path also checks relative moment errors.
The returned tables are convenient double-precision views; `high_precision`
retains MPFR nodes and weights, and the evaluator uses those values and compares
the functional across successive precision levels. User functions must support
MPFR arguments; do not coerce them to doubles inside an endpoint-sensitive
functional. Some R functions only accept ordinary numeric arguments.

Recognised boundary data may be labelled `boundary_at_input_precision`.
An inferred `numerical_boundary` is not enough to claim sharp bounds: the
bounds functions reject that status. General unrecognised singular cases can
remain unresolved. These checks are numerical, not interval-arithmetic or
input-uncertainty certificates. Use the manually triggered **Numerical stress
checks** workflow to run the longer 100-moment and high-power regressions.
