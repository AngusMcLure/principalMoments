# principalMoments

[![R-CMD-check](https://github.com/AngusMcLure/principalMoments/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AngusMcLure/principalMoments/actions/workflows/R-CMD-check.yaml)

Compute the two principal representations of probability distributions on a
known compact interval, given raw moments `m1, ..., mv` (`m0 = 1` is implicit).
The package has no runtime dependencies beyond R.

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

This is a double-precision implementation intended for low to moderate moment
orders. Raw moment matrices become ill-conditioned at high orders and near
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
