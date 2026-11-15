# Design and algorithm contract

## Scope and API

`principal_representations(moments, support = c(0, 1), tol = 1e-12)` accepts
probability raw moments m1 through mv. It returns sorted node/weight tables for
both principal representations, metadata, and numerical diagnostics.
`evaluate_representations()` accepts an integrand or a scalar distribution
functional. `functional_bounds()` adds the caller-asserted Markov-Krein condition.

A function nonlinear in x still has an expectation linear in the probability
measure. A functional nonlinear in the measure is different: principal
representations need not optimise it. Neither evaluation nor sorting two
computed values solves general moment-constrained optimisation.

Support must be a known finite interval with positive width. The first release
supports probability measures only, not arbitrary masses or unbounded supports.
Transform moments by the binomial formula for Y=(X-a)/(b-a), then work on [0,1].
Cancellation in this transformation is another source of error. Nodes are
mapped back only after normalised moment checks.

## Feasibility and boundary recovery

Write H_r^(j) = [m_(k+l+j)] for k,l=0,...,r-1, with m0=1.
For v=2n, feasibility requires H_(n+1)^(0) and H_n^(1)-H_n^(2)
to be positive semidefinite. For v=2n+1, it requires H_(n+1)^(1)
and H_(n+1)^(0)-H_(n+1)^(1) to be positive semidefinite.
These are necessary and sufficient truncated Hausdorff conditions in exact
arithmetic. Eigenvalues below -tol reject feasibility. Values at or below tol
trigger boundary recovery; eigenvalue ratios diagnose poor conditioning.

If a localising matrix is singular, its null vector defines a polynomial p.
The integral of p(x)^2 times its nonnegative multiplier is zero, so the unique
boundary distribution is supported at real roots of p plus zeros of the
multiplier (none, 0, 1, or both endpoints). Solve for weights against all
available moments using QR; verify all constraints. Multiple null vectors may
exist, but one annihilating polynomial and its multiplier suffice to obtain a
finite candidate support. Extra candidate nodes can have zero weight. Boundary weights at or below tol
in magnitude are set to zero, then all constraints are checked again.

Boundary status is tolerance dependent. If recovering the measure fails, do
not label an arbitrary lower-rank approximation a solution. Return an error.

## Interior principal representations

For an r-node Gaussian rule under moments t0,...,t_(2r-1), form
H=[t_(i+j)] and B=[t_(i+j+1)]. If H=R'R, the multiplication matrix is
J=R^(-T) B R^(-1). Its eigenvalues are nodes and t0 times squares of the
first components of its normalised eigenvectors are weights. This is the
symmetric orthogonal-polynomial quadrature construction (Golub-Welsch in an
orthonormal polynomial basis). Dense cost is cubic in r.

Use the following modified measures to obtain the interior nodes and weights:

| v | representation | multiplier | r |
|---|---|---|---|
| 2n-1 | lower | 1 | n |
| 2n-1 | upper | x(1-x) | n-1 |
| 2n | lower | x | n |
| 2n | upper | 1-x | n |

Divide modified quadrature weights by the multiplier at each interior node.
For Radau recover endpoint mass from total probability. For Lobatto recover
right mass from the first moment, then left mass from total probability.
The n=1 Lobatto rule consists only of endpoints. No interior Vandermonde solve
or monomial polynomial-root computation is used.

Every result must have nodes in support, nonnegative weights, mass one and
reproduce all normalised moments within tol. Only negative roundoff up to tol
is clipped; normalisation and clipping are followed by the full residual check.

## Validation

Tests use independent closed forms for v=1,2,3, polynomial exactness, uniform
and beta moments, seeded discrete measures, analytic group-testing bounds,
positive and negative derivative cases, arbitrary compact supports, finite
boundary distributions and invalid inputs. Cross-platform CI runs the installed
package checks and vignette. Analytic fixtures provide an independent oracle
without requiring Python. No claim of formal numerical certification is made.

## Next numerical work

1. Modified-moment / orthogonal-basis input to reduce raw Hankel conditioning.
2. Canonical-moment or modified-Chebyshev recurrence backend (compare QuadMomPy's
   Wheeler and product-difference algorithms) with explicit rank diagnostics.
3. Optional arbitrary precision and stress benchmarks near Hausdorff boundaries.
4. Optional projection of estimated moments as a separate, explicitly requested
   statistical operation, never hidden inside the exact-moment API.
5. General integrand optimisation via polynomial envelopes as a separate API.

These extensions are not implemented in version 0.1.0.
