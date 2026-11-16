# Robustness work in version 0.3.0

This is an implementation follow-up to the preserved
[version-0.2.0 assessment](NUMERICAL-ROBUSTNESS.md). It is a substantial first
precision backend, not a claim that arbitrary moment problems are now uniformly
well-conditioned or certified.

## Implemented

- Decimal-string and Rmpfr inputs preserve digits from the beginning. Numeric
  inputs retain their 53-bit information content. Rmpfr is now a runtime dependency;
  Python is not required to use or test the R package.
- Long ill-conditioned double-input sequences fail explicitly instead of being
  automatically declared exact boundary distributions. Approximate numerical
  boundary rules have a distinct status, and the sharp-bound APIs reject
  unresolved boundary membership.
- The precision consecutive backend constructs orthogonal-polynomial Stieltjes
  recurrences in MPFR arithmetic. Symmetric-Jacobi eigenvalues are only starting
  guesses: the actual nodes are refined in MPFR and weights are computed from
  Christoffel sums. Modified measures give the endpoint-constrained rules.
- Generalized rules use analytic MPFR Newton corrections with strictly ordered
  interior nodes and positive weights. Existing double fits are warm starts
  only. Otherwise continuation first changes the exponents under the uniform
  measure and then changes the moments through convex combinations approaching
  the supplied vector. Very small positive parameters use logarithmic Newton
  coordinates with a step limit. There is no clipping in the accepted precision
  solve. Stopping uses the requested accuracy, rather than demanding gratuitous
  digits simply because working precision is high.
- Success requires a second working-precision calculation, stability of nodes
  and weights, absolute and relative reproduction checks, and (for numeric
  inputs on this path) a local input-sensitivity check. A bounded Newton work
  budget prevents unlimited continuation retries at any one precision.
- Returned `high_precision` and `precision_reference` fields preserve both sets
  of MPFR rules. Evaluation uses these instead of the rounded display tables
  and checks the requested functional separately across precision levels.
- Arithmetic-order transformations promote to MPFR automatically when a double
  node transformation would no longer reproduce the input moments. This keeps
  a small nonzero prevalence even when the display node rounds to one.

## Evidence

The committed `tests/testthat/fixtures/precision-*.csv` files freeze independent
references from the original 100-digit calculations, plus Gauss-Jacobi fixtures
for 50 and 100 uniform moments. The latter were generated using SciPy's
classical Jacobi rules, independently of the R recurrence implementation.
The fixtures are plain CSV, so CI does not depend on Python or SciPy.

| Case | Version 0.2.0 | Version 0.3.0 with preserved input precision |
|---|---|---|
| Uniform Q, 17 moments, pool size 25 | Incorrect collapsed interval | Recovers the reference noncollapsed interval |
| Beta-generated prevalence 1%, 17 moments, size 25 | Incorrect interval excluding truth | Reference nodes, weights and bounds recovered |
| Beta-generated prevalence 0.1%, 16 moments, size 50 | Incorrect interval excluding truth | Reference nodes, weights and bounds recovered |
| Uniform Q, 50 moments | Not a reliable principal representation | Nodes agree with independent Jacobi rules to about 1.1e-16; bound error about 4.3e-15 |
| Uniform Q, 100 moments | Not reliable | Nodes agree to about 1.1e-16; bound error about 9.5e-16 |
| Known lower rule, orders 10,40,50 | Did not converge | Recovered from the high-precision fixture |
| Known lower rule, orders 25,26,51 | Prevalence error 0.00255 | Reference nodes/weights/value recovered within 1e-8 |
| Known lower rule, orders 2,101,103 | Prevalence error 0.0218 | Recovered from the high-precision fixture |
| Known lower rule, orders 100,101,201 | Severe small-moment acceptance error | Recovered from preserved high-precision input |

Initial recurrence benchmarks on the assessment machine took about 46 seconds
for 50 moments and 186 seconds for 100 moments, at working precisions of 528
and 928 bits with 512-bit input. The final initializer honours supplied precision
and reserves room for a second calculation: that input now starts at 512 bits
and refines at 1024 bits. Timings are illustrative, not performance promises. High precision is materially more expensive than the
old double path. The difficult high-power cases use an explicit `maxit=300` budget.
The 100-moment and high-power cases are in the separate manual
stress workflow; the normal matrix includes the 50-moment and former-failure
regressions.

For the 17-moment uniform case the accurate interval is approximately
[3.82446949%, 4.71995651%]. Passing rounded double moments without adequate
information is deliberately allowed to return an unresolved-input error;
raising working precision is not presented as recovery of discarded digits.

Some previous double-input boundary tests now distinguish recovery of a
numerical boundary representation from an exact sharp-bound claim. They still
check the recovered nodes, weights and evaluations against their formulas, and
also require the sharp-bound API to reject unresolved boundary status. Known
point masses and endpoint mixtures have analytic paths at input precision.

## Use

```r
# Preserve precision during moment generation, not only in the solver:
m <- 1 / Rmpfr::mpfr(2:51, 512)
r <- principal_representations(m)
r$diagnostics
functional_bounds(r, function(q) -expm1(log1p(-q)/25),
                  derivative = "positive")

# Or supply sufficiently accurate decimal strings for selected moment orders:
r <- generalised_principal_representations(moment_strings, c(10, 40, 50))
b <- power_moment_bounds(r, 1)
c(minimum = 1 - unname(b[2]), maximum = 1 - unname(b[1]))
```

`tol` controls absolute reproduction error. `accuracy` controls precision-level
stability of nodes/weights and functional evaluations; the precision backend
also requires relative reproduction errors below `accuracy/100`. `precision`
sets the initial bit count and `max_precision` caps refinement. `maxit` controls
Newton work in generalized solves. The total per-precision Newton budget is
four times `maxit`; the bound includes correction retries and both sides.

`high_precision` values should be used through the evaluator for functions that
are sensitive near endpoints. The public data frames are double-precision views.
User functions must support MPFR arguments without coercing away their precision.

## Limits and remaining work

- This is numerical precision refinement, not interval-arithmetic certification.
  Local inverse-Jacobian input sensitivity is a rejection diagnostic, not a
  rigorous global uncertainty enclosure. No dual-envelope certificates are claimed.
- Raw-moment conditioning is not magically removed: recurrence inner products
  still use polynomial coefficients and precise moments. Direct modified-moment
  and location-count-probability input remain future work.
- Genuine singular boundary problems beyond recognized point masses/endpoint
  mixtures may remain unresolved. The precision path currently requires [0,1];
  rescale before losing precision for other supports.
- Larger generalized problems can still be expensive or unresolved. There is
  **no advertised generalized dimension limit** comparable to the tested
  100-moment uniform consecutive example. Arbitrary-precision dense elimination
  in R is a performance bottleneck, and more specialized generalized quadrature
  algorithms remain worthwhile. A 23-selected-moment exploratory fixture
  (five size-2 and five size-3 pools, 25 individuals) did not complete within
  685 elapsed seconds and was stopped. It is not counted as a successful
  large-dimension result. After the accuracy-based stopping change, a repeat
  also remained unresolved within a 240-second external timeout.
- The input-collision example in the original report is still fundamental:
  if two distributions have the same rounded double moments, their original
  distinct prevalences cannot be reconstructed from those rounded inputs.
- Modified/count-basis input, rigorous uncertainty sets, certified dual bounds,
  and faster high-dimensional generalized linear algebra remain later stages.

Run the ordinary checks with `rcmdcheck::rcmdcheck()`. For the extended cases,
set `PRINCIPAL_MOMENTS_STRESS=true` when running tests, or trigger the
**Numerical stress checks** GitHub Actions workflow manually.


The uniform timing sweep and the difficult mixed fixture can be reproduced
from the repository root with
`Rscript benchmarks/robustness/precision-check.R` or, on Linux,
`timeout 240 Rscript benchmarks/robustness/precision-check.R --mixed23`.
The timeout is an external benchmark limit, not a convergence certificate.
