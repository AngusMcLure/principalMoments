# External GeneralizedGauss.jl experiments

These optional experiments support [the algorithm review](../../GENERALIZED-QUADRATURE-REVIEW.md).
They are not R package dependencies, R tests, or a supported solver backend.
`compare.jl` calls the Julia package's public API without modifying its source.
The additional flags only change the basis, supplied moments and solver options.

## Setup

Use Julia 1.11.7 and check out the reviewed fork:

```sh
git clone https://github.com/alexbercik/GeneralizedGauss.jl /tmp/GeneralizedGauss-review
git -C /tmp/GeneralizedGauss-review checkout cad739ca6f621cfae0b290b56da8e7b28ad2ab69
julia --project=/tmp/GeneralizedGauss-review -e 'using Pkg; Pkg.instantiate()'
julia --project=/tmp/GeneralizedGauss-review -e 'include("/tmp/GeneralizedGauss-review/test/runtests.jl")'
```

The full upstream suite passed 155 assertions in this review. Installation and
first-use compilation are substantial; do them before timing experiments.
Dependency versions are in `environment.txt`.

## Cases

All cases include normalization and pass explicit moments. No case integrates
against an assumed target density inside `compute_gauss_rule`.

- `radau`: orders 2,5, uniform-distribution moments. The lower rule has node
  z=(1/2)^(1/3), interior weight (1/3)/z^2 and the remaining weight at zero.
  It also runs as an initial compilation/warm-up case.
- `theory50`: orders 10,40,50, X~Beta(8.91,0.09), the mixed-size example in
  the theory document. Moments are products of exact decimal/rational factors
  computed at working precision.
- `mixed50`: orders 10,40,50, prescribed lower nodes (1/3,2/3) and the same
  weights. This matches the mathematical fixture behind regression case 198.
- `high201`: orders 100,101,201 with the same prescribed lower rule, matching
  the mathematical fixture behind case 216.
- `mixed23`: orders {2,3,...,23,25}, prescribed lower nodes i/13 and weights
  i/78, i=1,...,12. This is the fixture from `precision-check.R --mixed23`.

## Running

Run from the R repository root. Configurations have three colon-separated
fields: case name, `qr` or `raw`, and `strict` or `hybrid` intermediate solves.
Every configuration attempts both principal representations.

```sh
GG_BENCH_TOLERANCE=1e-100 julia --project=/tmp/GeneralizedGauss-review \
  benchmarks/generalized-gauss/compare.jl theory50:raw:strict mixed50:qr:strict high201:qr:strict

GG_BENCH_TOLERANCE=1e-100 GG_BENCH_ENDPOINT=right \
  julia --project=/tmp/GeneralizedGauss-review \
  benchmarks/generalized-gauss/compare.jl mixed50:qr:strict

GG_BENCH_TOLERANCE=1e-100 GG_BENCH_SCALE=first GG_BENCH_ENDPOINT=right \
  julia --project=/tmp/GeneralizedGauss-review \
  benchmarks/generalized-gauss/compare.jl theory50:qr:strict mixed50:qr:strict high201:qr:strict

GG_BENCH_BITS=768 GG_BENCH_TOLERANCE=1e-100 GG_BENCH_SCALE=first GG_BENCH_ENDPOINT=right \
  julia --project=/tmp/GeneralizedGauss-review \
  benchmarks/generalized-gauss/compare.jl theory50:qr:strict
```

Defaults are 512 bits, tolerance 1e-40, unscaled powers and the fork's natural
endpoint pairing. `GG_BENCH_SCALE=first` uses y=x^lambda_1 and maps the returned
nodes back to x before checking original moments. `GG_BENCH_ENDPOINT=right`
sets `add_endpoint=:right` for even basis length; odd lengths retain the required
natural pairing. QR always transforms supplied moments by the returned matrix.

For expensive trials, use a process limit. On Linux, for example:

```sh
GG_BENCH_TOLERANCE=1e-40 GG_BENCH_ENDPOINT=right \
  timeout --kill-after=5s 120s julia --project=/tmp/GeneralizedGauss-review \
  benchmarks/generalized-gauss/compare.jl mixed23:raw:hybrid mixed23:qr:hybrid
```

Repeat with `strict`, `GG_BENCH_SCALE=first`, and/or `qr` to investigate other
paths. `--kill-after` matters: one Julia process did not exit promptly on TERM.
A whole-process limit includes loading, compilation and earlier configurations;
an unfinished configuration was not necessarily given the full time limit.

## Interpreting results

`RESULT` records maximum absolute and relative residuals of the *original*
moments, including normalization; known lower node/weight error where available;
admissibility, prevalence, and node count. The final script also independently
checks the expected endpoint pattern and prints `screened` for admissible rules
with absolute residual <=1e-9, relative residual <=1e-10 and known-rule error
<=1e-8 where available. A screen is not a rigorous accuracy certificate.

`RULE` preserves decimal nodes and weights for precision-repeat comparisons.
`FAIL` records a caught solver error. A `START` with no subsequent result means
the external process was stopped or failed before completing that side.
`results.json` preserves the actual exploratory outputs; `upstream-tests.txt`
preserves the upstream summaries. Earlier logs precede the final reporting-only
addition of the endpoint-pattern/screen checks and the integer-power evaluation
optimization, so not every result has a `screened` field.

Several runs overlapped and timings include some JIT compilation. They should
not be used for speed ratios. No tested 23-moment configuration produced a
validated pair. These observations do not show that every parameter setting
fails, nor that the canonical algorithm cannot solve this fixture after changes.
