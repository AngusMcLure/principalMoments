# Review of GeneralizedGauss.jl for mixed pool sizes

Reviewed 16 November 2026 against principalMoments `e534928` (version 0.3.0).
This is an algorithm review and an external benchmark, not an implementation
change or a new supported-dimension claim.

## Recommendation

Prototype a prescribed-moment, prefix/canonical continuation backend in R,
using the extended GeneralizedGauss algorithm as the main reference. These
methods apply directly, but the fork is not a reliable replacement on all our
fixtures and the experiments do not establish a larger supported dimension.
Combine its continuation and triangular QR ideas with the package's existing MPFR,
boundary, admissibility and accuracy checks. Retain the present solver as a
fallback during validation. Julia is useful as an independent benchmark, but
need not become a dependency of the R package.

## Material reviewed

- Huybrechs, *On the computation of Gaussian quadrature rules for Chebyshev
  sets of linearly independent functions*, [arXiv:1710.11244](https://arxiv.org/abs/1710.11244),
  including the conditioning discussion in section 5.
- Bercik, Patrascu and Zingg, *Construction and Optimization of
  Summation-By-Parts Operators for General Function Spaces Using an Improved
  Generalized Gaussian Quadrature Algorithm*,
  [arXiv:2607.08934v1](https://arxiv.org/abs/2607.08934v1), including appendices A and B.
- [Original repository at a120f4e](https://github.com/daanhb/GeneralizedGauss.jl/tree/a120f4ee890b1dec29cb94ec9e6db78924924d13):
  every Julia source and test file, README, project metadata and license.
- [Extension at cad739c](https://github.com/alexbercik/GeneralizedGauss.jl/tree/cad739ca6f621cfae0b290b56da8e7b28ad2ab69):
  every Julia source, test and example file, README, project metadata and license.
  Dependencies' implementations were outside the review scope.
- The mixed-size identification section of the current
  [theory document](https://github.com/AngusMcLure/MX-PVS-sample-size/blob/clustered-group-testing-theory/clustered-group-testing-theory.Rmd),
  and our current precision solver and robustness reports.

Both Julia repositories use MIT licenses. Any source translated into the R
package must retain the relevant notices and attribution.

## Why the methods apply

For mixed pool sizes our basis is

\[
\phi(x)=(1,x^{\lambda_1},\ldots,x^{\lambda_d}),\qquad
c=(1,m_{\lambda_1},\ldots,m_{\lambda_d}).
\]

The sorted powers, including every prefix, form a complete Chebyshev system.
Quadrature equations are exactly our prescribed moment equations:
\(\sum_i w_i\phi_j(x_i)=c_j\). The algorithms can take these moments
directly: there is no need to recover a density, assume beta mixing, or
numerically integrate against the unknown distribution. The *current original
code* also accepts explicit moments; this is not solely a feature of the fork.

Let K=d+1 count the normalization function as well as nonconstant moments.

| Moment count | Lower principal representation | Upper principal representation |
|---|---|---|
| d=2n-1, K=2n | n interior nodes (Gauss) | Both endpoints and n-1 interior nodes (Lobatto) |
| d=2n, K=2n+1 | Left endpoint and n interior nodes (left Radau) | Right endpoint and n interior nodes (right Radau) |

These are exactly the patterns in the theory document. Computing them does
not change the parity-dependent orientation of prevalence bounds. The target
remains \(p=1-\sum_iw_ix_i\), evaluated on both rules using the existing
`power_moment_bounds` orientation logic.

## What changes algorithmically

Our `mp_selected` currently deforms the entire exponent vector from consecutive
powers, then deforms the entire moment vector from a uniform reference to the
target. Every correction works with the final system dimension. It already
uses analytic derivatives, exponent scaling and arbitrary precision, but
initialization and dense R-level MPFR elimination remain expensive.

Huybrechs instead adds constraints successively. Between principal rules, a
canonical rule fixes one node at a moving parameter xi. The next moment has
a monotone residual along the appropriate canonical family. Interlacing gives
the admissible interval. A sign-changing bracket supplies the next principal
rule's starting point. The scalar outer structure is valuable, but each
canonical evaluation still requires a nonlinear solve; this is not a reduction
of the entire calculation to cheap scalar arithmetic.

The extension implements adaptive continuation steps, saved brackets and
bounded secant refinement when a principal solve fails. It supports final
Lobatto and both Radau rules as well as Gauss. The original public driver is
restricted to an even number of functions and a final lower/open rule, even
though its internal representation types cover more cases. Its main loop also
has commented-out convergence assertions. The fork is the more useful
implementation reference for our two-sided bounds.

The fork's PDE/SBP construction and optimization results do not establish a
speedup or a maximum dimension for our moment inversion problem.

## Basis conditioning

The extension forms a weighted evaluation matrix
\(A_{qj}=\sqrt{\omega_q}\phi_j(t_q)\), then performs unpivoted
Householder QR, A=QR. Its lower-triangular transformation B=R^{-T} gives
\(\psi=B\phi\). We must also transform the supplied moments: \(c_\psi=Bc\).
Otherwise we solve a different moment problem.

Triangular structure preserves the successive prefix spans used in
continuation. QR avoids the extra condition-number squaring caused by forming
the Gram matrix A'A first. A known, positive *reference* measure, such as uniform
measure, is sufficient for constructing this transformation; it need not be
the unknown prevalence distribution. Quadrature resolution must be checked
for high powers rather than chosen solely from the number of functions.

This preconditioning does not make the inverse problem well-conditioned,
recover digits lost from double inputs, or remove the need for MPFR. Applying
B to rounded moments can amplify their uncertainty. Working precision must
cover both the transformation and the solve. A future implementation should
batch basis/derivative evaluation instead of repeatedly evaluating every
transformed function as an independent closure.

## Safeguards to keep

1. **Initialization and endpoint choice.** For our powers the first one-node
   rule is available analytically. In the unscaled theory example its node
   is approximately 0.99300349. The fork's scalar Newton step leaves the
   support, then its fixed interior-only bracket search reaches only
   0.98979592 from its initial midpoint and misses this root. Our existing
   substitution y=x^lambda_1 makes the first basis function linear and
   avoids this problem. Also avoid releasing a node at zero with every
   basis derivative zero, which occurs for unscaled powers above one;
   the right-anchored path works better on the tested unscaled examples.
   Rescaling introduces fractional powers, so corrections must also prevent
   invalid trial evaluations at negative nodes.
2. **Boundary handling.** Canonical continuation is designed for interior
   moment vectors. Point masses, endpoint mixtures, vanishing weights and
   coalescing nodes need separate handling. Do not replace existing analytic
   paths with an interior-only solver.
3. **Admissibility.** Retain positive-weight, ordered-node, support checks in
   every accepted correction. The fork's analytic NLsolve path is not a
   constrained positivity/order solver.
4. **Accuracy independent of solver flags.** The fork defaults to relaxed
   intermediate tolerances and permits final residual slack through
   `lost_digits`. For QR bases it raises the effective slack to at least its
   estimated digit loss, even if the caller requests `lost_digits=0`.
   Recheck all original moments with absolute and relative tolerances, repeat
   at higher precision, and check nodes, weights and the requested functional.
5. **Basis diagnostics are not proofs.** Random determinant sampling and
   numerically estimated Wronskian/Lipschitz bounds are useful diagnostics.
   For our power basis the analytic Chebyshev property is already known.
   A Wronskian can vanish at zero for valid selected powers; numerical
   diagnostics must not turn that into a false rejection.
6. **Use analytic derivatives.** Finite differences and derivative-free
   OrthoMADS add little for known powers. The MADS implementation searches
   Float64 node coordinates even when evaluations use BigFloat, so it is not
   the backend for high-precision node recovery.
7. **Overall continuation budget.** The fork limits retries for each proposed
   canonical step, but its sweep can accept extremely small steps, retain that
   step size and keep iterating. It needs a total iteration/time budget and
   progress checks, not only a per-step retry budget. Several experiments had
   to be stopped externally.

## Proposed implementation sequence

1. Add an internal canonical-continuation backend behind the current public
   API. Reuse input validation, arithmetic-order shortcuts, MPFR arithmetic,
   boundary handling and result metadata. Use analytic first derivatives.
2. Build a reusable triangular QR transformation at working precision, pass
   transformed *supplied* moments, and validate against the original powers.
3. Use adaptive canonical predictor/corrector steps with retained monotone
   brackets, safeguarded scalar refinement, and explicit work budgets.
   A tangent predictor from the canonical Jacobian is a possible later
   improvement over copying the previous seed.
4. Share the prefix ladder when computing both final representations.
   For even K, a Radau rule for K-1 functions can seed the final Gauss and
   Lobatto branches. For odd K, the Gauss rule for K-1 functions can seed
   left and right Radau branches. Avoid recomputing two complete ladders.
5. Keep the current backend as a fallback until all low-order closed forms,
   arithmetic-order reductions, mixed-power fixtures, boundary cases and
   independent reference checks pass. Add larger tests only when measured
   accuracy and runtime support them. Consider compiled MPFR linear algebra
   after profiling the algorithm; the current R elimination can otherwise
   hide improvements in continuation.

No universal pool-count or exponent limit follows from either paper.
Distinct moment count, largest exponent, exponent spacing, proximity to the
moment-space boundary and input precision all matter. A successful 23-moment
fixture would be evidence for that fixture, not a general 23-moment guarantee.

## External experiments

The extension's complete local test suite passed **155/155 assertions** under
Julia 1.11.7. Additional tests used supplied moments at 512-bit precision;
the beta example was repeated at 768 bits. These trials checked admissibility,
all original moment equations, and known nodes/weights where available.

| Experiment | Result |
|---|---|
| Uniform moments at orders 2,5 | Both Radau rules recovered; the lower rule agrees with its closed form |
| Theory beta example, orders 10,40,50 | Unscaled initialization fails; power scaling plus QR and right anchoring recovers both rules |
| Same beta example, 512 versus 768 bits | Maximum node/weight difference below 1.5e-153; prevalence interval [0.90020594951424%, 2.8513322822592%] |
| Known lower rule, orders 10,40,50, nodes and weights (1/3,2/3) | Default lower path fails; unscaled QR with right anchoring recovers both rules; known lower node/weight error below 6e-118 |
| Same known rule after power scaling | Lower rule recovered, but the final Lobatto solve fails at the requested tolerance |
| Known lower rule, orders 100,101,201 | Lower solver still fails in the tested configurations; obtaining only the upper rule is insufficient for our API |
| 23 selected moments from five size-2 and five size-3 pools | No validated pair obtained in the bounded trials, with/without QR, power scaling, strict or hybrid continuation; no improvement over our unresolved case established |

The last fixture has orders {2,3,...,23,25}, twelve known lower nodes i/13,
and weights i/78 for i=1,...,12. It is the same mathematical fixture used
in the earlier R stress assessment. External trial budgets were 90, 120 or 180
seconds for whole Julia processes, including loading and JIT compilation;
some processes needed forced stops after termination signals did not exit
promptly. Multiple experiments overlapped. These are bounded exploratory
trials, not controlled speed comparisons or convergence impossibility claims.
The very tight 1e-100 tolerance in some trials was deliberately chosen to
expose inaccurate recovery in tiny-moment examples; it is not a recommended
default for our 1e-8 accuracy API.

Results and reproduction instructions are recorded alongside the optional
[Julia benchmark](benchmarks/generalized-gauss/README.md). These experiments
do not add Julia to the R package's dependencies or normal CI. The observed
failures are reasons to adapt the method with our safeguards and retain the
existing backend, rather than porting the fork without validation.
