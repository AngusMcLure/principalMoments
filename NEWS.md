# principalMoments 0.3.0

* Add decimal-string and MPFR inputs, adaptive working-precision refinement,
  high-precision orthogonal-polynomial recurrences, refined nodes and Christoffel weights.
* Replace final mixed-power fitting with positive ordered MPFR Newton corrections
  and exponent/moment continuation. Double fits are warm starts only.
* Fail safely on unresolved long double-input ranks and uncertain sharp-bound
  claims; diagnose local input sensitivity separately from arithmetic precision.
* Check relative moment errors and functional stability; retain MPFR values
  beyond the double-precision display tables.
* Add independent Jacobi, former-failure and high-power regressions, with a
  separate manually triggered numerical stress workflow.

# principalMoments 0.2.0

* Add selected-power principal representations and power-moment bounds on
  [0,1], including the mixed-pool-size theory examples and equal-size reductions.

# principalMoments 0.1.0

* Initial principal representations, feasibility diagnostics and functional evaluation.

* Expand pooled-prevalence closed-form tests across one to three moments,
  extreme means, point masses, Bernoulli boundaries and third-moment endpoints.
* Use the direct one-node Gaussian formula to avoid amplified rounding errors
  for prevalence integrands near the right endpoint.
