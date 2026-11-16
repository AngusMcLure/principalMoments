#' Evaluate a function at both principal representations
#'
#' @param representations Object returned by [principal_representations()].
#' @param fun For `type = "integrand"`, a vectorised function of nodes returning
#'   one finite value per node. For `type = "distribution"`, a function of
#'   `(nodes, weights, ...)` returning one finite numeric value.
#' @param ... Additional arguments to `fun`.
#' @param type Evaluate the integral of an integrand or a distribution functional.
#' @return Named numeric vector with `lower` and `upper` evaluations. These names
#'   identify the representations, not necessarily the ordering of the values.
#' @details Arbitrary functions are accepted, but their extrema need not be
#'   attained by principal representations. In particular, nonlinear functions
#'   of a distribution (such as variance) require separate extremality arguments.
#' @export
#' @examples
#' reps <- principal_representations(1 / (2:4))
#' evaluate_representations(reps, exp)
evaluate_representations <- function(representations, fun, ...,
                                     type = c("integrand", "distribution")) {
  type <- match.arg(type)
  if (!inherits(representations, "principal_representations")) {
    stop("`representations` must come from principal_representations().", call. = FALSE)
  }
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  rules <- if (is.null(representations$high_precision)) {
    representations[c("lower", "upper")]
  } else {
    requireNamespace("Rmpfr", quietly = TRUE)
    representations$high_precision
  }
  evaluate_rule <- function(rule) {
    if (type == "integrand") {
      values <- fun(rule$nodes, ...)
      if ((!is.numeric(values) && !is_mp(values)) || length(values) != length(rule$nodes) ||
        any(!is.finite(values))) {
        stop("The integrand must return one finite numeric value per node.", call. = FALSE)
      }
      mp_num(sum(rule$weights * values))
    } else {
      value <- fun(rule$nodes, rule$weights, ...)
      if ((!is.numeric(value) && !is_mp(value)) || length(value) != 1L || !is.finite(value)) {
        stop("The distribution functional must return one finite numeric value.",
          call. = FALSE
        )
      }
      mp_num(value)
    }
  }
  values <- vapply(rules, evaluate_rule, numeric(1))
  if (!is.null(representations$precision_reference)) {
    previous <- vapply(representations$precision_reference, evaluate_rule, numeric(1))
    if (any(abs(previous - values) > representations$diagnostics$accuracy)) {
      stop("Functional evaluation is unresolved across precision refinement; increase precision.",
        call. = FALSE
      )
    }
  }
  values
}

#' Sharp bounds under an asserted derivative-sign condition
#'
#' @inheritParams evaluate_representations
#' @param derivative Sign of the (v+1)th derivative of the integrand on the
#'   interior of the support: `"positive"` or `"negative"`. This is a mathematical
#'   assertion supplied by the caller, not checked numerically by the package.
#' @return Named numeric vector `minimum`, `maximum`.
#' @details Requires a continuous integrand on the closed support and a
#'   (v+1)th derivative of the stated sign on its interior, or an equivalent
#'   appropriately oriented Chebyshev-system condition. These conditions are
#'   sufficient, not necessary. A derivative of constant nonnegative or
#'   nonpositive sign also permits equality. Without such a condition, use
#'   [evaluate_representations()] and do not interpret its values as bounds.
#' @export
#' @examples
#' reps <- principal_representations(1 / (2:4))
#' # Every derivative of exp is positive.
#' functional_bounds(reps, exp, derivative = "positive")
functional_bounds <- function(representations, fun, ...,
                              derivative = c("positive", "negative")) {
  if (identical(representations$diagnostics$status, "numerical_boundary")) {
    stop("Unresolved boundary membership: cannot claim sharp bounds at input precision.", call. = FALSE)
  }
  if (inherits(representations, "generalised_principal_representations")) {
    stop("Use power_moment_bounds() for selected powers; the ordinary derivative rule does not apply.",
      call. = FALSE
    )
  }
  if (missing(derivative)) {
    stop("Explicitly assert the derivative sign using `derivative`.", call. = FALSE)
  }
  derivative <- match.arg(derivative)
  values <- evaluate_representations(representations, fun, ...)
  if (derivative == "negative") values <- rev(values)
  if (values[1] > values[2] + 1e-10 * max(1, abs(values))) {
    stop("Evaluations contradict the asserted derivative sign.", call. = FALSE)
  }
  stats::setNames(unname(values), c("minimum", "maximum"))
}
