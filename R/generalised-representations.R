generalised_double <- function(moments, orders,
                               tol = 1e-9, maxit = 300L) {
  valid <- function(x) {
    is.numeric(x) && !is.complex(x) &&
      is.null(dim(x)) && length(x) > 0L && all(is.finite(x))
  }
  if (!valid(moments) || !valid(orders) || length(moments) != length(orders) ||
    any(orders <= 0) || any(diff(orders) <= 0)) {
    stop("Supply finite moments and equally long, strictly increasing positive orders.",
      call. = FALSE
    )
  }
  if (!valid(tol) || length(tol) != 1L || tol <= 0 || tol >= 1 ||
    !valid(maxit) || length(maxit) != 1L || maxit < 1 || maxit != floor(maxit)) {
    stop("Invalid `tol` or `maxit`.", call. = FALSE)
  }
  if (any(moments < 0 | moments > 1) || any(diff(moments) > tol)) {
    stop("Infeasible power moments on [0,1].", call. = FALSE)
  }
  d <- length(orders)
  point <- moments[1]^(1 / orders[1])
  predicted <- point^orders
  is_point <- d >= 2L && all(abs(moments - predicted) <=
    16 * .Machine$double.eps * pmax(abs(moments), abs(predicted), 1e-300))
  if (all(moments == moments[1]) || is_point) {
    rule <- if (is_point) {
      data.frame(nodes = point, weights = 1)
    } else {
      data.frame(nodes = c(0, 1), weights = c(1 - moments[1], moments[1]))
    }
    # One constraint alone identifies a unique measure only at endpoints.
    if (d == 1L && moments[1] > 0 && moments[1] < 1) {
      result <- principal_representations(moments, tol = tol)
      result$lower$nodes <- result$lower$nodes^(1 / orders[1])
    } else {
      rule <- rule[rule$weights > 0, , drop = FALSE]
      rownames(rule) <- NULL
      result <- list(
        lower = rule, upper = rule, moments = moments,
        support = c(0, 1), order = d, boundary = TRUE
      )
    }
  } else if (all(orders / orders[1] == seq_len(d))) {
    result <- principal_representations(moments, tol = tol)
    result$lower$nodes <- result$lower$nodes^(1 / orders[1])
    result$upper$nodes <- result$upper$nodes^(1 / orders[1])
  } else {
    target <- c(1, moments)
    exponents <- c(0, orders / orders[1])
    # Work in Y=X^orders[1], avoiding fractional derivatives at zero.
    lower <- selected_power_rule(target, exponents, FALSE, tol, maxit)
    upper <- selected_power_rule(target, exponents, TRUE, tol, maxit)
    lower$nodes <- lower$nodes^(1 / orders[1])
    upper$nodes <- upper$nodes^(1 / orders[1])
    result <- list(
      lower = lower, upper = upper, moments = moments,
      support = c(0, 1), order = d, boundary = FALSE
    )
  }
  residual <- function(rule) {
    max(abs(drop(outer(
      c(0, orders), rule$nodes,
      function(k, x) x^k
    ) %*%
      rule$weights) - c(1, moments)))
  }
  errors <- c(lower = residual(result$lower), upper = residual(result$upper))
  if (any(errors > tol)) stop("Power moment reproduction failed.", call. = FALSE)
  result$orders <- orders
  result$diagnostics$residuals <- errors
  result$boundary <- isTRUE(all.equal(result$lower, result$upper,
    tolerance = tol * 10
  ))
  class(result) <- c("generalised_principal_representations", "principal_representations")
  result
}

selected_power_rule <- function(target, exponents, upper, tol, maxit) {
  d <- length(exponents) - 1L
  r <- if (d %% 2L) (d + 1L) / 2L - as.integer(upper) else d / 2L
  ends <- if (d %% 2L) {
    if (upper) c(0, 1) else numeric()
  } else {
    if (upper) 1 else 0
  }
  n <- r + length(ends)
  # Parameters: all weights followed by the r unfixed nodes.
  unpack <- function(z) c(z[n + seq_len(r)], ends)
  basis <- function(x) outer(exponents, x, function(k, y) y^k)
  scale <- pmax(target, 1e-4)
  for (power in c(1, 2, 0.5, 4, 0.25)) {
    x <- (seq_len(r) / (r + 1))^power
    z <- c(rep(1 / n, n), x)
    for (iteration in seq_len(maxit)) {
      x <- unpack(z)
      f <- drop(basis(x) %*% z[seq_len(n)]) - target
      if (max(abs(f)) <= tol * 0.05) break
      deriv <- outer(exponents, x[seq_len(r)], function(k, y) {
        ifelse(k == 0, 0, k * y^pmax(k - 1, 0))
      })
      jac <- cbind(basis(x), sweep(deriv, 2, z[seq_len(r)], `*`)) / scale
      sv <- svd(jac)
      # Levenberg damping with line search, retaining bounds on all parameters.
      improved <- FALSE
      for (damping in c(0, 1e-10, 1e-7, 1e-4, 0.01)) {
        factors <- sv$d / (sv$d^2 + max(damping, 1e-28))
        step <- -drop(sv$v %*% (factors * crossprod(sv$u, f / scale)))
        for (alpha in 2^-(0:25)) {
          trial <- pmax(0, pmin(1, z + alpha * step))
          f_new <- drop(basis(unpack(trial)) %*% trial[seq_len(n)]) - target
          if (sum((f_new / scale)^2) < sum((f / scale)^2)) {
            z <- trial
            improved <- TRUE
            break
          }
        }
        if (improved) break
      }
      if (!improved) break
    }
    x <- unpack(z)
    w <- z[seq_len(n)]
    # Retain genuine small masses; drop zero-weight nodes and merge duplicates.
    idx <- order(x)
    rule <- data.frame(nodes = x[idx], weights = w[idx])
    rule <- rule[rule$weights > 0, , drop = FALSE]
    rule <- stats::aggregate(weights ~ nodes, rule, sum)
    rule <- rule[, c("nodes", "weights")]
    rownames(rule) <- NULL
    if (abs(sum(rule$weights) - 1) <= tol) {
      rule$weights <- rule$weights / sum(rule$weights)
      if (max(abs(drop(basis(rule$nodes) %*% rule$weights) - target)) <= tol) {
        return(rule)
      }
    }
  }
  stop(
    paste(
      "Generalised quadrature did not converge at `tol`; moments may be",
      "infeasible or ill-conditioned. Try rescaling orders or increasing maxit."
    ),
    call. = FALSE
  )
}

#' Sharp bounds for another power moment
#'
#' @param representations Output of [generalised_principal_representations()].
#' @param exponent Finite positive target exponent.
#' @return Named vector `minimum`, `maximum` for E(X^exponent).
#' @details The orientation is (-1)^k where k is the number of supplied orders
#'   greater than the target exponent. An already supplied order is identified.
#'   In the mixed-pool theory X=1-P, so prevalence bounds are
#'   `c(minimum = 1 - bounds[2], maximum = 1 - bounds[1])` at exponent one.
#' @export
power_moment_bounds <- function(representations, exponent) {
  if (!inherits(representations, "generalised_principal_representations") ||
    !is.numeric(exponent) || is.complex(exponent) || length(exponent) != 1L ||
    !is.finite(exponent) || exponent <= 0) {
    stop("Supply generalised representations and a positive finite exponent.", call. = FALSE)
  }
  if (identical(representations$diagnostics$status, "numerical_boundary")) {
    stop("Unresolved boundary membership: cannot claim sharp bounds at input precision.", call. = FALSE)
  }
  values <- evaluate_representations(representations, function(x) x^exponent)
  if (sum(representations$orders > exponent) %% 2L) values <- rev(values)
  stats::setNames(unname(values), c("minimum", "maximum"))
}
