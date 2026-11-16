#' Compute the two principal representations
#'
#' Construct probability measures on a compact interval matching truncated raw
#' moments. The zeroth moment is one and must not be included in `moments`.
#'
#' @param moments Numeric vector containing E(X), ..., E(X^v), with v >= 1.
#' @param support Finite increasing endpoints of the support interval.
#' @param tol Positive absolute tolerance on normalised moments, used for
#'   feasibility, numerical rank and moment reproduction (default 1e-12).
#'
#' @details For odd v = 2n-1, the lower rule is Gauss and the upper rule is
#'   Gauss-Lobatto. For even v = 2n, the rules are left and right Gauss-Radau.
#'   Boundary sequences have a unique representing measure, returned twice.
#'   A Cholesky-whitened multiplication matrix supplies interior quadrature
#'   nodes and weights. Singular localising matrices supply boundary support.
#'
#'   Raw moment inversion is ill-conditioned at high order or near the boundary.
#'   Results are checked against all supplied normalised moments. A diagnostic
#'   warning is issued for poorly conditioned interior problems; unresolved
#'   numerical problems cause an error rather than silent moment truncation.
#'   No noisy-moment projection is performed. Classification is numerical at
#'   `tol`; it is not an exact-arithmetic certificate of boundary membership.
#'
#' @return An object of class `principal_representations` with `lower` and
#'   `upper` data frames (`nodes`, `weights`), original `moments`, `support`,
#'   `order`, `boundary`, and `diagnostics` (normalised residuals, smallest
#'   feasibility eigenvalue, and minimum localising-matrix eigenvalue ratio).
#' @references Karlin and Studden (1966), Tchebycheff Systems.
#'   Gautschi (2004), Orthogonal Polynomials: Computation and Approximation.
#' @export
#' @examples
#' # First three moments of Uniform(0, 1)
#' reps <- principal_representations(1 / (2:4))
#' reps$lower
#' reps$upper
principal_representations <- function(moments, support = c(0, 1), tol = 1e-12) {
  if (!is.numeric(moments) || is.complex(moments) || !is.null(dim(moments)) ||
      length(moments) < 1L || any(!is.finite(moments))) {
    stop("`moments` must be a non-empty finite numeric vector of m1, ..., mv.",
         call. = FALSE)
  }
  if (!is.numeric(support) || is.complex(support) || length(support) != 2L ||
      any(!is.finite(support)) || support[1] >= support[2]) {
    stop("`support` must contain two finite increasing endpoints.", call. = FALSE)
  }
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L || !is.finite(tol) ||
      tol <= 0 || tol >= 1) {
    stop("`tol` must be a finite number strictly between zero and one.",
         call. = FALSE)
  }
  v <- length(moments)
  m <- normalise_moments(c(1, moments), support)
  if (any(!is.finite(m))) {
    stop("Moment scaling overflowed; rescale the problem.", call. = FALSE)
  }
  checks <- localising_matrices(m)
  spectra <- lapply(checks, function(z) eigen(z$matrix, symmetric = TRUE))
  smallest <- min(vapply(spectra, function(z) min(z$values), numeric(1)))
  if (smallest < -tol) {
    stop("Infeasible moments on the specified support (negative localising eigenvalue).",
         call. = FALSE)
  }
  ratios <- vapply(spectra, function(z) {
    max(0, min(z$values)) / max(max(z$values), .Machine$double.eps)
  }, numeric(1))
  boundary <- smallest <= tol
  if (boundary) {
    j <- which(vapply(spectra, function(z) min(z$values), numeric(1)) <= tol)[1]
    z <- spectra[[j]]
    polynomial <- z$vectors[, which.min(z$values)]
    # A null localising polynomial vanishes wherever its multiplier is positive.
    nodes <- c(real_roots(polynomial, tol), checks[[j]]$endpoints)
    nodes <- sort(nodes[nodes >= -tol & nodes <= 1 + tol])
    nodes <- pmax(0, pmin(1, nodes))
    nodes <- nodes[c(TRUE, diff(nodes) > sqrt(.Machine$double.eps))]
    if (!length(nodes)) numerical_failure()
    design <- outer(0:v, nodes, function(k, x) x^k)
    weights <- tryCatch(qr.solve(design, m, tol = tol),
                        error = function(e) numerical_failure())
    weights[abs(weights) <= tol] <- 0
    lower <- upper <- checked_rule(nodes, weights, m, tol)
  } else {
    if (min(ratios) < 1e-9) {
      warning("Raw moment matrices are poorly conditioned; inspect residuals and sensitivity to `tol`.",
              call. = FALSE)
    }
    if (v %% 2L == 1L) {
      n <- (v + 1L) / 2L
      g <- gaussian_rule(m, n)
      lower <- checked_rule(g$nodes, g$weights, m, tol)
      r <- n - 1L
      if (r > 0L) {
        t <- m[2:(2 * r + 1)] - m[3:(2 * r + 2)]
        g <- gaussian_rule(t, r)
        w <- g$weights / (g$nodes * (1 - g$nodes))
        right <- m[2] - sum(w * g$nodes)
        upper <- checked_rule(c(0, g$nodes, 1),
                              c(1 - sum(w) - right, w, right), m, tol)
      } else {
        upper <- checked_rule(c(0, 1), c(1 - m[2], m[2]), m, tol)
      }
    } else {
      n <- v / 2L
      g <- gaussian_rule(m[2:(2 * n + 1)], n)
      w <- g$weights / g$nodes
      lower <- checked_rule(c(0, g$nodes), c(1 - sum(w), w), m, tol)
      g <- gaussian_rule(m[1:(2 * n)] - m[2:(2 * n + 1)], n)
      w <- g$weights / (1 - g$nodes)
      upper <- checked_rule(c(g$nodes, 1), c(w, 1 - sum(w)), m, tol)
    }
  }
  residual <- function(rule) {
    max(abs(drop(outer(0:v, rule$nodes, function(k, x) x^k) %*%
                   rule$weights) - m))
  }
  diagnostics <- list(
    residuals = c(lower = residual(lower), upper = residual(upper)),
    min_eigenvalue = smallest, min_eigenvalue_ratio = min(ratios)
  )
  for (name in c("lower", "upper")) {
    rule <- get(name)
    rule$nodes <- support[1] + diff(support) * rule$nodes
    assign(name, rule)
  }
  structure(list(lower = lower, upper = upper, moments = moments,
                 support = support, order = v, boundary = boundary,
                 diagnostics = diagnostics), class = "principal_representations")
}

normalise_moments <- function(m, support) {
  a <- support[1]
  width <- diff(support)
  vapply(0:(length(m) - 1L), function(k) {
    i <- 0:k
    sum(choose(k, i) * (-a)^(k - i) * m[i + 1L]) / width^k
  }, numeric(1))
}

moment_matrix <- function(m, n, shift = 0L) {
  outer(0:(n - 1L), 0:(n - 1L), function(i, j) m[i + j + shift + 1L])
}

localising_matrices <- function(m) {
  v <- length(m) - 1L
  if (v %% 2L == 0L) {
    n <- v / 2L
    list(list(matrix = moment_matrix(m, n + 1L), endpoints = numeric()),
         list(matrix = moment_matrix(m, n, 1L) - moment_matrix(m, n, 2L),
              endpoints = c(0, 1)))
  } else {
    n <- (v + 1L) / 2L
    list(list(matrix = moment_matrix(m, n, 1L), endpoints = 0),
         list(matrix = moment_matrix(m, n) - moment_matrix(m, n, 1L),
              endpoints = 1))
  }
}

# Golub-Welsch in a Cholesky-orthonormalised monomial basis.
# This avoids polynomial roots and Vandermonde solves for interior sequences.
gaussian_rule <- function(m, n) {
  # Avoid an unnecessary Cholesky/eigen round trip for a one-node rule.
  # A one-ulp node error near one is amplified by endpoint-singular integrands.
  if (n == 1L) return(list(nodes = m[2] / m[1], weights = m[1]))
  h <- moment_matrix(m, n)
  h1 <- moment_matrix(m, n, 1L)
  inverse <- tryCatch(backsolve(chol(h), diag(n)),
                      error = function(e) numerical_failure())
  jacobi <- crossprod(inverse, h1 %*% inverse)
  eig <- eigen((jacobi + t(jacobi)) / 2, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes = eig$values[ord], weights = m[1] * eig$vectors[1, ord]^2)
}

real_roots <- function(coefficients, tol) {
  while (length(coefficients) > 1L && abs(tail(coefficients, 1)) < tol) {
    coefficients <- head(coefficients, -1L)
  }
  if (length(coefficients) == 1L) return(numeric())
  roots <- polyroot(coefficients)
  Re(roots[abs(Im(roots)) < sqrt(tol)])
}

checked_rule <- function(nodes, weights, m, tol) {
  if (any(!is.finite(c(nodes, weights))) || any(nodes < -tol) ||
      any(nodes > 1 + tol) || any(weights < -tol)) numerical_failure()
  nodes <- pmax(0, pmin(1, nodes))
  weights <- pmax(0, weights)
  weights <- weights / sum(weights)
  residual <- max(abs(drop(outer(0:(length(m) - 1L), nodes,
                                 function(k, x) x^k) %*% weights) - m))
  if (!is.finite(residual) || residual > tol) numerical_failure()
  keep <- weights > 0
  data.frame(nodes = nodes[keep], weights = weights[keep])
}

numerical_failure <- function() {
  stop(paste("Moment inversion could not reproduce all constraints at `tol`.",
             "The sequence may be infeasible or numerically ill-conditioned;",
             "rescale, reduce order, or assess tolerance sensitivity."), call. = FALSE)
}
