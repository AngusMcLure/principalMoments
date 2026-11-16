# Internal multiprecision arithmetic. MPFR values are retained until display.
is_mp <- function(x) inherits(x, "mpfr")
mp_num <- function(x) if (is_mp(x)) Rmpfr::asNumeric(x) else x
mp_value <- function(x, bits) {
  if (is_mp(x)) Rmpfr::roundMpfr(x, bits) else Rmpfr::mpfr(x, bits)
}
mp_zeros <- function(n, bits) Rmpfr::mpfr(rep(0, n), bits)
mp_outer <- function(x, y, fun = `*`) Rmpfr::outer(x, y, FUN = fun)

# Pivoted elimination supports vector or matrix RHS, without converting to double.
mp_solve <- function(a, b) {
  n <- nrow(a)
  vector <- is.null(dim(b))
  if (vector) dim(b) <- c(n, 1L)
  for (k in seq_len(n)) {
    pivot <- k - 1L + which(abs(a[k:n, k]) == max(abs(a[k:n, k])))[1L]
    if (a[pivot, k] == 0) stop("Unresolved singular moment system.", call. = FALSE)
    if (pivot != k) {
      a[c(k, pivot), ] <- a[c(pivot, k), ]
      b[c(k, pivot), ] <- b[c(pivot, k), ]
    }
    if (k < n) {
      for (i in (k + 1L):n) {
        factor <- a[i, k] / a[k, k]
        a[i, k:n] <- a[i, k:n] - factor * a[k, k:n]
        b[i, ] <- b[i, ] - factor * b[k, ]
      }
    }
  }
  x <- b * 0
  for (k in n:1L) {
    value <- b[k, ]
    if (k < n) {
      value <- value - a[k, (k + 1L):n, drop = FALSE] %*%
        x[(k + 1L):n, , drop = FALSE]
    }
    x[k, ] <- value / a[k, k]
  }
  if (vector) methods::as(x, "mpfr") else x
}

# Stieltjes recurrence from moments, evaluated entirely at working precision.
mp_gauss <- function(m, n, bits) {
  if (!n) {
    return(list(nodes = mp_zeros(0, bits), weights = mp_zeros(0, bits)))
  }
  if (n == 1L) {
    return(list(nodes = m[2] / m[1], weights = m[1]))
  }
  if (m[1] <= 0) stop("Infeasible or unresolved modified moments.", call. = FALSE)
  inner <- function(p, q, shift = 0L) {
    indices <- outer(seq_along(p), seq_along(q), `+`) - 1L + shift
    sum(mp_outer(p, q) * m[indices])
  }
  alpha <- beta <- mp_zeros(n, bits)
  p <- mp_value(1, bits) / sqrt(m[1])
  previous <- mp_zeros(1, bits)
  for (j in seq_len(n)) {
    alpha[j] <- inner(p, p, 1L)
    if (j < n) {
      xp <- c(mp_value(0, bits), p)
      padded <- c(p, mp_value(0, bits))
      old <- c(previous, mp_zeros(length(xp) - length(previous), bits))
      next_p <- xp - alpha[j] * padded - beta[j] * old
      norm2 <- inner(next_p, next_p)
      if (norm2 <= 0) {
        stop("Infeasible or unresolved recurrence; use more accurate moments.",
          call. = FALSE
        )
      }
      beta[j + 1L] <- sqrt(norm2)
      previous <- p
      p <- next_p / beta[j + 1L]
    }
  }
  jac <- diag(mp_num(alpha), n)
  for (j in 2:n) jac[j, j - 1L] <- jac[j - 1L, j] <- mp_num(beta[j])
  nodes <- mp_value(sort(eigen(jac, symmetric = TRUE, only.values = TRUE)$values), bits)
  evaluate <- function(x) {
    a <- da <- mp_value(0, bits)
    b <- 1 / sqrt(m[1])
    db <- mp_value(0, bits)
    squares <- b * b
    for (j in seq_len(n)) {
      f <- (x - alpha[j]) * b - beta[j] * a
      df <- b + (x - alpha[j]) * db - beta[j] * da
      if (j < n) {
        a <- b
        da <- db
        b <- f / beta[j + 1L]
        db <- df / beta[j + 1L]
        squares <- squares + b * b
      }
    }
    list(f = f, df = df, weight = 1 / squares)
  }
  weights <- mp_zeros(n, bits)
  eps <- mp_value(2, bits)^(-bits + 24L)
  for (i in seq_len(n)) {
    for (j in 1:30) {
      val <- evaluate(nodes[i])
      correction <- val$f / val$df
      nodes[i] <- nodes[i] - correction
      if (abs(correction) < eps) break
    }
    weights[i] <- evaluate(nodes[i])$weight
  }
  if (any(nodes <= 0 | nodes >= 1) || any(diff(nodes) <= 0) || any(weights <= 0)) {
    stop("Unresolved quadrature support at working precision.", call. = FALSE)
  }
  list(nodes = nodes, weights = weights)
}

mp_consecutive <- function(m, bits) {
  d <- length(m) - 1L
  if (d %% 2L) {
    n <- (d + 1L) / 2L
    lower <- mp_gauss(m, n, bits)
    r <- n - 1L
    t <- if (r) m[2:(2 * r + 1)] - m[3:(2 * r + 2)] else mp_zeros(0, bits)
    upper <- mp_gauss(t, r, bits)
    w <- upper$weights / (upper$nodes * (1 - upper$nodes))
    right <- m[2] - sum(w * upper$nodes)
    upper <- list(
      nodes = c(mp_value(0, bits), upper$nodes, mp_value(1, bits)),
      weights = c(1 - sum(w) - right, w, right)
    )
  } else {
    n <- d / 2L
    lower <- mp_gauss(m[2:(2 * n + 1)], n, bits)
    w <- lower$weights / lower$nodes
    lower <- list(nodes = c(mp_value(0, bits), lower$nodes), weights = c(1 - sum(w), w))
    upper <- mp_gauss(m[1:(2 * n)] - m[2:(2 * n + 1)], n, bits)
    w <- upper$weights / (1 - upper$nodes)
    upper <- list(nodes = c(upper$nodes, mp_value(1, bits)), weights = c(w, 1 - sum(w)))
  }
  list(lower = lower, upper = upper)
}

mp_basis <- function(orders, nodes) mp_outer(orders, nodes, function(k, x) x^k)

# Newton correction with ordered support and positive weights. No clipping.
mp_correct <- function(z, orders, target, ends, r, bits, maxit, budget = NULL, accuracy = 1e-8) {
  n <- length(z) - r
  scale <- abs(target)
  scale[scale == 0] <- mp_value(2, bits)^(-bits / 2)
  unpack <- function(z) c(z[n + seq_len(r)], ends)
  admissible <- function(z) {
    x <- z[n + seq_len(r)]
    all(z[seq_len(n)] > 0) && all(x > 0 & x < 1) && all(diff(x) > 0)
  }
  goal <- mp_value(accuracy, bits) / 100
  for (iteration in seq_len(maxit)) {
    if (!is.null(budget)) {
      budget$remaining <- budget$remaining - 1L
      if (budget$remaining < 0) stop("Continuation work budget exhausted; increase maxit or reduce dimension.", call. = FALSE)
    }
    x <- unpack(z)
    b <- mp_basis(orders, x)
    f <- methods::as(b %*% z[seq_len(n)], "mpfr") - target
    derivative <- mp_outer(orders, x[seq_len(r)], function(k, y) {
      k * y^Rmpfr::pmax(k - 1, 0)
    })
    for (j in seq_len(r)) derivative[, j] <- derivative[, j] * z[j]
    jac <- Rmpfr::cbind(b, derivative)
    logarithmic <- min(z) < mp_value("1e-8", bits)
    if (logarithmic) for (j in seq_along(z)) jac[, j] <- jac[, j] * z[j]
    step <- tryCatch(mp_solve(jac, -f), error = function(e) NULL)
    if (is.null(step)) {
      return(NULL)
    }
    if (max(abs(f / scale)) < goal && max(abs(step)) < goal) {
      return(z)
    }
    improved <- FALSE
    alpha <- mp_value(1, bits)
    if (logarithmic && max(abs(step)) > 10) alpha <- 10 / max(abs(step))
    for (j in 1:70) {
      candidate <- if (logarithmic) z * exp(alpha * step) else z + alpha * step
      if (admissible(candidate)) {
        new_f <- methods::as(mp_basis(orders, unpack(candidate)) %*% candidate[seq_len(n)], "mpfr") - target
        if (max(abs(new_f / scale)) < max(abs(f / scale))) {
          z <- candidate
          improved <- TRUE
          break
        }
      }
      alpha <- alpha / 2
    }
    if (!improved) {
      return(NULL)
    }
  }
  NULL
}

mp_selected <- function(moments, orders, bits, maxit, warm = NULL, accuracy = 1e-8) {
  d <- length(orders)
  e <- c(mp_value(0, bits), mp_value(orders, bits) / orders[1])
  target <- c(mp_value(1, bits), moments)
  start <- NULL
  out <- list()
  budget <- new.env(parent = emptyenv())
  budget$remaining <- 4L * maxit
  for (side in c("lower", "upper")) {
    upper <- side == "upper"
    r <- if (d %% 2L) (d + 1L) / 2L - as.integer(upper) else d / 2L
    ends <- if (d %% 2L) {
      if (upper) c(0, 1) else numeric()
    } else {
      as.numeric(upper)
    }
    ends <- mp_value(ends, bits)
    if (!is.null(warm)) {
      candidate <- warm[[side]]
      interior <- which(candidate$nodes > 0 & candidate$nodes < 1)
      endpoint <- which(candidate$nodes == 0 | candidate$nodes == 1)
      if (length(interior) == r && length(candidate$weights) == r + length(ends) &&
        all(candidate$weights > 0)) {
        initial <- c(
          mp_value(candidate$weights[c(interior, endpoint)], bits),
          mp_value(candidate$nodes[interior], bits)^orders[1]
        )
        corrected <- mp_correct(initial, e, target, ends, r, bits, min(maxit, 20L), budget, accuracy)
        if (!is.null(corrected)) {
          n <- length(corrected) - r
          nodes <- c(corrected[n + seq_len(r)], ends)
          idx <- order(nodes)
          out[[side]] <- list(
            nodes = nodes[idx]^(1 / mp_value(orders[1], bits)),
            weights = corrected[seq_len(n)][idx]
          )
          next
        }
      }
    }
    if (is.null(start)) start <- mp_consecutive(1 / mp_value(1:(d + 1L), bits), bits)
    rule <- start[[side]]
    interiors <- which(rule$nodes > 0 & rule$nodes < 1)
    endpoints <- which(rule$nodes == 0 | rule$nodes == 1)
    z <- c(rule$weights[c(interiors, endpoints)], rule$nodes[interiors])
    # Stage 1: deform consecutive exponents to the selected system for Uniform Y.
    # Stage 2: deform the uniform moment vector to the supplied target vector.
    for (stage in 1:2) {
      t <- 0
      dt <- if (stage == 1) 0.25 else 2
      base_end <- max(16, ceiling(mp_num(-log2(min(target[target > 0])))) + 4)
      end_t <- if (stage == 1) 1 else max(base_end, ceiling(bits * 0.6))
      next_probe <- base_end
      steps <- 0L
      while (t < end_t) {
        steps <- steps + 1L
        if (steps > 256L || dt < 2^-25) {
          stop("Unresolved generalized continuation; input may be ill-conditioned or on a boundary.",
            call. = FALSE
          )
        }
        next_t <- min(end_t, t + dt)
        exponents <- if (stage == 1) (1 - next_t) * mp_value(0:d, bits) + next_t * e else e
        mix <- if (next_t == end_t) mp_value(0, bits) else mp_value(2, bits)^(-next_t)
        wanted <- if (stage == 1) {
          1 / (1 + exponents)
        } else {
          mix / (1 + e) + (1 - mix) * target
        }
        next_z <- mp_correct(z, exponents, wanted, ends, r, bits, maxit, budget, accuracy)
        if (is.null(next_z)) {
          dt <- dt / 2
        } else {
          z <- next_z
          t <- next_t
          if (stage == 2 && t >= next_probe && t < end_t) {
            final_z <- mp_correct(z, e, target, ends, r, bits, min(maxit, 20L), budget, accuracy)
            if (!is.null(final_z)) {
              z <- final_z
              break
            }
            next_probe <- t + 16
          }
          dt <- min(if (stage == 1) 0.5 else 4, dt * 1.5)
        }
      }
    }
    n <- length(z) - r
    nodes <- c(z[n + seq_len(r)], ends)
    idx <- order(nodes)
    out[[side]] <- list(nodes = nodes[idx]^(1 / mp_value(orders[1], bits)), weights = z[seq_len(n)][idx])
  }
  out
}

# First-order sensitivity to one-ulp relative uncertainty in numeric input.
# This is a rejection diagnostic, not a rigorous global error enclosure.
mp_sensitivity <- function(rules, orders, moments, bits) {
  radius <- c(mp_value(0, bits), abs(moments) * mp_value(2, bits)^(-52L))
  estimates <- numeric()
  for (rule in rules) {
    interior <- which(rule$nodes > 0 & rule$nodes < 1)
    n <- length(rule$nodes)
    r <- length(interior)
    if (n + r != length(orders) + 1L) next
    k <- mp_value(c(0, orders), bits)
    jac <- mp_basis(k, rule$nodes)
    derivative <- mp_outer(k, rule$nodes[interior], function(k, x) k * x^(k - 1))
    for (j in seq_along(interior)) derivative[, j] <- derivative[, j] * rule$weights[interior[j]]
    jac <- Rmpfr::cbind(jac, derivative)
    identity <- mp_value(diag(length(k)), bits)
    dim(identity) <- c(length(k), length(k))
    inverse <- mp_solve(jac, identity)
    estimates <- c(estimates, mp_num(max(abs(inverse) %*% radius)))
  }
  if (length(estimates)) max(estimates) else 0
}
