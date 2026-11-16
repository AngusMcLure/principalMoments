# Analytic fixtures from the one-, two- and three-moment formulas in the
# clustered group-testing theory note; independent of the inversion backend.
prevalence_integrand <- function(q, s) -expm1(log1p(-q) / s)

expect_prevalence_bounds <- function(reps, lower, upper) {
  for (s in c(1, 2, 5, 25, 100)) {
    expected <- c(
      sum(lower$weights * prevalence_integrand(lower$nodes, s)),
      sum(upper$weights * prevalence_integrand(upper$nodes, s))
    )
    if (identical(reps$diagnostics$status, "numerical_boundary")) {
      # Rounded flat sequences determine an approximate boundary rule, but
      # must not be advertised as exact sharp bounds without resolving rank.
      expect_error(functional_bounds(reps, prevalence_integrand,
        s = s,
        derivative = "positive"
      ), "Unresolved boundary")
      actual <- evaluate_representations(reps, prevalence_integrand, s = s)
    } else {
      actual <- functional_bounds(reps, prevalence_integrand,
        s = s,
        derivative = "positive"
      )
    }
    expect_equal(unname(actual), expected,
      tolerance = 1e-8,
      info = paste(
        "moments", paste(reps$moments, collapse = ","),
        "pool size", s
      )
    )
  }
}

expect_closed_rule <- function(actual, nodes, weights) {
  expect_equal(actual$nodes, nodes, tolerance = 1e-8)
  expect_equal(actual$weights, weights, tolerance = 1e-8)
}

test_that("one-moment prevalence bounds hold across rare and common positives", {
  for (m1 in c(1e-8, 1e-4, 0.01, 0.15, 0.5, 0.85, 0.99, 1 - 1e-4, 1 - 1e-8)) {
    reps <- principal_representations(m1)
    expect_prevalence_bounds(
      reps, list(nodes = m1, weights = 1),
      list(nodes = c(0, 1), weights = c(1 - m1, m1))
    )
  }
})

test_that("two-moment closed forms hold throughout the feasible moment region", {
  for (m1 in c(1e-8, 1e-4, 0.01, 0.15, 0.5, 0.85, 0.99, 1 - 1e-4, 1 - 1e-8)) {
    # The feasible interval is m1^2 <= m2 <= m1. Fractions near both
    # ends exercise nearly constant and nearly Bernoulli mixing laws.
    for (fraction in c(0.001, 0.05, 0.5, 0.95, 0.999)) {
      m2 <- m1^2 + fraction * m1 * (1 - m1)
      if ((m1 < 1e-6 || m1 > 1 - 1e-6) && fraction <= 0.05) {
        expect_warning(
          reps <- principal_representations(c(m1, m2)),
          "poorly conditioned"
        )
        expect_lt(reps$diagnostics$min_eigenvalue_ratio, 1e-9)
      } else {
        reps <- principal_representations(c(m1, m2))
      }
      b <- m2 / m1
      a <- (m1 - m2) / (1 - m1)
      lower <- list(nodes = c(0, b), weights = c(1 - m1 / b, m1 / b))
      upper <- list(
        nodes = c(a, 1),
        weights = c((1 - m1) / (1 - a), (m1 - a) / (1 - a))
      )
      expect_false(reps$boundary)
      expect_closed_rule(reps$lower, lower$nodes, lower$weights)
      expect_closed_rule(reps$upper, upper$nodes, upper$weights)
      expect_prevalence_bounds(reps, lower, upper)
    }
  }
})

test_that("point masses including all mass at zero collapse the prevalence bounds", {
  for (v in 1:3) {
    for (point in c(0, 0.2, 0.9, 1)) {
      # One moment identifies a point mass only at endpoints; use >=2
      # moments for the interior point-mass cases.
      if (v == 1 && point > 0 && point < 1) next
      reps <- principal_representations(point^seq_len(v))
      expect_true(reps$boundary)
      expect_identical(reps$lower, reps$upper)
      expect_closed_rule(reps$lower, point, 1)
      rule <- list(nodes = point, weights = 1)
      expect_prevalence_bounds(reps, rule, rule)
    }
  }
})

test_that("Bernoulli boundary moments identify prevalence independently of pool size", {
  for (v in 2:3) {
    for (mass in c(1e-4, 0.01, 0.5, 0.99, 1 - 1e-4)) {
      reps <- principal_representations(rep(mass, v))
      expect_true(reps$boundary)
      expect_identical(reps$lower, reps$upper)
      rule <- list(nodes = c(0, 1), weights = c(1 - mass, mass))
      expect_closed_rule(reps$lower, rule$nodes, rule$weights)
      expect_prevalence_bounds(reps, rule, rule)
    }
  }
})

test_that("three-moment closed forms hold for asymmetric interior moment sequences", {
  for (m1 in c(0.01, 0.15, 0.5, 0.85, 0.99)) {
    for (fraction in c(0.05, 0.3, 0.8, 0.95)) {
      m2 <- m1^2 + fraction * m1 * (1 - m1)
      # Hausdorff bounds on the third moment, given the first two.
      m3_min <- m2^2 / m1
      m3_max <- m2 - (m1 - m2)^2 / (1 - m1)
      for (third_fraction in c(0.1, 0.5, 0.9)) {
        m3 <- m3_min + third_fraction * (m3_max - m3_min)
        reps <- principal_representations(c(m1, m2, m3))
        a <- (m3 - m1 * m2) / (m2 - m1^2)
        b <- a * m1 - m2
        x <- (a + c(-1, 1) * sqrt(a^2 - 4 * b)) / 2
        w <- c((x[2] - m1) / diff(x), (m1 - x[1]) / diff(x))
        c_node <- (m2 - m3) / (m1 - m2)
        w_c <- (m1 - m2) / (c_node * (1 - c_node))
        w_right <- m1 - c_node * w_c
        lower <- list(nodes = x, weights = w)
        upper <- list(
          nodes = c(0, c_node, 1),
          weights = c(1 - w_c - w_right, w_c, w_right)
        )
        expect_false(reps$boundary)
        expect_closed_rule(reps$lower, lower$nodes, lower$weights)
        expect_closed_rule(reps$upper, upper$nodes, upper$weights)
        expect_prevalence_bounds(reps, lower, upper)
      }
    }
  }
})

test_that("third-moment endpoints recover their unique two-atom measures", {
  for (m1 in c(0.01, 0.2, 0.5, 0.8, 0.99)) {
    for (fraction in c(0.05, 0.5, 0.95)) {
      m2 <- m1^2 + fraction * m1 * (1 - m1)
      b <- m2 / m1
      a <- (m1 - m2) / (1 - m1)
      cases <- list(
        list(
          m3 = m2^2 / m1,
          nodes = c(0, b), weights = c(1 - m1 / b, m1 / b)
        ),
        list(
          m3 = m2 - (m1 - m2)^2 / (1 - m1),
          nodes = c(a, 1),
          weights = c((1 - m1) / (1 - a), (m1 - a) / (1 - a))
        )
      )
      for (case in cases) {
        reps <- principal_representations(c(m1, m2, case$m3))
        expect_true(reps$boundary)
        expect_identical(reps$lower, reps$upper)
        expect_closed_rule(reps$lower, case$nodes, case$weights)
        expect_prevalence_bounds(reps, case, case)
      }
    }
  }
})
