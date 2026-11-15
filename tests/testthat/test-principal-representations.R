moments_of <- function(nodes, weights, v) {
  vapply(seq_len(v), function(k) sum(weights * nodes^k), numeric(1))
}

expect_matches <- function(reps, moments, tolerance = 1e-10) {
  for (rule in reps[c("lower", "upper")]) {
    expect_true(all(diff(rule$nodes) > 0))
    expect_true(all(rule$nodes >= reps$support[1]))
    expect_true(all(rule$nodes <= reps$support[2]))
    expect_true(all(rule$weights >= 0))
    expect_equal(sum(rule$weights), 1, tolerance = tolerance)
    expect_equal(moments_of(rule$nodes, rule$weights, length(moments)),
                 moments, tolerance = tolerance)
  }
}

test_that("one-moment rules have the known closed forms", {
  reps <- principal_representations(0.3)
  expect_equal(reps$lower$nodes, 0.3)
  expect_equal(reps$lower$weights, 1)
  expect_equal(reps$upper$nodes, c(0, 1))
  expect_equal(reps$upper$weights, c(0.7, 0.3))
  expect_false(reps$boundary)
})

test_that("two-moment rules agree with analytic Radau formulas", {
  m1 <- 0.3
  m2 <- 0.16
  reps <- principal_representations(c(m1, m2))
  b <- m2 / m1
  a <- (m1 - m2) / (1 - m1)
  expect_equal(reps$lower$nodes, c(0, b))
  expect_equal(reps$lower$weights, c(1 - m1 / b, m1 / b))
  expect_equal(reps$upper$nodes, c(a, 1))
  expect_equal(reps$upper$weights, c((1 - m1) / (1 - a),
                                   (m1 - a) / (1 - a)))
})

test_that("three-moment uniform rules agree with classical quadrature", {
  reps <- principal_representations(1 / (2:4))
  expect_equal(reps$lower$nodes, (1 + c(-1, 1) / sqrt(3)) / 2)
  expect_equal(reps$lower$weights, c(0.5, 0.5))
  expect_equal(reps$upper$nodes, c(0, 0.5, 1))
  expect_equal(reps$upper$weights, c(1, 4, 1) / 6)
})

test_that("all four parity rules reproduce every supplied moment", {
  for (v in 1:10) {
    m <- 1 / (seq_len(v) + 1)
    reps <- principal_representations(m)
    expect_matches(reps, m)
    n <- ceiling(v / 2)
    expect_equal(nrow(reps$lower), if (v %% 2) n else n + 1)
    expect_equal(nrow(reps$upper), n + 1)
    expect_lt(max(reps$diagnostics$residuals), 1e-12)
  }
})

test_that("beta and discrete distributions are enclosed for positive derivatives", {
  set.seed(47)
  for (v in 1:8) {
    # Beta(2, 3) moments and exact MGF integral computed independently
    k <- seq_len(v)
    m <- exp(lbeta(2 + k, 3) - lbeta(2, 3))
    reps <- principal_representations(m)
    expect_matches(reps, m)
    bounds <- functional_bounds(reps, exp, derivative = "positive")
    truth <- integrate(function(x) exp(x) * dbeta(x, 2, 3), 0, 1)$value
    expect_lte(bounds[1], truth + 1e-10)
    expect_gte(bounds[2], truth - 1e-10)
    x <- sort(runif(20))
    w <- rexp(20)
    w <- w / sum(w)
    m <- moments_of(x, w, v)
    reps <- principal_representations(m)
    expect_matches(reps, m)
    bounds <- functional_bounds(reps, exp, derivative = "positive")
    expect_lte(bounds[1], sum(w * exp(x)) + 1e-10)
    expect_gte(bounds[2], sum(w * exp(x)) - 1e-10)
  }
})

test_that("boundary measures are recovered with all constraints checked", {
  cases <- list(list(x = 0, w = 1), list(x = 1, w = 1),
                list(x = 0.3, w = 1), list(x = c(0, 1), w = c(0.7, 0.3)),
                list(x = c(0.2, 0.8), w = c(0.4, 0.6)),
                list(x = c(0, 0.4, 1), w = c(0.2, 0.5, 0.3)))
  for (case in cases) {
    for (v in 6:7) {
      m <- moments_of(case$x, case$w, v)
      reps <- principal_representations(m)
      expect_true(reps$boundary)
      expect_equal(reps$lower, reps$upper)
      expect_matches(reps, m)
      expect_equal(reps$lower$nodes, case$x, tolerance = 1e-7)
      expect_equal(reps$lower$weights, case$w, tolerance = 1e-7)
    }
  }
  expect_true(principal_representations(0)$boundary)
  expect_true(principal_representations(1)$boundary)
  expect_true(principal_representations(c(0.4, 0.16))$boundary)
  expect_true(principal_representations(c(0.4, 0.4))$boundary)
})

test_that("general compact supports transform correctly", {
  # Uniform(-2, 3), tested without using the package's transformation
  for (v in 1:6) {
    k <- seq_len(v)
    m <- (3^(k + 1) - (-2)^(k + 1)) / (5 * (k + 1))
    reps <- principal_representations(m, support = c(-2, 3))
    expect_matches(reps, m, tolerance = 1e-8)
  }
})

test_that("invalid and infeasible moments fail explicitly", {
  for (m in list(numeric(), NA_real_, Inf, "a", matrix(1))) {
    expect_error(principal_representations(m), "moments")
  }
  for (s in list(c(1, 0), c(0, Inf), 0)) {
    expect_error(principal_representations(0.5, support = s), "support")
  }
  for (tol in list(0, -1, NA_real_, Inf, c(1, 2))) {
    expect_error(principal_representations(0.5, tol = tol), "tol")
  }
  for (m in list(-0.1, 1.1, c(0.5, 0.1), c(0.2, 0.3),
                 c(0.5, 1 / 3, 0.01))) {
    expect_error(principal_representations(m), "Infeasible")
  }
})
