test_that("polynomials through degree v are identified", {
  reps <- principal_representations(1 / (2:5))
  for (k in 1:4) {
    values <- evaluate_representations(reps, function(x) x^k)
    expect_equal(unname(values), rep(1 / (k + 1), 2))
  }
})

test_that("negative derivative reverses the bounds", {
  reps <- principal_representations(1 / (2:4))
  positive <- functional_bounds(reps, exp, derivative = "positive")
  negative <- functional_bounds(reps, function(x) -exp(x), derivative = "negative")
  expect_equal(unname(negative), -rev(unname(positive)))
  expect_error(functional_bounds(reps, exp), "Explicitly")
  expect_error(functional_bounds(reps, exp, derivative = "negative"), "contradict")
})

test_that("the pooled prevalence integrand agrees with closed-form bounds", {
  m1 <- 0.15
  m2 <- 0.06
  s <- 25
  f <- function(q, s) -expm1(log1p(-q) / s)
  reps <- principal_representations(c(m1, m2))
  bounds <- functional_bounds(reps, f, s = s, derivative = "positive")
  b <- m2 / m1
  y <- (1 - 2 * m1 + m2) / (1 - m1)
  expect_equal(unname(bounds), c(m1 / b * f(b, s),
                                1 - (1 - m1) * y^(1 / s - 1)))
})

test_that("distribution functionals and arguments are supported", {
  reps <- principal_representations(0.3)
  variance <- function(nodes, weights, multiplier) {
    multiplier * (sum(weights * nodes^2) - sum(weights * nodes)^2)
  }
  values <- evaluate_representations(reps, variance, multiplier = 2,
                                     type = "distribution")
  expect_equal(unname(values), c(0, 0.42))
})

test_that("malformed function results are rejected", {
  reps <- principal_representations(0.3)
  expect_error(evaluate_representations(NULL, exp), "representations")
  expect_error(evaluate_representations(reps, 1), "function")
  expect_error(evaluate_representations(reps, function(x) NA_real_), "finite")
  expect_error(evaluate_representations(reps, function(x) 1), "per node")
  expect_error(evaluate_representations(reps, function(x, w) c(1, 2),
                                        type = "distribution"), "one finite")
})
