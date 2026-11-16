power_moments <- function(x, w, orders) {
  vapply(orders, function(k) sum(w * x^k), numeric(1))
}

expect_power_constraints <- function(reps, m, orders, tol = 1e-8) {
  for (rule in reps[c("lower", "upper")]) {
    expect_true(all(rule$weights >= 0))
    expect_equal(sum(rule$weights), 1, tolerance = tol)
    expect_true(all(rule$nodes >= 0 & rule$nodes <= 1))
    expect_equal(power_moments(rule$nodes, rule$weights, orders), m,
      tolerance = tol
    )
  }
}

test_that("selected mixed-pool powers reproduce the theory example", {
  for (orders in list(
    c(2, 3, 5), c(5, 6, 11), c(2, 3, 4, 5),
    c(1.5, 2.8, 4.1), c(2, 3, 5, 7, 9)
  )) {
    m <- 1 / (orders + 1) # Uniform X, hence true prevalence 1/2
    reps <- generalised_principal_representations(m, orders)
    expect_power_constraints(reps, m, orders)
    bounds <- power_moment_bounds(reps, 1)
    expect_lte(bounds[1], 0.5 + 1e-8)
    expect_gte(bounds[2], 0.5 - 1e-8)
    expect_equal(nrow(reps$lower), ceiling(length(orders) / 2) +
      as.integer(length(orders) %% 2 == 0))
    expect_equal(nrow(reps$upper), ceiling(length(orders) / 2) + 1)
    expect_error(
      functional_bounds(reps, exp, derivative = "positive"),
      "ordinary derivative rule"
    )
    for (k in orders) {
      expect_equal(unname(power_moment_bounds(reps, k)),
        rep(1 / (k + 1), 2),
        tolerance = 1e-8
      )
    }
  }
})

test_that("known nonuniform principal measures are recovered independently", {
  cases <- list(
    list(orders = c(2, 3, 5), x = c(0.2, 0.8), w = c(0.3, 0.7), side = "lower"),
    list(orders = c(2, 3, 5), x = c(0, 0.4, 1), w = c(0.2, 0.5, 0.3), side = "upper"),
    list(orders = c(2, 5), x = c(0, 0.6), w = c(0.4, 0.6), side = "lower"),
    list(orders = c(2, 5), x = c(0.3, 1), w = c(0.7, 0.3), side = "upper")
  )
  for (case in cases) {
    m <- power_moments(case$x, case$w, case$orders)
    reps <- generalised_principal_representations(m, case$orders)
    expect_power_constraints(reps, m, case$orders)
    expect_equal(reps[[case$side]]$nodes, case$x, tolerance = 1e-6)
    expect_equal(reps[[case$side]]$weights, case$w, tolerance = 1e-6)
  }
})

test_that("selected-power intervals narrow for nested joint-count information", {
  x <- c(0.1, 0.4, 0.8, 0.95)
  w <- c(0.2, 0.3, 0.1, 0.4)
  intervals <- lapply(list(c(2, 3), c(2, 3, 5), c(2, 3, 5, 7)), function(k) {
    reps <- generalised_principal_representations(power_moments(x, w, k), k)
    power_moment_bounds(reps, 1)
  })
  for (i in 1:2) {
    expect_gte(intervals[[i + 1]][1], intervals[[i]][1] - 1e-8)
    expect_lte(intervals[[i + 1]][2], intervals[[i]][2] + 1e-8)
  }
})

local({
  # Run the unchanged original test files through the new public entry point.
  # Thus closed forms, extremes and boundary cases cannot drift out of sync.
  env <- new.env(parent = environment())
  env$principal_representations <- function(moments, support = c(0, 1), tol = 1e-12) {
    if (!identical(support, c(0, 1))) {
      return(principalMoments::principal_representations(moments, support, tol))
    }
    result <- generalised_principal_representations(moments, seq_along(moments), tol)
    # The original functional_bounds assertions use the ordinary derivative
    # rule, which is valid here precisely because these orders are consecutive.
    class(result) <- "principal_representations"
    result
  }
  for (file in c(
    "test-evaluation.R", "test-principal-representations.R",
    "test-prevalence-closed-forms.R"
  )) {
    sys.source(file, envir = env)
  }
})

test_that("equal pool sizes recover the original prevalence intervals", {
  for (s in c(2, 5, 25, 100)) {
    for (v in 1:3) {
      # Uniform Q=1-X^s: E[(1-Q)^k] = 1/(k+1).
      reps <- generalised_principal_representations(
        1 / (seq_len(v) + 1),
        s * seq_len(v)
      )
      old <- principal_representations(1 / (seq_len(v) + 1))
      expected <- functional_bounds(old, function(q, s) -expm1(log1p(-q) / s),
        s = s,
        derivative = "positive"
      )
      b <- power_moment_bounds(reps, 1)
      expect_equal(c(1 - unname(b[2]), 1 - unname(b[1])), unname(expected),
        tolerance = 1e-8
      )
    }
  }
})

test_that("generalised input validation and failure are explicit", {
  for (orders in list(c(2, 2), c(3, 2), c(0, 2), c(-1, 2), c(NA, 2))) {
    expect_error(generalised_principal_representations(c(0.5, 0.3), orders))
  }
  expect_error(
    generalised_principal_representations(c(0.2, 0.3), c(2, 5)),
    "Infeasible"
  )
  expect_error(generalised_principal_representations(c(0.5, 0.01), c(2, 3),
    maxit = 2
  ), "did not converge")
  expect_error(generalised_principal_representations(0.5, 2, tol = 0))
  expect_error(power_moment_bounds(NULL, 1))
})


test_that("nonconsecutive boundary constraints recover unique measures", {
  for (orders in list(c(2, 3, 5), c(2, 3, 4, 5))) {
    for (case in list(
      list(x = 0, w = 1), list(x = 1, w = 1),
      list(x = 0.4, w = 1),
      list(x = c(0, 1), w = c(0.7, 0.3))
    )) {
      m <- power_moments(case$x, case$w, orders)
      reps <- generalised_principal_representations(m, orders)
      expect_power_constraints(reps, m, orders)
      for (k in c(0.5, 1, 1.5, 6)) {
        expect_equal(unname(power_moment_bounds(reps, k)),
          rep(sum(case$w * case$x^k), 2),
          tolerance = 1e-5
        )
      }
    }
  }
})

test_that("equal-pool reduction preserves extreme two-moment prevalence bounds", {
  for (m1 in c(1e-4, 0.01, 0.5, 0.99, 1 - 1e-4)) {
    for (fraction in c(0, 0.05, 0.5, 0.95, 1)) {
      m2 <- m1^2 + fraction * m1 * (1 - m1)
      old <- principal_representations(c(m1, m2))
      # Generate the equivalent moments without subtractive cancellation.
      negative_moments <- c(1 - m1, (1 - m1)^2 + fraction * m1 * (1 - m1))
      if (fraction == 1) negative_moments <- rep(1 - m1, 2)
      for (s in c(2, 5, 25, 100)) {
        reps <- generalised_principal_representations(negative_moments, s * 1:2)
        b <- power_moment_bounds(reps, 1)
        expected <- functional_bounds(old, function(q) -expm1(log1p(-q) / s),
          derivative = "positive"
        )
        expect_equal(c(1 - unname(b[2]), 1 - unname(b[1])), unname(expected),
          tolerance = 1e-6
        )
      }
    }
  }
})

test_that("mixed-size example has the independent golden-ratio upper rule", {
  reps <- generalised_principal_representations(1 / c(3, 4, 6), c(2, 3, 5))
  c_node <- (sqrt(5) - 1) / 2
  w_c <- (1 / 12) / (c_node^2 * (1 - c_node))
  w_right <- 1 / 3 - w_c * c_node^2
  expect_equal(reps$upper$nodes, c(0, c_node, 1), tolerance = 1e-8)
  expect_equal(reps$upper$weights, c(1 - w_c - w_right, w_c, w_right),
    tolerance = 1e-8
  )
  # Targets on either side of, and between, supplied exponents test orientation.
  for (t in c(0.5, 1, 2.5, 4, 8)) {
    b <- power_moment_bounds(reps, t)
    expect_lte(b[1], 1 / (t + 1) + 1e-8)
    expect_gte(b[2], 1 / (t + 1) - 1e-8)
  }
})
