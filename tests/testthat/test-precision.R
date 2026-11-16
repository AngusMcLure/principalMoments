precision_fixture <- function(id) {
  x <- read.csv(test_path("fixtures", "precision-inputs.csv"), colClasses = "character")
  r <- read.csv(test_path("fixtures", "precision-reference.csv"), colClasses = "character")
  list(input = x[x$id == id, ], reference = r[r$id == id, ])
}

expect_reference <- function(result, reference, fun, tolerance = 1e-9) {
  values <- evaluate_representations(result, fun)
  for (i in seq_len(nrow(reference))) {
    ref <- reference[i, ]
    expect_equal(unname(values[ref$side]), as.numeric(ref$value), tolerance = tolerance)
    nodes <- as.numeric(strsplit(ref$nodes, ";", fixed = TRUE)[[1]])
    weights <- as.numeric(strsplit(ref$weights, ";", fixed = TRUE)[[1]])
    expect_equal(result[[ref$side]]$nodes, nodes, tolerance = tolerance)
    expect_equal(result[[ref$side]]$weights, weights, tolerance = tolerance)
  }
}

test_that("long double inputs fail safely instead of claiming point identification", {
  expect_error(principal_representations(1 / (2:18)), "Unresolved moment rank")
  expect_error(principal_representations(1 / (2:31)), "Unresolved moment rank")
  expect_error(
    generalised_principal_representations(1 / (2:18), 25 * (1:17)),
    "Unresolved moment rank"
  )
  expect_error(
    principal_representations(1 / (2:18), precision = 192),
    "Unresolved input precision"
  )
})

test_that("MPFR recurrence recovers the previously collapsed prevalence intervals", {
  for (id in c("17", "133", "161")) {
    fixture <- precision_fixture(id)
    result <- principal_representations(fixture$input$moment)
    s <- if (id == "161") 50 else 25
    expect_false(result$boundary)
    expect_identical(result$diagnostics$status, "precision_checked")
    expect_gt(result$diagnostics$working_precision, 100)
    expect_reference(result, fixture$reference, function(q) -expm1(log1p(-q) / s))
    b <- functional_bounds(result, function(q) -expm1(log1p(-q) / s), derivative = "positive")
    expect_lt(b[1], b[2])
    if (id != "17") {
      truth <- if (id == "161") 0.001 else 0.01
      expect_lte(b[1], truth)
      expect_gte(b[2], truth)
    }
  }
})

test_that("fifty moment MPFR nodes agree with independent Jacobi quadrature", {
  fixture <- precision_fixture("uniform50")
  result <- principal_representations(fixture$input$moment)
  expect_false(result$boundary)
  expect_equal(nrow(result$lower), 26L)
  expect_reference(result, fixture$reference, function(q) -expm1(log1p(-q) / 25), 1e-10)
  expect_lt(max(result$diagnostics$refinement_change), 1e-10)
})

test_that("continuation repairs mixed-size failures and relative stopping errors", {
  for (id in c("198", "206")) {
    fixture <- precision_fixture(id)
    result <- generalised_principal_representations(
      fixture$input$moment,
      as.numeric(fixture$input$order)
    )
    expect_reference(result, fixture$reference, function(x) 1 - x, 1e-8)
    expect_false(result$boundary)
  }
})

test_that("precise evaluations retain endpoint information beyond double precision", {
  # Mapping a node of X^s back to X rounds the display node to one, while
  # its MPFR value still carries a nonzero prevalence.
  result <- generalised_principal_representations("0.5", 1e20)
  expect_equal(result$lower$nodes, 1)
  value <- evaluate_representations(result, function(x) 1 - x)
  expect_gt(value[1], 0)
  expect_equal(unname(value[1]), log(2) / 1e20, tolerance = 1e-10)
})

test_that("controls and unresolved boundaries cannot masquerade as sharp bounds", {
  expect_error(principal_representations(0.5, accuracy = 0), "accuracy")
  expect_error(principal_representations("0.5", max_precision = 64), "max_precision")
  expect_error(principal_representations("0.5", precision = Inf), "precision")
  expect_error(principal_representations(c("0.5", "nonsense")), "moments")
  expect_error(principal_representations("0.5", support = c(1, 2)), "support")
  reps <- principal_representations(0.3)
  reps$diagnostics$status <- "numerical_boundary"
  expect_error(functional_bounds(reps, exp, derivative = "positive"), "Unresolved boundary")
})

if (identical(Sys.getenv("PRINCIPAL_MOMENTS_STRESS"), "true")) {
  test_that("one hundred moments recover independent principal rules", {
    fixture <- precision_fixture("uniform100")
    result <- principal_representations(fixture$input$moment)
    expect_reference(result, fixture$reference, function(q) -expm1(log1p(-q) / 25), 1e-10)
  })
  test_that("high-power mixed fixtures require preserved input information", {
    for (id in c("214", "216")) {
      fixture <- precision_fixture(id)
      result <- generalised_principal_representations(
        fixture$input$moment,
        as.numeric(fixture$input$order),
        maxit = 300L
      )
      expect_reference(result, fixture$reference, function(x) 1 - x, 1e-8)
    }
  })
}

test_that("functional stability is checked separately from node stability", {
  reps <- generalised_principal_representations("0.5", 1e20)
  reps$precision_reference$lower$nodes <- 1 - Rmpfr::mpfr("1e-10", 256)
  expect_error(
    evaluate_representations(reps, function(x) -expm1(log1p(-x) / 100)),
    "Functional evaluation is unresolved"
  )
})

test_that("power transformations automatically preserve otherwise lost nodes", {
  reps <- generalised_principal_representations(0.5, 1e20)
  expect_false(is.null(reps$high_precision))
  value <- evaluate_representations(reps, function(x) 1 - x)
  expect_gt(value[1], 0)
  expect_equal(unname(value[1]), log(2) / 1e20, tolerance = 1e-10)
})
