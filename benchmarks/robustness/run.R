# Run after generate.py, from the repository root. No solver code is changed.
for (file in list.files("R", full.names = TRUE)) sys.source(file, .GlobalEnv)
root <- "benchmarks/robustness/results"
cases <- read.csv(file.path(root, "cases.csv"), stringsAsFactors = FALSE)
inputs <- read.csv(file.path(root, "inputs.csv"))
references <- read.csv(file.path(root, "reference.csv"), stringsAsFactors = FALSE)
rows <- list()
rules <- list()
run_case <- function(case, tol = if (case$api == "fixed") 1e-12 else 1e-9,
                     tag = "default") {
  x <- inputs[inputs$id == case$id, ]
  warnings <- character()
  start <- proc.time()[[3]]
  setTimeLimit(cpu = 12, elapsed = 15, transient = TRUE)
  result <- tryCatch(withCallingHandlers({
    if (case$api == "fixed") principal_representations(x$moment, tol = tol) else
      generalised_principal_representations(x$moment, x$order, tol = tol)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
  elapsed <- proc.time()[[3]] - start
  row <- cbind(case, data.frame(tag = tag, tol = tol, elapsed = elapsed,
    status = if (inherits(result, "error")) "error" else "returned",
    warning = paste(unique(warnings), collapse = " | "),
    error = if (inherits(result, "error")) conditionMessage(result) else "",
    boundary = NA, nlower = NA, nupper = NA, residual = NA,
    lower_value = NA_real_, upper_value = NA_real_, reference_error = NA_real_,
    node_error = NA_real_, weight_error = NA_real_, covers_truth = NA))
  if (!inherits(result, "error")) {
    row$boundary <- result$boundary
    row$nlower <- nrow(result$lower)
    row$nupper <- nrow(result$upper)
    row$residual <- max(result$diagnostics$residuals)
    f <- if (case$api == "fixed") function(q) -expm1(log1p(-q) / case$s) else
      function(x) 1 - x
    values <- evaluate_representations(result, f)
    row$lower_value <- values[1]
    row$upper_value <- values[2]
    reference <- references[references$id == case$id, ]
    node_errors <- weight_errors <- errors <- numeric()
    for (j in seq_len(nrow(reference))) {
      ref <- reference[j, ]
      rule <- result[[ref$side]]
      errors <- c(errors, abs(values[ref$side] - as.numeric(ref$value)))
      nodes <- as.numeric(strsplit(ref$nodes, ";", fixed = TRUE)[[1]])
      weights <- as.numeric(strsplit(ref$weights, ";", fixed = TRUE)[[1]])
      if (nrow(rule) == length(nodes)) {
        node_errors <- c(node_errors, max(abs(rule$nodes - nodes)))
        weight_errors <- c(weight_errors, max(abs(rule$weights - weights)))
      } else {
        node_errors <- c(node_errors, Inf)
        weight_errors <- c(weight_errors, Inf)
      }
    }
    if (length(errors)) row$reference_error <- max(errors)
    if (length(node_errors)) row$node_error <- max(node_errors)
    if (length(weight_errors)) row$weight_error <- max(weight_errors)
    if (case$family == "mixed_beta") {
      interval <- if (case$d %% 2) values else rev(values)
      row$covers_truth <- interval[1] <= 0.01 + 1e-8 && interval[2] >= 0.01 - 1e-8
    }
    rules[[paste(case$id, tag, sep = "_")]] <<- result
  }
  rows[[length(rows) + 1L]] <<- row
  write.csv(do.call(rbind, rows), file.path(root, "runs.csv"), row.names = FALSE)
  saveRDS(rules, file.path(root, "rules.rds"))
  message(case$id, ": ", case$family, " ", case$label, " d=", case$d, " ", row$status,
          " reference_error=", signif(row$reference_error, 3), " time=", round(elapsed, 2))
}
for (i in seq_len(nrow(cases))) run_case(cases[i, ])
# Tolerance is not a substitute for precision: test both lower and higher
# numerical-rank thresholds on identical inputs, without changing the code.
for (i in which(cases$family == "fixed" & cases$label == "Q_uniform" &
                cases$d %in% c(10,12,14,16,18,20))) {
  for (tol in c(1e-9,1e-15,1e-18)) run_case(cases[i, ], tol, paste0("tol_",tol))
}
writeLines(capture.output(sessionInfo()), file.path(root, "session.txt"))
