# Follow-up: tolerance sensitivity for representative mixed cases.
# run.R defines run_case and the existing fixture objects; source only its
# definitions above the first execution loop to avoid rerunning the main sweep.
code <- readLines("benchmarks/robustness/run.R")
end <- grep("^for \\(i in seq_len\\(nrow\\(cases\\)\\)\\)", code)[1]
eval(parse(text = code[seq_len(end - 1)]), .GlobalEnv)
rows <- lapply(seq_len(nrow(read.csv(file.path(root,"runs.csv")))), function(i) {
  read.csv(file.path(root,"runs.csv"), stringsAsFactors = FALSE)[i, ]
})
rules <- readRDS(file.path(root,"rules.rds"))
which_cases <- which(cases$family == "mixed_known" &
                      cases$label %in% c("[10, 40]x[1, 1]_lower",
                         "[10, 40]x[2, 2]_lower", "[25, 26]x[1, 1]_lower",
                         "[2, 101]x[1, 1]_lower", "[100, 101]x[1, 1]_lower"))
for (i in which_cases) {
  for (tol in c(1e-12, 1e-15, 1e-20, 1e-30)) {
    run_case(cases[i, ], tol, paste0("mixed_tol_", tol))
  }
}
# Low order on a general interval: cancellation occurs before moment inversion.
# Independent raw moments for Uniform(a,a+1) are evaluated here using the
# binomial identity with only positive summands, avoiding subtraction.
interval_rows <- list()
for (a in c(0, 10, 1000, 1e6)) {
  for (v in c(2,4,8,12)) {
    m <- vapply(1:v, function(k) {
      j <- 0:k
      sum(choose(k,j) * a^(k-j)/(j+1))
    }, numeric(1))
    warnings <- character()
    answer <- tryCatch(withCallingHandlers(principal_representations(m, c(a,a+1)),
                          warning=function(w) {warnings <<- c(warnings,conditionMessage(w));
                                               invokeRestart("muffleWarning")}),
                       error=function(e) e)
    interval_rows[[length(interval_rows)+1]] <- data.frame(a=a,v=v,
      status=if(inherits(answer,"error")) conditionMessage(answer) else "returned",
      boundary=if(inherits(answer,"error")) NA else answer$boundary)
  }
}
write.csv(do.call(rbind, interval_rows),file.path(root,"intervals.csv"),row.names=FALSE)
