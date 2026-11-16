# From the repository root, using Rmpfr. These are timing/accuracy benchmarks,
# not guarantees for all moment vectors. Use an external process timeout for
# expensive mixed problems (for example GNU timeout on Linux).
# Rscript benchmarks/robustness/precision-check.R
# timeout 240 Rscript benchmarks/robustness/precision-check.R --mixed23
for (file in list.files("R",full.names=TRUE)) sys.source(file,.GlobalEnv)
if ("--mixed23" %in% commandArgs(TRUE)) {
  n <- 5L
  orders <- sort(unique(as.vector(outer(0:n,0:n,function(i,j)2*i+3*j))))[-1]
  d <- length(orders)
  r <- (d+1L)/2L
  nodes <- Rmpfr::mpfr(seq_len(r),512)/(r+1L)
  weights <- Rmpfr::mpfr(seq_len(r),512)/sum(seq_len(r))
  m <- methods::as(mp_basis(mp_value(orders,512),nodes)%*%weights,"mpfr")
  print(system.time(result <- generalised_principal_representations(
    m,orders,max_precision=768
  )))
  print(result$diagnostics)
  print(c(node_error=mp_num(max(abs(result$high_precision$lower$nodes-nodes))),
          weight_error=mp_num(max(abs(result$high_precision$lower$weights-weights)))))
} else {
  for(v in c(17,50,100)) {
    cat("Uniform moments:",v,"\n")
    m <- 1/Rmpfr::mpfr(2:(v+1L),512)
    print(system.time(result <- principal_representations(m)))
    print(result$diagnostics)
    print(functional_bounds(result,function(q)-expm1(log1p(-q)/25),
                            derivative="positive"))
  }
}
