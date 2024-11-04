check.strong.surr <-
function(df, k, PTE.results) {
  closest.index <- sapply(df$W, function(w) {which.min(abs(PTE.results$w.values - w))})
  return(1*(PTE.results$R.w.s[closest.index] > k))
}
