get.mu.hat.0 <-
function(s, A.s0, A.y0) {
  kernel <- function(x, h) {return(dnorm(x / h) / h)}
  h.0 <- bw.nrd(A.s0)*length(A.s0)^(-0.2)
  return(sum(kernel(A.s0 - s, h.0) * A.y0 / sum(kernel(A.s0 - s, h.0))))
}
