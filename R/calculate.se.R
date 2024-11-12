calculate.se <-
function(y1.weak, y0.weak, y1.strong, y0.strong) {
  n.b.treat <- length(y1.weak) + length(y1.strong)
  n.b.control <- length(y0.weak) + length(y0.strong)
  
  n1w <- length(y1.weak)
  n1s <- length(y1.strong)
  n0w <- length(y0.weak)
  n0s <- length(y0.strong)
  
  # Set empirical means and variances as 0, then fill in for nonempty groups
  mu1w <- 0; mu1s <- 0; mu0w <- 0; mu0s <- 0
  s1.2 <- 0; s2.2 <- 0; s3.2 <- 0; s4.2 <- 0
  if(n1w > 0) {mu1w <- mean(y1.weak); s1.2 <- var(y1.weak)}
  if(n1s > 0) {mu1s <- mean(y1.strong); s2.2 <- var(y1.strong)}
  if(n0w > 0) {mu0w <- mean(y0.weak); s3.2 <- var(y0.weak)}
  if(n0s > 0) {mu0s <- mean(y0.strong); s4.2 <- var(y0.strong)}
  
  # Make corrections for n = 0 or n = 1 groups so se.delta.P won't return NA
  if(n.b.treat == 0) {n.b.treat <- 1}
  if(n.b.control == 0) {n.b.control <- 1}
  if(is.na(s1.2)) {s1.2 <- 0}
  if(is.na(s2.2)) {s2.2 <- 0}
  if(is.na(s3.2)) {s3.2 <- 0}
  if(is.na(s4.2)) {s4.2 <- 0}
  
  var1 <- (1 / n.b.treat) * ((n1w / n.b.treat) * s1.2 + (n1s / n.b.treat) * s2.2 +
                                n1w * n1s / (n.b.treat^2) * (mu1w - mu1s)^2)
  
  var0  = (1 / n.b.control) * ((n0w / n.b.control) * s3.2 + (n0s / n.b.control) * s4.2 +
                                 n0w * n0s /(n.b.control^2) * (mu0w-mu0s)^2)
  
  se.delta.P <- sqrt(var1 +var0)
  return(se.delta.P)
}
