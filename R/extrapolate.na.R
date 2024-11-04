extrapolate.na <-
function(y.strong, s.strong) {
  if (sum(is.na(y.strong)) != 0) {
    ind = which(is.na(y.strong))
    sub.nona = cbind(s.strong, y.strong)[-ind,]
    for (i in ind) {
      mat <- cbind(abs(sub.nona[,1] - s.strong[i]), sub.nona[,2])
      mmm <- which(mat[,1] == min(mat[,1]))[1]
      y.strong[i] <- sub.nona[mmm,2]
    }
  }
  return(y.strong)
}
