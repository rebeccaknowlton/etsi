etsi.main <-
function(Study.A, Study.B, kappa = NULL) {
  
  # Check if both Study.A and Study.B are dataframes
  if (!is.data.frame(Study.A) | !is.data.frame(Study.B)) {
    stop("Both Study.A and Study.B must be dataframes.")
  }
  
  # Validate required columns
  required.cols <- c("A", "Y", "S", "W")
  if (!all(required.cols %in% names(Study.A)) | !all(required.cols %in% names(Study.B))) {
    stop("Both dataframes must contain columns A, Y, S, and W.")
  }
  
  # Check if the delta column is present in Study.A and Study.B. If not, create it based on kappa.
  if (!("delta" %in% colnames(Study.A)) | !("delta" %in% colnames(Study.B))) {
    if (is.null(kappa)) {
      stop("Either the delta column must be present in both dataframes, or a kappa threshold must be provided.")
    } else {
      PTE.results <- hetsurr::hetsurr.fun(y1 = Study.A$Y[Study.A$A==1],
                                          y0 = Study.A$Y[Study.A$A==0], 
                                          s1 = Study.A$S[Study.A$A==1], 
                                          s0 = Study.A$S[Study.A$A==0], 
                                          w1 = Study.A$W[Study.A$A==1], 
                                          w0 = Study.A$W[Study.A$A==0])
      Study.A$delta <- check.strong.surr(Study.A, kappa, PTE.results)
      Study.B$delta <- check.strong.surr(Study.B, kappa, PTE.results)
    }
  }
  
  y1.weak <- Study.B$Y[(Study.B$A == 1) & (Study.B$delta == 0)]
  y0.weak <- Study.B$Y[(Study.B$A == 0) & (Study.B$delta == 0)]
  s1.strong <- Study.B$S[(Study.B$A == 1) & (Study.B$delta == 1)]
  s0.strong <- Study.B$S[(Study.B$A == 0) & (Study.B$delta == 1)]
  A.s0.strong <- Study.A$S[(Study.A$A == 0) & (Study.A$delta == 1)]
  A.y0.strong <- Study.A$Y[(Study.A$A == 0) & (Study.A$delta == 1)]
  
  # For strong surrogate, predict tilde Y based on Study A control group
  y1.strong <- unlist(lapply(s1.strong, get.mu.hat.0, A.s0 = A.s0.strong, A.y0 = A.y0.strong))
  y0.strong <- unlist(lapply(s0.strong, get.mu.hat.0, A.s0 = A.s0.strong, A.y0 = A.y0.strong))
  
  #Extrapolate if NA values 
  y1.strong <- extrapolate.na(y1.strong, s1.strong)
  y0.strong <- extrapolate.na(y0.strong, s0.strong)
  
  # Calculate and return estimates
  delta.P <- mean(c(y1.weak, y1.strong)) - mean(c(y0.weak, y0.strong))
  se.delta.P <- calculate.se(y1.weak = y1.weak, y1.strong = y1.strong, y0.weak = y0.weak, y0.strong = y0.strong) 
  p.value <- 2 * (1 - pnorm(abs(delta.P / se.delta.P)))
  return.list <- list("delta.P" = delta.P,
                      "se.delta.P" = se.delta.P,
                      "p.value" = p.value)
  
  return(return.list)  
}
