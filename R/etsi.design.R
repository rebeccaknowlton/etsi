etsi.design <-
function(Study.A, n.b0 = NULL, n.b1 = NULL, psi = NULL, w.range = NULL, kappa = NULL,
                        desired.power = NULL, iterations = 100) {
  
  # Check if Study.A is a dataframe
  if (!is.data.frame(Study.A)) {
    stop("Study.A must be a dataframe.")
  }
  
  # Validate required columns
  required.cols <- c("A", "Y", "S", "W")
  if (!all(required.cols %in% names(Study.A))) {
    stop("Study.A must contain columns A, Y, S, and W.")
  }
  
  # Check that kappa provided or Study.A has delta column
  if (!("delta" %in% names(Study.A)) & is.null(kappa)) {
    stop("Please provide either a 'delta' column in Study.A or a 'kappa' value for determining strong surrogacy.")
  }
  
  # Check that user provided either 'desired.power' or both 'n.b0' and 'n.b1'
  if (is.null(desired.power) & (is.null(n.b0) | is.null(n.b1))) {
    stop("Please specify either 'desired.power' or both 'n.b0' and 'n.b1'.")
  } else if (!is.null(desired.power) & !is.null(n.b0) & !is.null(n.b1)) {
    stop("Please specify only one of 'desired.power' or both 'n.b0' and 'n.b1,' not both.")
  }
  
  # if w.range provided, subset Study.A to match
  if(!is.null(w.range)) {Study.A = Study.A[Study.A$W >= w.range[1] & Study.A$W <= w.range[2] ,]}
  
  # If delta column not present in Study.A, create it based on kappa
  if (!("delta" %in% colnames(Study.A))) {
    PTE.results <- hetsurr::hetsurr.fun(y1 = Study.A$Y[Study.A$A==1],
                                        y0 = Study.A$Y[Study.A$A==0], 
                                        s1 = Study.A$S[Study.A$A==1], 
                                        s0 = Study.A$S[Study.A$A==0], 
                                        w1 = Study.A$W[Study.A$A==1], 
                                        w0 = Study.A$W[Study.A$A==0])
    Study.A$delta <- check.strong.surr(Study.A, kappa, PTE.results)
  }
  
  # Generalized Cross-Validation
  holdout <- 0.5
  results <- rep(NA, iterations)
  for (i in 1:iterations) {
    
    # Split Study.A into train and test
    train.index <- sample(1:nrow(Study.A), round(nrow(Study.A) * holdout))
    Study.A.train <- Study.A[train.index, ]
    Study.A.test <- Study.A[-train.index, ]
    
    # Estimate pi from Study.A
    pi.b <- sum(Study.A.test$delta) / length(Study.A.test$delta)
    
    # Partition test set based on surrogate strength
    y1.weak <- Study.A.test$Y[(Study.A.test$A == 1) & (Study.A.test$delta == 0)]
    y0.weak <- Study.A.test$Y[(Study.A.test$A == 0) & (Study.A.test$delta == 0)]
    s1.strong <- Study.A.test$S[(Study.A.test$A == 1) & (Study.A.test$delta == 1)]
    s0.strong <- Study.A.test$S[(Study.A.test$A == 0) & (Study.A.test$delta == 1)]
    
    # Predict tilde Y for strong surrogates
    y1.strong <- unlist(lapply(s1.strong, get.mu.hat.0, 
                               A.s0 = Study.A.train$S[(Study.A.train$A == 0) & (Study.A.train$delta == 1)], 
                               A.y0 = Study.A.train$Y[(Study.A.train$A == 0) & (Study.A.train$delta == 1)]))
    y0.strong <- unlist(lapply(s0.strong, get.mu.hat.0, 
                               A.s0 = Study.A.train$S[(Study.A.train$A == 0) & (Study.A.train$delta == 1)], 
                               A.y0 = Study.A.train$Y[(Study.A.train$A == 0) & (Study.A.train$delta == 1)]))
  
    #Extrapolate if NA values 
    y1.strong <- extrapolate.na(y1.strong, s1.strong)
    y0.strong <- extrapolate.na(y0.strong, s0.strong)
    
    # Calculate delta.weak and delta.strong
    delta.all <- mean(Study.A.test$Y[Study.A.test$A==1]) - mean(Study.A.test$Y[Study.A.test$A==0])
    if ((is.null(y1.weak) | is.null(y0.weak))) {delta.weak <- 0} else {
      delta.weak <- mean(y1.weak) - mean(y0.weak)
      if (!is.null(psi)) {
        delta.weak <- psi * (delta.weak / delta.all)
      }
    }
    if ((is.null(y1.strong)) | (is.null(y0.strong))) {delta.strong <- 0} else {
      delta.strong <- mean(y1.strong) - mean(y0.strong)
      if (!is.null(psi)) {
        delta.strong <- psi * (delta.strong / delta.all)
      }
    }
    
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
    
    # Make corrections for n = 1 groups so empirical vars won't be NA
    if(is.na(s1.2)) {s1.2 <- 0}
    if(is.na(s2.2)) {s2.2 <- 0}
    if(is.na(s3.2)) {s3.2 <- 0}
    if(is.na(s4.2)) {s4.2 <- 0}
    
    # Calculate desired output
    if (is.null(desired.power)) {
      power <-  1 - pnorm(1.96 - ((1-pi.b)*delta.weak + pi.b*delta.strong) / (sqrt((1/n.b1) * ((1-pi.b)*s1.2 + pi.b*s2.2 + pi.b*(1-pi.b)*(mu1w - mu1s)^2) + (1/n.b0)*((1-pi.b)*s3.2 + pi.b*s4.2 + pi.b*(1-pi.b)*(mu0w-mu0s)^2))))
      results[i] <- power
    } else {
      sample.size <- ((1.96 - qnorm(1-desired.power)) * sqrt((1-pi.b)*s1.2 + pi.b*s2.2 + pi.b*(1-pi.b)*(mu1w-mu1s)^2 + (1-pi.b)*s3.2 + pi.b*s4.2 + pi.b*(1-pi.b)*(mu0w-mu0s)^2) / ((1-pi.b)*delta.weak + pi.b*delta.strong))^2
      results[i] <- sample.size
    }
  }
  if (is.null(desired.power)) {
    ans <- list("power" = mean(results))
  } else {
    ans <- list("sample.size" = mean(results))
  }
  return(ans)
}
