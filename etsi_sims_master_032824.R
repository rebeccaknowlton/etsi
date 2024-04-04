#setting <- 1
n <- 1000 # sample size for each group

num.sim <- 100

get.parameters <- function(setting){
  if(setting ==1){
    w.lower.A <- 2
    w.upper.A <- 6
    w.lower.B <- 4
    w.upper.B <- 6
    s0.scale <- 0.5  
    s0.shape <- 0.5
    s1.scale <- 1
    s1.shape <- 1
    beta0 <- 5
    beta1 <- 0.1
    beta2 <- 2.6
    beta3 <- 0.1
    beta4 <- 0.1
    beta5 <- 0.6
    sd.y <- 6
    return(list("w.lower.A" = w.lower.A, "w.upper.A" = w.upper.A, "w.lower.B" = w.lower.B, "w.upper.B" = w.upper.B, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
  if(setting ==2){
    w.lower.A <- 2
    w.upper.A <- 6
    w.lower.B <- 4
    w.upper.B <- 6
    s0.scale <- 0.5  
    s0.shape <- 0.5
    s1.scale <- 1
    s1.shape <- 1
    beta0 <- 5
    beta1 <- 0.1
    beta2 <- 1
    beta3 <- 0.1
    beta4 <- 0.1
    beta5 <- 0.1
    sd.y <- 6
    return(list("w.lower.A" = w.lower.A, "w.upper.A" = w.upper.A, "w.lower.B" = w.lower.B, "w.upper.B" = w.upper.B, "s0.scale" = s0.scale, "s0.shape" = s0.shape, "s1.scale" = s1.scale, "s1.shape" = s1.shape, "beta0" = beta0, "beta1" = beta1, "beta2" = beta2, "beta3" = beta3, "beta4" = beta4, "beta5" = beta5, "sd.y" = sd.y))
  }
}

# plot of R vs W to see what level of heterogeneity is present. Need to mess with numbers to get wider range but smaller treatment effect or power will just = 1
#w.grid <- runif(1000, 2, 6)
#plot(w.grid, get.truth(2, w.grid)$R)
#range(get.truth(2,w.grid)$R)[2] - range(get.truth(2,w.grid)$R)[1]
#tmp.data = gen.data(1000,2)
#est.delta.B(tmp.data$Y[tmp.data$A==1],tmp.data$Y[tmp.data$A==0])

gen.data <- function(n, setting, study="A"){
  params <- get.parameters(setting=setting)
  if(setting == 1) {
    if (study == "A") {
      data.temp <- data.frame(A = c(rep(1,n),rep(0,n)),
                              W = runif(n*2, params$w.lower.A, params$w.upper.A))  
    } else if (study == "B") {
      data.temp <- data.frame(A = c(rep(1,n),rep(0,n)),
                              W = runif(n*2, params$w.lower.B, params$w.upper.B))
    }

    data.temp$S[data.temp$A==1] <- rgamma(n, shape = params$s1.shape, scale = params$s1.scale)
    data.temp$S[data.temp$A==0] <- rgamma(n, shape = params$s0.shape, scale = params$s0.scale)
    data.temp$Y <- params$beta0 + params$beta1 * data.temp$A + params$beta2 * data.temp$S + params$beta3 * data.temp$A * data.temp$S +
      data.temp$W * params$beta4  +  data.temp$W *data.temp$A * params$beta5 + rnorm(n*2, 0, params$sd.y)
    return(data.temp)
  }
  if(setting == 2) {
    if (study == "A") {
      data.temp <- data.frame(A = c(rep(1,n),rep(0,n)),
                              W = runif(n*2, params$w.lower.A, params$w.upper.A))  
    } else if (study == "B") {
      data.temp <- data.frame(A = c(rep(1,n),rep(0,n)),
                              W = runif(n*2, params$w.lower.B, params$w.upper.B))
    }
    
    data.temp$S[data.temp$A==1] <- rgamma(n, shape = params$s1.shape, scale = params$s1.scale)
    data.temp$S[data.temp$A==0] <- rgamma(n, shape = params$s0.shape, scale = params$s0.scale)
    data.temp$Y <- params$beta0 + params$beta1 * data.temp$A + params$beta2 * data.temp$S + params$beta3 * data.temp$A * data.temp$S +
      data.temp$W * params$beta4  +  (data.temp$W^2) *data.temp$A * params$beta5 + rnorm(n*2, 0, params$sd.y)
    return(data.temp)
  }
}

get.truth <- function(setting, grid){
  if(setting ==1){
    params <- get.parameters(setting=1)
    delta.s <- params$beta1+params$beta3*params$s0.shape*params$s0.scale + grid*params$beta5
    delta <- params$beta1+(params$beta2+params$beta3)*(params$s1.shape*params$s1.scale) - (params$beta2)*(params$s0.shape*params$s0.scale) + grid*params$beta5
    R = 1-delta.s/delta
    return(list("delta.s" = delta.s, "delta" = delta, "R" = R))
  }
  if(setting ==2){
    params <- get.parameters(setting=2)
    delta.s <- params$beta1+params$beta3*params$s0.shape*params$s0.scale + (grid^2)*params$beta5
    delta <- params$beta1+(params$beta2+params$beta3)*(params$s1.shape*params$s1.scale) - (params$beta2)*(params$s0.shape*params$s0.scale) + (grid^2)*params$beta5
    R = 1-delta.s/delta
    return(list("delta.s" = delta.s, "delta" = delta, "R" = R))
  }
}

# generate study A

set.seed(1) # should be the same for all parallel runs
study.A <- gen.data(n, setting, "A")
control.A <- study.A[study.A$A == 0,]
treat.A <- study.A[study.A$A == 1,]

# kernel smoothing
kernel <- function(x, h) {return(dnorm(x / h) / h)}
h.0 <- bw.nrd(control.A$S)*length(control.A$S)^(-0.2)
get.mu.hat.0 <- function(s) {return(sum(kernel(control.A$S - s, h.0) * control.A$Y / sum(kernel(control.A$S - s, h.0))))}

get.mu.hat.0.boot <- function(control.A, s) {
  control.A.boot <- control.A[sample(1:n, n, replace = TRUE),]
  h.0 <- bw.nrd(control.A.boot$S)*length(control.A.boot$S)^(-0.2)
  return(sum(kernel(control.A.boot$S - s, h.0) * control.A.boot$Y / sum(kernel(control.A.boot$S - s, h.0))))
}

est.delta.B <- function(y1, y0) {
  delta.B <- mean(y1) - mean(y0)
  se.delta.B <- sqrt((1 / n) * var(y1) + (1 / n) * var(y0))    
  return(list("delta.B" = delta.B, "se.delta.B" = se.delta.B))
}

est.delta.B.boot <- function(s1, s0) {
  # predict y values based on control group from study A
  y1 <- unlist(lapply(s1, get.mu.hat.0))
  y0 <- unlist(lapply(s0, get.mu.hat.0))
  if (sum(is.na(y1)) != 0) {
    ind = which(is.na(y1))
    sub.nona = cbind(s1, y1)[-ind,]
    for (i in ind) {
      mat <- cbind(abs(sub.nona[,1] - s1[i]), sub.nona[,2])
      mmm <- which(mat[,1] == min(mat[,1]))[1]
      y1[i] <- sub.nona[mmm,2]
    }
  }
  if (sum(is.na(y0)) != 0) {
    ind = which(is.na(y0))
    sub.nona = cbind(s0, y0)[-ind,]
    for (i in ind) {
      mat <- cbind(abs(sub.nona[,1] - s0[i]), sub.nona[,2])
      mmm <- which(mat[,1] == min(mat[,1]))[1]
      y0[i] <- sub.nona[mmm,2]
    }
  }
  
  delta.B <- mean(y1) - mean(y0)
  num.boot <- 200
  boot.delta.B <- rep(NA, num.boot)
  y1.tmp <- rep(NA, n)
  y0.tmp <- rep(NA, n)
  for (j in 1:num.boot) {
    # need to resample a
    boot.s1 <- s1[sample(1:n, n, replace = TRUE)]
    boot.s0 <- s0[sample(1:n, n, replace = TRUE)]
    boot.control.A <- control.A[sample(1:n, n, replace = TRUE),]
    for (idx in 1:n) {
      y1.tmp[idx] <- get.mu.hat.0.boot(boot.control.A, boot.s1[idx])
      y0.tmp[idx] <- get.mu.hat.0.boot(boot.control.A, boot.s0[idx])
    }
    if (sum(is.na(y1.tmp)) != 0) {
      ind = which(is.na(y1.tmp))
      sub.nona = cbind(boot.s1, y1.tmp)[-ind,]
      for (i in ind) {
        mat <- cbind(abs(sub.nona[,1] - boot.s1[i]), sub.nona[,2])
        mmm <- which(mat[,1] == min(mat[,1]))[1]
        y1.tmp[i] <- sub.nona[mmm,2]
      }
    }
    if (sum(is.na(y0.tmp)) != 0) {
      ind = which(is.na(y0.tmp))
      sub.nona = cbind(boot.s0, y0.tmp)[-ind,]
      for (i in ind) {
        mat <- cbind(abs(sub.nona[,1] - boot.s0[i]), sub.nona[,2])
        mmm <- which(mat[,1] == min(mat[,1]))[1]
        y0.tmp[i] <- sub.nona[mmm,2]
      }
    }
    boot.delta.B[j] <- mean(y1.tmp) - mean(y0.tmp)
  }
  se.delta.B <- sd(boot.delta.B)
  return(list("delta.B" = delta.B, "se.delta.B" = se.delta.B))
}

# generate study B with same data generation setup, and test H_0: \delta_B = 0 using the predicted outcome mu.hat

y.obs.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                              "se.delta.B" = rep(NA, num.sim),
                              "z" = rep(NA, num.sim))
y.pred.control.boot.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                                            "se.delta.B" = rep(NA, num.sim),
                                            "z" = rep(NA, num.sim))

#THIS IS WHAT MAKES EACH PARALLEL VERSION DIFFERENT
set.seed(parallel.num*100)


for (i in 1:num.sim) {
  study.B <- gen.data(n, setting, "B")
  
  # save estimates from observed y's
  results.temp <- est.delta.B(y1 = study.B$Y[study.B$A==1], y0 = study.B$Y[study.B$A==0])
  y.obs.estimates$delta.B[i] <- results.temp$delta.B
  y.obs.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.obs.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  # save estimates from predicted y's using control group mu.hat and bootstrapped variance
  results.temp <- est.delta.B.boot(s1 = study.B$S[study.B$A == 1], s0 = study.B$S[study.B$A == 0])
  y.pred.control.boot.estimates$delta.B[i] <- results.temp$delta.B
  y.pred.control.boot.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.pred.control.boot.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  print(i)
}


# compile and save summary table 
df <- data.frame("avg.true.delta.B" = rep(get.truth(setting, c(5))$delta, 2),
                 "avg.est.delta.B" = c(mean(y.obs.estimates$delta.B), mean(y.pred.control.boot.estimates$delta.B)),
                 "prop.null.reject" = c(sum(abs(y.obs.estimates$z) > 1.96) / num.sim, sum(abs(y.pred.control.boot.estimates$z) > 1.96) / num.sim))
df <- t(df) 
colnames(df) <- c("observed.y", "predict.y.control.boot")
write.table(df, paste("etsi.output.summary",setting, "_032824",parallel.num,".txt",sep=""), quote = FALSE, row.names = TRUE)

# more detailed output table for troubleshooting
est.results <- cbind(y.obs.estimates, y.pred.control.boot.estimates)
colnames(est.results) <- c("delta.B.obs","se.delta.B.obs","z.obs","delta.B.pred.con.boot","se.delta.B.pred.con.boot","z.pred.con.boot")
write.table(est.results, paste("etsi.output.all",setting, "_032824",parallel.num,".txt",sep=""), quote = FALSE, row.names = FALSE)


