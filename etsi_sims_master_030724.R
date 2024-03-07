### SET NUMBERS HERE - only touch this! ###

setting <- 1
n <- 1000 # sample size for each group

############################################

num.sim <- 100
# Set following parameters for simulation settings:
# S = s.intercept + s.treat * G + N(0, sd.s)
# Y = beta0 + beta1 * S + beta2 * G + N(0, sd.y)

get.parameters <- function(setting){
  if(setting ==1){
    # surrogate is perfect, R = 1, delta = 1
    s.intercept <- 3
    s.treat <- 1
    sd.s <- 3
    beta0 <- 5
    beta1 <- 1
    beta2 <- 0
    sd.y <- 6
    return(list("s.intercept" = s.intercept, "s.treat" = s.treat, "sd.s" = sd.s, "beta0"=beta0, "beta1"=beta1, "beta2" = beta2, "sd.y"=sd.y))
  }
  if(setting == 2){
    # surrogate is imperfect, R = 0.5, delta = 1
    s.intercept <- 3
    s.treat <- 1
    sd.s <- 3
    beta0 <- 5
    beta1 <- 0.5
    beta2 <- 0.5
    sd.y <- 6
    return(list("s.intercept" = s.intercept, "s.treat" = s.treat, "sd.s" = sd.s, "beta0"=beta0, "beta1"=beta1, "beta2" = beta2, "sd.y"=sd.y))
  }
  if(setting == 3){
    # not a surrogate but there is a treatment effect, R = 0, delta = 1
    s.intercept <- 3
    s.treat <- 0
    sd.s <- 3
    beta0 <- 5
    beta1 <- 0
    beta2 <- 1
    sd.y <- 6
    return(list("s.intercept" = s.intercept, "s.treat" = s.treat, "sd.s" = sd.s, "beta0"=beta0, "beta1"=beta1, "beta2" = beta2, "sd.y"=sd.y))
  }
  if(setting == 4){
    # not a surrogate, no treatment effect
    s.intercept <- 3
    s.treat <- 0
    sd.s <- 3
    beta0 <- 5
    beta1 <- 0
    beta2 <- 0
    sd.y <- 6
    return(list("s.intercept" = s.intercept, "s.treat" = s.treat, "sd.s" = sd.s, "beta0"=beta0, "beta1"=beta1, "beta2" = beta2, "sd.y"=sd.y))
  }
  if(setting == 5){
    # not is a "surrogate" as far as predicting Y, but no treatment effect
    s.intercept <- 3
    s.treat <- 0
    sd.s <- 3
    beta0 <- 5
    beta1 <- 1
    beta2 <- 0
    sd.y <- 6
    return(list("s.intercept" = s.intercept, "s.treat" = s.treat, "sd.s" = sd.s, "beta0"=beta0, "beta1"=beta1, "beta2" = beta2, "sd.y"=sd.y))
  }
}

gen.data = function(n,setting){
  params = get.parameters(setting=setting)
  data.temp <- data.frame(G = append(rep(1, n), rep(0, n)))
  data.temp$S <- params$s.intercept + params$s.treat * data.temp$G + rnorm(2*n, 0, params$sd.s)
  data.temp$Y <- params$beta0 + params$beta1 * data.temp$S + params$beta2 * data.temp$G + rnorm(2*n, 0, params$sd.y)
  return(data.temp)
}  

get.truth = function(setting){
  params = get.parameters(setting)
  delta.s = params$beta2
  delta = params$beta1 * params$s.treat + params$beta2
  R = 1-delta.s/delta
  return(list("delta.s" = delta.s, "delta" = delta, "R" = R))
}


# generate study A

set.seed(1) # should be the same for all parallel runs
study.A <- gen.data(n, setting)
control.A <- study.A[study.A$G == 0,]
treat.A <- study.A[study.A$G == 1,]

# kernel smoothing to get mu.hat(s) in the treatment and control group for study A

kernel <- function(x, h) {return(dnorm(x / h) / h)}
h.0 <- bw.nrd(control.A$S)*length(control.A$S)^(-0.2)
h.1 <- bw.nrd(treat.A$S)*length(treat.A$S)^(-0.2)

get.mu.hat.0 <- function(s) {return(sum(kernel(control.A$S - s, h.0) * control.A$Y / sum(kernel(control.A$S - s, h.0))))}
get.mu.hat.1 <- function(s) {return(sum(kernel(treat.A$S - s, h.1) * treat.A$Y / sum(kernel(treat.A$S - s, h.1))))}

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
y.pred.separate.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                                        "se.delta.B" = rep(NA, num.sim),
                                        "z" = rep(NA, num.sim))
y.pred.control.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                                       "se.delta.B" = rep(NA, num.sim),
                                       "z" = rep(NA, num.sim))
y.pred.control.boot.estimates <- data.frame("delta.B" = rep(NA, num.sim),
                                            "se.delta.B" = rep(NA, num.sim),
                                            "z" = rep(NA, num.sim))

#THIS IS WHAT MAKES EACH PARALLEL VERSION DIFFERENT
set.seed(parallel.num*100)


for (i in 1:num.sim) {
  start_time <- Sys.time()
  study.B <- gen.data(n, setting)
  
  # save estimates from observed y's
  results.temp <- est.delta.B(y1 = study.B$Y[study.B$G==1], y0 = study.B$Y[study.B$G==0])
  y.obs.estimates$delta.B[i] <- results.temp$delta.B
  y.obs.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.obs.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  # get predicted values from mu.hat functions, separately for control/treatment groups
  study.B$Y.pred.separate[study.B$G == 0] <- unlist(lapply(study.B$S[study.B$G == 0], get.mu.hat.0))
  study.B$Y.pred.separate[study.B$G == 1] <- unlist(lapply(study.B$S[study.B$G == 1], get.mu.hat.1))
  # save estimates from predicted y's using separate mu.hats
  results.temp <- est.delta.B(y1 = study.B$Y.pred.separate[study.B$G == 1], y0 = study.B$Y.pred.separate[study.B$G == 0])
  y.pred.separate.estimates$delta.B[i] <- results.temp$delta.B
  y.pred.separate.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.pred.separate.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  # get predicted values from mu.hat using control group only
  study.B$Y.pred.control <- unlist(lapply(study.B$S, get.mu.hat.0))
  # save estimates from predicted y's using control group mu.hat
  results.temp <- est.delta.B(y1 = study.B$Y.pred.control[study.B$G == 1], y0 = study.B$Y.pred.control[study.B$G == 0])
  y.pred.control.estimates$delta.B[i] <- results.temp$delta.B
  y.pred.control.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.pred.control.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  # save estimates from predicted y's using control group mu.hat and bootstrapped variance
  results.temp <- est.delta.B.boot(s1 = study.B$S[study.B$G == 1], s0 = study.B$S[study.B$G == 0])
  y.pred.control.boot.estimates$delta.B[i] <- results.temp$delta.B
  y.pred.control.boot.estimates$se.delta.B[i] <- results.temp$se.delta.B
  y.pred.control.boot.estimates$z[i] <- results.temp$delta.B / results.temp$se.delta.B
  
  print(i)
  end_time <- Sys.time()
  print(end_time-start_time)
}

# compile and save table 
df <- data.frame("true.delta.B" = rep(get.truth(setting)$delta, 4),
                 "avg.est.delta.B" = c(mean(y.obs.estimates$delta.B), mean(y.pred.separate.estimates$delta.B), mean(y.pred.control.estimates$delta.B), mean(y.pred.control.boot.estimates$delta.B)),
                 "prop.null.reject" = c(sum(abs(y.obs.estimates$z) > 1.96) / num.sim, sum(abs(y.pred.separate.estimates$z) > 1.96) / num.sim, sum(abs(y.pred.control.estimates$z) > 1.96) / num.sim, sum(abs(y.pred.control.boot.estimates$z) > 1.96) / num.sim))
df <- t(df) 
colnames(df) <- c("observed.y", "predict.y.separate", "predict.y.control", "predict.y.control.boot")
write.table(df, paste("etsi.outputfile",setting, "_030724",parallel.num,".txt",sep=""), quote = FALSE, row.names = TRUE)
