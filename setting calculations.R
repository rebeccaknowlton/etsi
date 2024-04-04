#####################################################
### this is the totally parametric/linear version ###
#####################################################

possible.beta = seq(from = 0.1, to = 8, by = 0.1)
beta.grid = data.frame(beta1 = possible.beta, 
                       beta2 = possible.beta,
                       beta3 = possible.beta,
                       beta5 = possible.beta)
big.beta.grid = data.frame(expand.grid(beta.grid))
big.beta.grid$delta.lower = big.beta.grid$beta1 + (big.beta.grid$beta2 + big.beta.grid$beta3) * 1 - big.beta.grid$beta2 * 0.25 + big.beta.grid$beta5 * 2
big.beta.grid$delta.upper = big.beta.grid$beta1 + (big.beta.grid$beta2 + big.beta.grid$beta3) * 1 - big.beta.grid$beta2 * 0.25 + big.beta.grid$beta5 * 6
big.beta.grid$delta.s.lower = big.beta.grid$beta1 + big.beta.grid$beta3 * 0.25 + big.beta.grid$beta5 * 2
big.beta.grid$delta.s.upper = big.beta.grid$beta1 + big.beta.grid$beta3 * 0.25 + big.beta.grid$beta5 * 6
big.beta.grid$R.lower = 1 - big.beta.grid$delta.s.lower / big.beta.grid$delta.lower 
big.beta.grid$R.upper = 1 - big.beta.grid$delta.s.upper / big.beta.grid$delta.upper

# what's the largest possible difference in R?
idx = which(big.beta.grid$R.upper - big.beta.grid$R.lower == min(big.beta.grid$R.upper - big.beta.grid$R.lower))
big.beta.grid[idx,] 
# about 0.26

# okay, so filter on parts of the grid that have R varying by at least 0.25, the higher value of R > 0.6, and find the minimum overall treatment effect (so power isn't too high)
idx = which((big.beta.grid$R.upper - big.beta.grid$R.lower < -0.25) & (big.beta.grid$R.lower > 0.6))
big.beta.grid[idx,][which( big.beta.grid[idx,]$delta.lower == min(big.beta.grid[idx,]$delta.lower)),]

# these would be the parameters then
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

delta = beta1 + (beta2 + beta3) * s1.scale * s1.shape - beta2 * s0.scale * s0.shape + beta5 * c(2, 6)
delta.s = beta1 + beta3 * s0.scale * s0.shape + beta5 * c(2, 6)
R = 1 - delta.s / delta

delta 
delta.s
R

######################################################################
### can I get something better if instead beta5 is for a term W^2? ###
######################################################################

possible.beta = seq(from = 0.1, to = 8, by = 0.1)
beta.grid = data.frame(beta1 = possible.beta, 
                       beta2 = possible.beta,
                       beta3 = possible.beta,
                       beta5 = possible.beta)
big.beta.grid = data.frame(expand.grid(beta.grid))
big.beta.grid$delta.lower = big.beta.grid$beta1 + (big.beta.grid$beta2 + big.beta.grid$beta3) * 1 - big.beta.grid$beta2 * 0.25 + big.beta.grid$beta5 * 4
big.beta.grid$delta.upper = big.beta.grid$beta1 + (big.beta.grid$beta2 + big.beta.grid$beta3) * 1 - big.beta.grid$beta2 * 0.25 + big.beta.grid$beta5 * 36
big.beta.grid$delta.s.lower = big.beta.grid$beta1 + big.beta.grid$beta3 * 0.25 + big.beta.grid$beta5 * 4
big.beta.grid$delta.s.upper = big.beta.grid$beta1 + big.beta.grid$beta3 * 0.25 + big.beta.grid$beta5 * 36
big.beta.grid$R.lower = 1 - big.beta.grid$delta.s.lower / big.beta.grid$delta.lower 
big.beta.grid$R.upper = 1 - big.beta.grid$delta.s.upper / big.beta.grid$delta.upper

# what's the largest possible difference in R?
idx = which(big.beta.grid$R.upper - big.beta.grid$R.lower == min(big.beta.grid$R.upper - big.beta.grid$R.lower))
big.beta.grid[idx,] 
# about 0.49

# okay, so filter on parts of the grid that have R varying by at least 0.3, the higher value of R > 0.6, and find the minimum overall treatment effect (so power isn't too high)
idx = which((big.beta.grid$R.upper - big.beta.grid$R.lower < -0.3) & (big.beta.grid$R.lower > 0.6))
big.beta.grid[idx,][which( big.beta.grid[idx,]$delta.lower == min(big.beta.grid[idx,]$delta.lower)),]


# these would be the parameters then
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


delta = beta1 + (beta2 + beta3) * s1.scale * s1.shape - beta2 * s0.scale * s0.shape + beta5 * c(4, 36)
delta.s = beta1 + beta3 * s0.scale * s0.shape + beta5 * c(4, 36)
R = 1 - delta.s / delta

delta 
delta.s
R

