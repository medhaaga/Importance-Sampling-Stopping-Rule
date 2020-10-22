###############################################################
#### The code plots ESS comparison plots with varying #########
#### sigma and rho using determinant and trace for low ########
#### and high correlation of the target distribution. #########
###############################################################

set.seed(1)
library(mvtnorm)

p <- 2
N <- 1e4

var.is <- function(Sigma, Lambda){
  foo <- 2*solve(Lambda) - solve(Sigma)
  return((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo)))) * solve(foo))
}

obj.det <- function(par, Lambda){
  sigma <- par[1]
  rho <- par[2]
  Sigma <- matrix(c(sigma, rho*sigma, rho*sigma, sigma), 2, 2)
  foo <- 2*solve(Lambda) - solve(Sigma)
  return((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo))))^2 / det(foo))
}

obj.trace <- function(par, Lambda){
  sigma <- par[1]
  rho <- par[2]
  Sigma <- matrix(c(sigma, rho*sigma, rho*sigma, sigma), 2, 2)
  foo <- 2*solve(Lambda) - solve(Sigma)
  return((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo)))) * sum(diag(solve(foo))))
}

###########################################################
############ Using determinants ###########################
###########################################################

############### Low correlation of target #################

lambda <- .1
Lambda <- matrix(c(1, lambda, lambda, 1), 2, 2)
optim(c(1,0), fn = obj.det, Lambda = Lambda, method = "BFGS")

############################################
###### Varying sigma #######################
############################################

prop_sigma <- seq((1/sqrt(2))+.1, 2*sqrt((p+1)/p), .1)  #proposal sigmas
ess <- numeric(length = length(prop_sigma))
ess.kong <- numeric(length = length(prop_sigma))
ess.boot <- numeric(length = length(prop_sigma))

###### Simple Importance Sampling ############

rho <- .1

for (y in 1:length(prop_sigma)){
  
  Sigma <- matrix(c(prop_sigma[y], rho* prop_sigma[y], rho* prop_sigma[y], prop_sigma[y]), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*det(var.f)*(1/det(true.var))
  
  B <- N
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * det(varIbar)/det( var(boot_means) + tcrossprod(colMeans(boot_means) - snis))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/det_lowCorr_sigma_IS.pdf")
plot(prop_sigma, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess/N, lty=2, col = "red")
lines(prop_sigma, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=prop_sigma[ess == max(ess)], col = "pink")
dev.off()

############################################
###### Varying rho #######################
############################################

prop_rho <- seq(0, 1, .1)  #proposal sigmas
ess <- numeric(length = length(prop_rho))
ess.kong <- numeric(length = length(prop_rho))
ess.boot <- numeric(length = length(prop_rho))

###### Simple Importance Sampling ############

sigma <- 1.5

for (y in 1:length(prop_rho)){
  
  Sigma <- matrix(c(sigma, prop_rho[y]*sigma, prop_rho[y]*sigma, sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*det(var.f)*(1/det(true.var))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * det(varIbar)/det( var(boot_means) + tcrossprod(colMeans(boot_means) - snis))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/det_lowCorr_rho_IS.pdf")
plot(prop_rho, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "rho", ylab = "ESS/N")
lines(prop_rho, ess/N, lty=2, col = "red")
lines(prop_rho, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=prop_rho[ess == max(ess)], col = "pink")


###########################################################
################ High correlation #########################
###########################################################

lambda <- .9
Lambda <- matrix(c(1, lambda, lambda, 1), 2, 2)
optim(par = c(1,.7), fn = obj.det, Lambda = Lambda)


############################################
###### Varying sigma #######################
############################################

prop_sigma <- seq((1/sqrt(2))+.3, 2*sqrt((p+1)/p), .1)  #proposal sigmas
ess <- numeric(length = length(prop_sigma))
ess.kong <- numeric(length = length(prop_sigma))
ess.boot <- numeric(length = length(prop_sigma))

###### Simple Importance Sampling ############

rho <- .9

for (y in 1:length(prop_sigma)){
  
  Sigma <- matrix(c(prop_sigma[y], rho*prop_sigma[y], rho*prop_sigma[y], prop_sigma[y]), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*det(var.f)*(1/det(true.var))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * det(varIbar)/det( var(boot_means) + tcrossprod(colMeans(boot_means) - snis))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/det_highCorr_sigma_IS.pdf")
plot(prop_sigma, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess/N, lty=2, col = "red")
lines(prop_sigma, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = prop_sigma[ess == max(ess)], col = "pink")
dev.off()

############################################
###### Varying rho #######################
############################################

prop_rho <- seq(0, 1, .1)  #proposal sigmas
ess <- numeric(length = length(prop_rho))
ess.kong <- numeric(length = length(prop_rho))
ess.boot <- numeric(length = length(prop_rho))

###### Simple Importance Sampling ############

sigma <- 1.5

for (y in 1:length(prop_rho)){
  
  Sigma <- matrix(c(sigma, prop_rho[y]*sigma, prop_rho[y]*sigma, sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*det(var.f)*(1/det(true.var))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * det(varIbar)/det( var(boot_means) + tcrossprod(colMeans(boot_means) - snis))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/det_highCorr_rho_IS.pdf")
plot(prop_rho, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "rho", ylab = "ESS/N")
lines(prop_rho, ess/N, lty=2, col = "red")
lines(prop_rho, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = prop_rho[ess == max(ess)], col = "pink")
dev.off()


#############################################
#############################################
########## Using trace ######################
#############################################
#############################################


############################################
######## Low correlation of target #########
############################################

lambda <- .1
Lambda <- matrix(c(1, lambda, lambda, 1), 2, 2)
optim(c(1.5,0.1), fn = obj.trace, Lambda = Lambda, method = "BFGS")

############################################
###### Varying sigma #######################
############################################

prop_sigma <- seq((1/sqrt(2))+.1, 2*sqrt((p+1)/p), .1)  #proposal sigmas
ess <- numeric(length = length(prop_sigma))
ess.kong <- numeric(length = length(prop_sigma))
ess.boot <- numeric(length = length(prop_sigma))

##############################################
###### Simple Importance Sampling ############
##############################################

rho <- .1

for (y in 1:length(prop_sigma)){
  
  Sigma <- matrix(c(prop_sigma[y], rho, rho, prop_sigma[y]), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*sum(diag(var.f))*(1/sum(diag(true.var)))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * sum(diag(varIbar))/sum(diag( var(boot_means) + tcrossprod(colMeans(boot_means) - snis)))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/trace_lowCorr_sigma_IS.pdf")
plot(prop_sigma, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess/N, lty=2, col = "red")
lines(prop_sigma, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=prop_sigma[ess == max(ess)], col = "pink")
dev.off()

############################################
###### Varying rho #######################
############################################

prop_rho <- seq(0, 1, .1)  #proposal sigmas
ess <- numeric(length = length(prop_rho))
ess.kong <- numeric(length = length(prop_rho))
ess.boot <- numeric(length = length(prop_rho))

##############################################
###### Simple Importance Sampling ############
##############################################

sigma <- 1.5

for (y in 1:length(prop_rho)){
  
  Sigma <- matrix(c(sigma, prop_rho[y], prop_rho[y], sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*sum(diag(var.f))*(1/sum(diag(true.var)))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * sum(diag(varIbar))/sum(diag( var(boot_means) + tcrossprod(colMeans(boot_means) - snis)))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/trace_lowCorr_rho_IS.pdf")
plot(prop_rho, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "rho", ylab = "ESS/N")
lines(prop_rho, ess/N, lty=2, col = "red")
lines(prop_rho, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=prop_rho[ess == max(ess)], col = "pink")
dev.off()

###########################################################
################ High correlation #########################
###########################################################

lambda <- .8
Lambda <- matrix(c(1, lambda, lambda, 1), 2, 2)
optim(c(1,.8), fn = obj.trace, Lambda = Lambda, method = "BFGS")

############################################
###### Varying sigma #######################
############################################

prop_sigma <- seq((1/sqrt(2))+.3, 2*sqrt((p+1)/p), .05)  #proposal sigmas
ess <- numeric(length = length(prop_sigma))
ess.kong <- numeric(length = length(prop_sigma))
ess.boot <- numeric(length = length(prop_sigma))

##############################################
###### Simple Importance Sampling ############
##############################################

rho <- .8

for (y in 1:length(prop_sigma)){
  
  Sigma <- matrix(c(prop_sigma[y], rho*prop_sigma[y], rho*prop_sigma[y], prop_sigma[y]), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*sum(diag(var.f))*(1/sum(diag(true.var)))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * sum(diag(varIbar))/sum(diag( var(boot_means) + tcrossprod(colMeans(boot_means) - snis)))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/trace_highCorr_sigma_IS.pdf")
plot(prop_sigma, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess/N, lty=2, col = "red")
lines(prop_sigma, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = prop_sigma[ess == max(ess)], col = "pink")
dev.off()

############################################
###### Varying rho #######################
############################################

prop_rho <- seq(.6, .939, .02)  #proposal sigmas
ess <- numeric(length = length(prop_rho))
ess.kong <- numeric(length = length(prop_rho))
ess.boot <- numeric(length = length(prop_rho))

##############################################
###### Simple Importance Sampling ############
##############################################

sigma <- 1.5

for (y in 1:length(prop_rho)){
  
  Sigma <- matrix(c(sigma, prop_rho[y]*sigma, prop_rho[y]*sigma, sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  # true Var_p(h(x))
  is <- rmvnorm(N, mean=rep(0,p), sigma = Sigma )  ## Importance samples
  weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
  norm_weights <- weights/sum(weights)
  snis <- colMeans(is*weights)
  
  ess.kong[y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
  ess[y] <- N*sum(diag(var.f))*(1/sum(diag(true.var)))
  
  B <- 500
  boot_means <- matrix(0, nrow = B, ncol = p)
  
  varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N)    ## empirical Var_p(h(x))
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
  }
  ess.boot[y] <- N * sum(diag(varIbar))/sum(diag( var(boot_means) + tcrossprod(colMeans(boot_means) - snis)))
  print(paste(ess[y], ess.kong[y], ess.boot[y]))
}

pdf(file = "Out_gaussian/trace_highCorr_rho_IS.pdf")
plot(prop_rho, ess.kong/N, type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.kong/N, ess.boot/N, ess/N), xlab = "rho", ylab = "ESS/N")
lines(prop_rho, ess/N, lty=2, col = "red")
lines(prop_rho, ess.boot/N, lty=3, col = "blue")
legend("topright", legend = c("Kong", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = prop_rho[ess == max(ess)], col = "pink")
dev.off()

