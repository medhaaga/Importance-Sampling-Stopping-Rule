set.seed(1)
library(mvtnorm)


## 81 points on the grid (-2,2)*(-2,2) are taken at a step length of .5 in both directions
prop_mu <- matrix(0, nrow = 81, ncol = 2)  #proposal mu 
prop_mu[,1] <- c(rep(-2,9), rep(-1.5,9), rep(-1,9), rep(-.5,9), rep(-0,9), rep(.5,9), rep(1,9), rep(1.5,9), rep(2,9))
prop_mu[,2] <- rep(seq(-2,2,.5), 9)

N <- c(1e2, 1e3, 1e4)    #Number of samples
ess <- matrix(0, nrow = 3, ncol = 81)
ess.snis <- matrix(0, nrow = 3, ncol = 81)
ess.boot <- matrix(0, nrow = 3, ncol = 81)


######### function calculates true N*Variance of IS estiamtor ###########
var.is <- function(mu){
  mu <- as.vector(mu)
  var.is <- as.numeric(exp(t(mu) %*% mu)) * (diag(c(1,1)) + (mu %*% t(mu)))
  return (var.is)
}

######### function calculates true N*Variance of SNIS estiamtor ###########
var.snis <- function(mu){
  mu <- as.vector(mu)
  var.snis <- as.numeric(exp(t(mu) %*% mu)) * (diag(c(1,1)) + (mu %*% t(mu)))
  return (var.snis)
}

##############################################
###### Simple Importance Sampling ############
##############################################

for (x in 1:length(N)){
  for (y in 1:nrow(prop_mu)){
    
    true.var.is <- var.is(prop_mu[y,])  ## N*true IS variance
    var.f <- diag(c(1,1))  # true Var_p(h(x))
    is <- rmvnorm(n=N[x], mean=prop_mu[y,], sigma = diag(2))   ## Importance samples
    weights <- dmvnorm(is, mean=c(0,0), sigma = diag(2))/dmvnorm(is, mean=prop_mu[y,], sigma= diag(2))   
    norm_weights <- weights/sum(weights)
    snis <- colMeans(is*weights)
    
    ess.snis[x,y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
    ess[x,y] <- N[x]*(1/det(true.var.is))
    
    B <- 500
    boot_means <- matrix(0, nrow = B, ncol = 2)
    
    varIbar <- cov.wt(x = is-snis, wt  = weights)$cov/(N[x])   ## empirical Var_p(h(x))
    for (i in 1:B){
      sample_idx <- sample(1:N[x], size= N[x], replace = TRUE)
      boot_means[i,] <- colMeans(is[sample_idx,] * weights[sample_idx])
    }
    ess.boot[x,y] <- N[x] * det(varIbar)/det(var(boot_means))
    print(paste(ess[x,y], ess.snis[x,y], ess.boot[x,y]))
  }
}
save(ess.snis, ess, ess.boot, file = paste("Out_gaussian/multi_gaussian-IS_m", m, ".Rdata", sep = ""))
load(file = paste("Out_gaussian/multi_gaussian-IS_m", m, ".Rdata", sep = ""))
# plots for x coordinate fixed at 0 and y varying from -2 to 2.

## N=1e2
pdf(file = paste("Out_gaussian/multi_gaussian-IS_m", m, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,3))
plot(prop_mu[37:45,2], ess.snis[1,37:45]/N[1], type = "l", main = "No of samples = 1e2"
     , ylim = range(ess.snis[1,37:45]/N[1], ess.boot[1,37:45]/N[1], ess[1,37:45]/N[1]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[1,37:45]/N[1], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[1,37:45]/N[1], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=0, col = "pink")


## N=1e3
plot(prop_mu[37:45,2], ess.snis[2,37:45]/N[2], type = "l", main = "No of samples = 1e3"
     , ylim = range(ess.snis[2,37:45]/N[2], ess.boot[2,37:45]/N[2], ess[2,37:45]/N[2]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[2,37:45]/N[2], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[2,37:45]/N[2], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=0, col = "pink")


## N=1e4
plot(prop_mu[37:45,2], ess.snis[3,37:45]/N[3], type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.snis[3,37:45]/N[3], ess.boot[3,37:45]/N[3], ess[3,37:45]/N[3]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[3,37:45]/N[3], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[3,37:45]/N[3], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = 0, col = "pink")
dev.off()




##################################################
########## Weighted Importance Sampling# #########
##################################################

for (x in 1:length(N)){
  for (y in 1:nrow(prop_mu)){
    
    true.var.is <- var.snis(prop_mu[y,])
    var.f <- diag(2)
    is <- rmvnorm(n=N[x], mean=prop_mu[y,], sigma = diag(2))   ## Importance samples
    weights <- dmvnorm(is, mean=c(0,0), sigma = diag(2))/dmvnorm(is, mean=prop_mu[y,], sigma= diag(2))   
    norm_weights <- weights/sum(weights)
    snis <- colMeans(is*weights)/sum(weights)
    
    ess.snis[x,y] <- 1/sum(norm_weights^2)
    ess[x,y] <- N[x]*(1/det(true.var.is))
    
    B <- 500
    boot_means <- matrix(0, nrow = B, ncol = 2)
    
    varIbar <- cov.wt(x = is-snis, wt  = norm_weights)$cov/(N[x])   ## empirical Var_p(h(x))
    for (i in 1:B){
      sample_idx <- sample(1:N[x], size= N[x], replace = TRUE)
      boot_means[i,] <- colSums(is[sample_idx,] * weights[sample_idx])/sum(weights[sample_idx])
    }
    ess.boot[x,y] <- N[x] * det(varIbar)/det(var(boot_means))
    print(paste(ess[x,y], ess.snis[x,y], ess.boot[x,y]))
  }
  
}
save(ess.snis, ess, ess.boot, file = paste("Out_gaussian/multi_gaussian-SNIS_m", m, ".Rdata", sep = ""))

load(file = paste("Out_gaussian/multi_gaussian-SNIS_m", m, ".Rdata", sep = ""))

# plots for x coordinate fixed at 0 and y varying from -2 to 2.

## N=1e2
pdf(file = paste("Out_gaussian/multi_gaussian-SNIS_m", m, ".pdf", sep = ""), height = 4, width = 10)

par(mfrow = c(1,3))
plot(prop_mu[37:45,2], ess.snis[1,37:45]/N[1], type = "l", main = "No of samples = 1e2"
     , ylim = range(ess.snis[1,37:45]/N[1], ess.boot[1,37:45]/N[1], ess[1,37:45]/N[1]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[1,37:45]/N[1], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[1,37:45]/N[1], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=0, col = "pink")


## N=1e3
plot(prop_mu[37:45,2], ess.snis[2,37:45]/N[2], type = "l", main = "No of samples = 1e3"
     , ylim = range(ess.snis[2,37:45]/N[2], ess.boot[2,37:45]/N[2], ess[2,37:45]/N[2]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[2,37:45]/N[2], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[2,37:45]/N[2], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=0, col = "pink")


## N=1e4
plot(prop_mu[37:45,2], ess.snis[3,37:45]/N[3], type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.snis[3,37:45]/N[3], ess.boot[3,37:45]/N[3], ess[3,37:45]/N[3]), xlab = "y coordinate", ylab = "ESS/N")
lines(prop_mu[37:45,2], ess[3,37:45]/N[3], lty=2, col = "red")
lines(prop_mu[37:45,2], ess.boot[3,37:45]/N[3], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = 0, col = "pink")
dev.off()


