set.seed(1)
rm(list= ls())
source("functions.R")
library(mvtnorm)
library(mcmcse)

############################################################
######### True ESS/n vs rho for diff settings @#############
################# in bivariate Normal case##################
############################################################


p <- 2
lambdas <- c(0.1, 0.8, 0.5)
rho_max <- c(.6, .8, .7)
Lambda <- diag(c(3,1))
Upsilon <- diag(c(3,3))

for(i in 1:3)
{
  Lambda[1,2] <- Lambda[2,1] <- lambdas[i]*sqrt(prod(diag(Lambda)))
  rhos <- seq(0, rho_max[i], .1)
  true.ess <- matrix(0, nrow = p+1, ncol = length(rhos))
  for (t in 1:length(rhos)){
    rho <- rhos[t]
    Upsilon[1,2] <- Upsilon[2,1] <- rho*sqrt(prod(diag(Upsilon)))
    print(det(2*solve(Lambda) - solve(Upsilon)))
    true.var <- var.is(Upsilon, Lambda)  ## N*true IS variance
    var.f <- Lambda  
    true.ess[1,t] <- (det(var.f)/det(true.var))^(1/p)
    true.ess[2:(p+1),t] <- diag(var.f)/diag(true.var)
  }
  save(true.ess, file = paste("trueESSvsRho_p2_setting", i, ".Rdata", sep=""))
}


#####################################################
###### Univarate and multivate true ESS vs rho  ######
########### for dimensions p = 2, 3, 10 ##############
#####################################################

dims <- c(2,3,10)

for (i in 1:3)
{
  p <- dims[i]
  
  sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
  lambda <- .6
  Lambda <- diag(sigma)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
      Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
    }
  }
  
  tau <- c(rep(2,ceiling(p/2)), rep(((p+1)/p), floor(p/2)))
  rho <- seq(.3, .7, .01)
  true.ess <- matrix(0, nrow = p+1, ncol = length(rho))
  for (t in 1:length(rho)){
    corr <- rho[t]
    Upsilon <- diag(tau)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        Upsilon[i,j] <- corr*sqrt(tau[i]*tau[j])
        Upsilon[j,i] <- corr*sqrt(tau[i]*tau[j])
      }
    }
    true.var <- var.is(Upsilon, Lambda)  ## N*true IS variance
    var.f <- Lambda  
    true.ess[1,t] <- (det(var.f)/det(true.var))^(1/p)
    true.ess[2:(p+1),t] <- diag(var.f)/diag(true.var)
  }
  save(true.ess, file = paste("trueESSvsRho_p", p, ".Rdata", sep=""))
}




##########################################################
##### Termination for diff epsilon in settings-1,2,3 #####
##########################################################

step <- 1e2
loop <- 20
dims <- c(2,3,10)
epsilons <- c(.1, .08, .06, .04, .03)

for (i in 1:3)
{
  p <- dims[i] # dimension
  print(paste("Dimension p=", p, sep=""))
  #target covariance
  sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
  lambda <- .6
  Lambda <- diag(sigma)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
      Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
    }
  }
  
  # proposal covariance
  tau <- c(rep(2,ceiling(p/2)), rep(((p+1)/p), floor(p/2)))
  rho <- .5
  Sigma <- diag(tau)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Sigma[i,j] <- rho*sqrt(tau[i]*tau[j])
      Sigma[j,i] <- rho*sqrt(tau[i]*tau[j])
    }
  }
  
  # Validity check
  foo <- 2*solve(Lambda) - solve(Sigma)
  min(eigen(foo)$values)
  
  ### Epsilon = 0.1
  min_ess <- minESS(p, eps = epsilons[1])
  N_min <- round(min_ess/2)
  all_ESS1 <- list()
  all_ESS1[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=0, p=p)
  for (i in 1:p){
    min_ess <- minESS(1, eps = epsilons[1])
    all_ESS1[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=i, p=p)
  }

  ### Epsilon = 0.08
  min_ess <- minESS(p, eps = epsilons[2])
  N_min <- round(min_ess/2)
  all_ESS2 <- list()
  all_ESS2[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=0, p=p)
  for (i in 1:p){
    min_ess <- minESS(1, eps = epsilons[2])
    all_ESS2[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=i, p=p)
  }
  
  ### Epsilon = 0.06
  min_ess <- minESS(p, eps = epsilons[3])
  N_min <- round(min_ess/2)
  all_ESS3 <- list()
  all_ESS3[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=0, p=p)
  for (i in 1:p){
    min_ess <- minESS(1, eps = epsilons[3])
    all_ESS3[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=i, p=p)
  }
  
  ### Epsilon = 0.04
  min_ess <- minESS(p, eps = epsilons[4])
  N_min <- round(min_ess/2)
  all_ESS4 <- list()
  all_ESS4[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=0, p=p)
  for (i in 1:p){
    min_ess <- minESS(1, eps = epsilons[4])
    all_ESS4[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, fun=fun, Lambda, h=i, p=p)
  }
  
  ### Epsilon = 0.03
  min_ess <- minESS(p, eps = epsilons[5])
  N_min <- round(min_ess/2)
  all_ESS5 <- list()
  all_ESS5[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, fun = fun, h=0, p=p)
  for (i in 1:p){
    min_ess <- minESS(1, eps = epsilons[5])
    all_ESS5[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, fun=fun, Lambda, h=i, p=p)
  }
  
  min_ess_eps <- apply(array(epsilons, dim = c(1,length(epsilons))), MARGIN =1, FUN=minESS, p=p, alpha=.05) # c(minESS(p, eps = .1), minESS(p, eps=.08), minESS(p, eps=.06), minESS(p, eps=.04), minESS(p, eps=.03))
  truth <- min_ess_eps/((det(Lambda)/det(var.is(Sigma, Lambda)))^(1/p))
  
  save(truth, all_ESS1, all_ESS2, all_ESS3, all_ESS4, all_ESS5, file = paste("ESS_eps_p", p, ".Rdata", sep = ""))
}


