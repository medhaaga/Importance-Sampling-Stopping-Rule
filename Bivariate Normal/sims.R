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
lambdas <- c(0.1, 0.5, 0.8)
rho_max <- c(.6, .7, .8)
Lambda <- diag(c(2,1))
Upsilon <- diag(c(2,2))
mu <- rep(2, p)

for(i in 1:3)
{
  Lambda[1,2] <- Lambda[2,1] <- lambdas[i]*sqrt(prod(diag(Lambda)))
  rhos <- seq(0, rho_max[i], .1)
  true.ess.uis <- matrix(0, nrow = p+1, ncol = length(rhos))
  true.ess.snis <- matrix(0, nrow = p+1, ncol = length(rhos))
  for (t in 1:length(rhos)){
    rho <- rhos[t]
    Upsilon[1,2] <- Upsilon[2,1] <- rho*sqrt(prod(diag(Upsilon)))
    print(det(2*solve(Lambda) - solve(Upsilon)))
    true.var.uis <- var.uis(Upsilon, Lambda, mu)
    true.var.snis <- var.snis(Upsilon, Lambda, mu)
    var.f <- Lambda  
    true.ess.uis[1,t] <- (det(var.f)/det(true.var.uis))^(1/p)
    true.ess.uis[2:(p+1),t] <- diag(var.f)/diag(true.var.uis)
    true.ess.snis[1,t] <- (det(var.f)/det(true.var.snis))^(1/p)
    true.ess.snis[2:(p+1),t] <- diag(var.f)/diag(true.var.snis)
  }
  save(true.ess.uis, true.ess.snis, file = paste("Out/trueESSvsRho_p2_setting", i, ".Rdata", sep=""))
}



##########################################################
##### Termination for diff epsilon in settings-1,2,3 #####
##########################################################

start.time <- Sys.time()
step <- 5e2
loop <- 1e2
dims <- c(2,10)
epsilons <- c(.06, .04, .02)
lambdas <- c(0.1, 0.5, 0.8)
rhos <- c(0.1, 0.5, 0.7)

for (a in 1:3)
{
  for (b in 1:2)
  {
    p <- dims[b]
    mu <- rep(2, p)
    print(paste("Dimension p=", p, " setting ", a, sep=""))
    
    #target covariance
    sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
    lambda <- lambdas[a]
    Lambda <- diag(sigma)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
        Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
      }
    }
    
    # proposal covariance
    tau <- rep(2, p)
    rho <- rhos[a]
    Upsilon <- diag(tau)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        Upsilon[i,j] <- rho*sqrt(tau[i]*tau[j])
        Upsilon[j,i] <- rho*sqrt(tau[i]*tau[j])
      }
    }
    
    # Validity check
    foo <- 2*solve(Lambda) - solve(Upsilon)
    min(eigen(foo)$values)
    
    ### Epsilon = 0.06
    min_ess <- minESS(p, eps = epsilons[1])
    N_min <- round(min_ess)
    all_ESS1 <- list()
    all_ESS1[[1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=0, p=p)
    
    ### Epsilon = 0.04
    min_ess <- minESS(p, eps = epsilons[2])
    N_min <- round(min_ess)
    all_ESS2 <- list()
    all_ESS2[[1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=0, p=p)
    
    ### Epsilon = 0.02
    min_ess <- minESS(p, eps = epsilons[3])
    N_min <- round(min_ess)
    all_ESS3 <- list()
    all_ESS3[[1]] <- is_ESS(min_ESS, step, loop, N_min, Upsilon, Lambda, mu, fun, h=0, p)
    
    min_ess_eps <- apply(array(epsilons, dim = c(1,length(epsilons))), MARGIN =1, FUN=minESS, p=p, alpha=.05) # c(minESS(p, eps = .1), minESS(p, eps=.08), minESS(p, eps=.06), minESS(p, eps=.04), minESS(p, eps=.03))
    truth.uis <- min_ess_eps/((det(Lambda)/det(var.uis(Upsilon, Lambda, mu)))^(1/p))
    truth.snis <- min_ess_eps/((det(Lambda)/det(var.snis(Upsilon, Lambda, mu)))^(1/p))
    
    save(truth.uis, truth.snis, all_ESS1, all_ESS2, all_ESS3, file = paste("Out/ESS_eps_p", p, "_setting", a, ".Rdata", sep = ""))
    
  }
}
end.time <- Sys.time()
print(end.time - start.time)
