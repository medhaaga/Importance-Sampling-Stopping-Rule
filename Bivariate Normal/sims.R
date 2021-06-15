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
  save(true.ess.uis, true.ess.snis, file = paste("trueESSvsRho_p2_setting", i, ".Rdata", sep=""))
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

step <- 5e2
loop <- 20
dims <- c(2,10)
epsilons <- c(.1, .06, .04, .02)
lambdas <- c(0.1, 0.8, 0.5)
rhos <- c(0.1, 0.7, 0.5)

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
    
    ### Epsilon = 0.1
    min_ess <- minESS(p, eps = epsilons[1])
    N_min <- round(min_ess/2)
    all_ESS1 <- list()
    all_ESS1[[1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=0, p=p)
    # for (i in 1:p){
    #   min_ess <- minESS(1, eps = epsilons[1])
    #   all_ESS1[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=i, p=p)
    # }
    # 
    ### Epsilon = 0.06
    min_ess <- minESS(p, eps = epsilons[2])
    N_min <- round(min_ess/2)
    all_ESS2 <- list()
    all_ESS2[[1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=0, p=p)
    # for (i in 1:p){
    #   min_ess <- minESS(1, eps = epsilons[2])
    #   all_ESS2[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=i, p=p)
    # }
    # 
    ### Epsilon = 0.04
    min_ess <- minESS(p, eps = epsilons[3])
    N_min <- round(min_ess/2)
    all_ESS3 <- list()
    all_ESS3[[1]] <- is_ESS(min_ess, step, loop, N_min, Upsilon, Lambda, mu, fun = fun, h=0, p=p)
    # for (i in 1:p){
    #   min_ess <- minESS(1, eps = epsilons[3])
    #   all_ESS3[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, Upsilon, Lambda, mu, fun, h=i, p=p)
    # }
    # 
    ### Epsilon = 0.02
    min_ess <- minESS(p, eps = epsilons[4])
    N_min <- round(min_ess/2)
    all_ESS4 <- list()
    all_ESS4[[1]] <- is_ESS(min_ESS, step, loop, N_min, Upsilon, Lambda, mu, fun, h=0, p)
    # for (i in 1:p){
    #   min_ess <- minESS(1, eps = epsilons[4])
    #   all_ESS4[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, Upsilon, Lambda, mu, fun, h=i, p)
    # }
    # 
    min_ess_eps <- apply(array(epsilons, dim = c(1,length(epsilons))), MARGIN =1, FUN=minESS, p=p, alpha=.05) # c(minESS(p, eps = .1), minESS(p, eps=.08), minESS(p, eps=.06), minESS(p, eps=.04), minESS(p, eps=.03))
    truth.uis <- min_ess_eps/((det(Lambda)/det(var.uis(Upsilon, Lambda, mu)))^(1/p))
    truth.snis <- min_ess_eps/((det(Lambda)/det(var.snis(Upsilon, Lambda, mu)))^(1/p))
    
    save(truth.uis, truth.snis, all_ESS1, all_ESS2, all_ESS3, all_ESS4, file = paste("ESS_eps_p", p, "_setting", a, ".Rdata", sep = ""))
    
  }
}
 

# p <- 1
# sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
# lambda <- .5
# Lambda <- diag(sigma)
# for (i in 1:(p-1)){
#   for (j in (i+1):p){
#     Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
#     Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
#   }
# }
# 
# Lambda <- matrix(1)
# # proposal covariance
# tau <- c(rep(2,ceiling(p/2)), rep(((p+1)/p), floor(p/2)))
# rho <- .5
# Upsilon <- diag(tau)
# for (i in 1:(p-1)){
#   for (j in (i+1):p){
#     Upsilon[i,j] <- rho*sqrt(tau[i]*tau[j])
#     Upsilon[j,i] <- rho*sqrt(tau[i]*tau[j])
#   }
# }
# Upsilon <- matrix(2)
# 
# # Validity check
# foo <- 2*solve(Lambda) - solve(Upsilon)
# min(eigen(foo)$values)
# 
# mu <- rep(2, p)
# 
# samp_size <- seq(1e2, 1e4, 1e3)
# 
# truevar.uis <- var.uis(Upsilon, Lambda, mu)
# truevar.snis <- var.snis(Upsilon, Lambda, mu)
# trueESS.uis <- (det(Lambda)/det(truevar.uis))^(1/p)
# trueESS.snis <- (det(Lambda)/det(truevar.snis))^(1/p)
# trueESS.uis.num <- det(Lambda)^(1/p)
# trueESS.uis.denom <- det(truevar.uis)^(1/p)
# trueESS.snis.num <- det(Lambda)^(1/p)
# trueESS.snis.denom <- det(truevar.snis)^(1/p)
# reps <- 1e2
# samps <- length(samp_size)
# ESS.kong <- matrix(0, ncol = samps, nrow = reps)
# ESS.snis <- matrix(0, ncol = samps, nrow = reps)
# ESS.uis <- matrix(0, ncol = samps, nrow = reps)
# ESS.snis.num <- matrix(0, ncol = samps, nrow = reps)
# ESS.uis.num <- matrix(0, ncol = samps, nrow = reps)
# ESS.snis.denom <- matrix(0, ncol = samps, nrow = reps)
# ESS.uis.denom <- matrix(0, ncol = samps, nrow = reps)
# 
# 
# for (r in 1:reps)
# {
#   print(r)
#   for (i in 1:samps)
#   {
#     is <- rmvnorm(n=samp_size[i], mean = mu, sigma = Upsilon)
#     run_weights <- dmvnorm(is, mean = mu, sigma = Lambda)/dmvnorm(is, mean = mu, sigma = Upsilon)
#     norm_weights <- run_weights/sum(run_weights)
#     snis <- colSums(is * norm_weights)
#     uis <- colSums(is * run_weights)/samp_size[i]
#     
#     ESS.kong[r,i] <- 1/(samp_size[i]*sum(norm_weights^2))
#     foo1 <- matrix(uis, nrow = samp_size[i], ncol = p, byrow = TRUE)
#     foo2 <- matrix(snis, nrow = samp_size[i], ncol = p, byrow = TRUE)
#     
#     varIbar.uis <- (t(is - foo1) %*% (run_weights * (is - foo1)))/(samp_size[i]^2)
#     varIbar.snis <- (t(is - foo2) %*% (norm_weights * (is - foo2)))/samp_size[i]
#     emp.var.uis <-  (t(run_weights*is - foo1) %*% (run_weights*is - foo1))/(samp_size[i]^2)
#     emp.var.snis <- (t(norm_weights*(is - foo2)) %*% (norm_weights * (is - foo2)))
#     ESS.uis.num[r,i] <- det(varIbar.uis)^(1/p)
#     ESS.uis.denom[r,i] <- det(emp.var.uis)^(1/p)
#     ESS.uis[r,i] <- ((det(varIbar.uis)/det(emp.var.uis))^(1/p))
#     ESS.snis.num[r,i] <- det(varIbar.snis)^(1/p)
#     ESS.snis.denom[r,i] < det(emp.var.snis)^(1/p)
#     ESS.snis[r,i] <- ((det(varIbar.snis)/det(emp.var.snis))^(1/p))
#   }
# }
# 
# save(ESS.kong, trueESS.uis, trueESS.snis, ESS.uis.num, ESS.uis.denom, ESS.uis, ESS.snis.num, ESS.snis.denom, ESS.snis, file = "uisVSsnis.Rdata")
