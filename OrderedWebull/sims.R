set.seed(1)
rm(list = ls())
library(mcmcse)
source("functions.R")


########## Fish dataset1 #################
t <- c(83.5, 91.0, 91.0, 97.0, 107.0, 109.5, 114.0, 115.41, 128.61, 133.53, 138.58, 140.0, 152.08, 155.1)
t <- (t-80)/100
failures <- c(0,6,3,3,2)
stress <- c(0, .3, .5, .7, 1)
m <- 3
n <- 14
r <- 14
p <- m+2

####### hyperparams for priors ###########
val <- .8
shape.a <- val 
shape.b <- val 
ord.a <- val 
ord.b <- val 
ord.A <- rep(1, m+1)
p_H1 <- p
p_H2 <- p-1

n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])


reps <- 5e2
samp_size <- seq(1e3, 1e4, 1e3)
samps <- length(samp_size)
kong <- matrix(0, nrow = samps, ncol = reps)
multESS1 <- matrix(0, nrow = samps, ncol = reps)
num.multESS1 <- matrix(0, nrow = samps, ncol = reps)
denom.multESS1 <- matrix(0, nrow = samps, ncol = reps)
uniESS1 <- array(0, dim = c(p, reps, samps))
multESS2 <- matrix(0, nrow = samps, ncol = reps)
num.multESS2 <- matrix(0, nrow = samps, ncol = reps)
denom.multESS2 <- matrix(0, nrow = samps, ncol = reps)
uniESS2 <- array(0, dim = c(4, reps, samps))

for (l in 1:reps)
{
  print(l)
  sim <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, samp_size[length(samp_size)])
  is.samp <- sim$samp
  run_weights <- sim$weights
  
  for (z in 1:samps)
  {
    M <- samp_size[z]
    H1 <- fun1(is.samp[1:M,])
    H2 <- fun2(is.samp[1:M,])
    norm_weights <- run_weights[1:M]/sum(run_weights[1:M])
    
    kong[z,l] <- 1/sum(norm_weights^2)
    
    snis1 <- colSums(H1 * norm_weights)
    snis2 <- colSums(H2 * norm_weights)
    
    varIbar1 <- (t(H1 - snis1) %*% (norm_weights * (H1 - snis1)))/M
    varIbar2 <- (t(H2 - snis2) %*% (norm_weights * (H2 - snis2)))/M
    
    emp.var1 <- (t(H1 - snis1) %*% (norm_weights^2 * (H1 - snis1)))
    emp.var2 <- (t(H2 - snis2) %*% (norm_weights^2 * (H2 - snis2)))
    
    multESS1[z,l] <- M*((det(varIbar1)/det(emp.var1))^(1/p_H1))
    num.multESS1[z,l] <- det(varIbar1)^(1/p_H1)
    denom.multESS1[z,l] <- det(emp.var1)^(1/p_H1)
    uniESS1[,l,z] <- M* diag(varIbar1)/ diag(emp.var1)
    
    multESS2[z,l] <- M*((det(varIbar2)/det(emp.var2))^(1/p_H2))
    num.multESS2[z,l] <- det(varIbar2)^(1/p_H2)
    denom.multESS2[z,l] <- det(emp.var2)^(1/p_H2)
    uniESS2[,l,z] <- M* diag(varIbar2)/ diag(emp.var2)
    
  }
}
save(kong, uniESS1, uniESS2, multESS1, multESS2, num.multESS1, num.multESS2, denom.multESS1, denom.multESS2, file = "ESSvsSampSize_objects_kong_uni_mult.Rdata")

###########################################################
############ Termination point at epsilon = 0.05 ##########
###########################################################

step <- 1e2
loop <- 10

e <- .05
min_ess <- minESS(p, eps = e)
N_min <- round(min_ess/2)
all_ESS_H1 <- list()
all_ESS_H2 <- list()
all_ESS_H1[[1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, fun=fun1, h=0, p=p_H1)
all_ESS_H2[[1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, fun=fun2, h=0, p=p_H2)

for (i in 1:p_H1){
  min_ess <- minESS(1, eps = e)
  N_min <- round(min_ess/2)
  all_ESS_H1[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, fun = fun1, h=i, p=p_H1)
}

for (i in 1:p_H2){
  min_ess <- minESS(1, eps = e)
  N_min <- round(min_ess/2)
  all_ESS_H2[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, fun = fun2, h=i, p=p_H2)
}
save(all_ESS_H1, all_ESS_H2, file = "ESS_eps_H1_H2.Rdata")


