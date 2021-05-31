set.seed(1)
rm(list = ls())
library(mvtnorm)
library(mcmcse)

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

is_ESS <- function(min_ESS, step, loop, N_min, Sigma, Lambda, fun, h=0, p){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  iter.true <- rep(N_min, loop)
  
  if (h==0){
    means.emp <- matrix(0, nrow = loop, ncol = p)
    means.kong <- matrix(0, nrow = loop, ncol = p)
    means.true <- matrix(0, nrow = loop, ncol = p)
  } else{
    means.emp <- matrix(0, nrow = loop, ncol = 1)
    means.kong <- matrix(0, nrow = loop, ncol = 1)
    means.true <- matrix(0, nrow = loop, ncol = 1)
  }
  ess.emp <- rep(0, loop)
  ess.kong <- rep(0, loop)
  ess.true <- rep(0, loop)
  
  for(t in 1:loop){
    
    print(t)
    is <- rmvnorm(N_min, mean=rep(0,p), sigma = Sigma )
    weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
    is <- fun(is)
    is.more <- is
    norm_weights <- weights/sum(weights)
    
    if (h==0)
      H <- as.matrix(is, ncol = p) 
    for (i in 1:p){
      if(h==i)
        H <- as.matrix(is[,i], ncol = 1)
    }
    
    snis <- colSums(H * norm_weights)
    
    true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
    var.f <- Lambda  
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[t] <- N_min*((det(varIbar)/det(emp.var))^(1/p))
    ess.kong[t]  <- 1/sum(norm_weights^2)  ##comment to not calculate kong
    ess.true[t] <- N_min*((det(var.f)/det(true.var))^(1/p))
    means.emp[t,] <- snis
    means.kong[t,] <- snis   ##comment to not calculate kong
    means.true[t,] <- snis
    
    while (ess.emp[t] <= min_ess ||ess.kong[t] <= min_ess || ess.true[t] <= min_ess){
      if(ess.emp[t]%%step==0) print(ess.emp[t])
      samp <- rmvnorm(step, mean=rep(0,p), sigma = Sigma )
      weights.more <- dmvnorm(samp, mean=rep(0,p), sigma=Lambda)/dmvnorm(samp, mean=rep(0,p), sigma = Sigma)   
      weights <- c(weights, weights.more)
      is.more <- rbind(is.more, fun(samp))
      if (h==0)
        H <- as.matrix(is.more, ncol = p) 
      for (i in 1:p){
        if(h==i)
          H <- as.matrix(is.more[,i], ncol = 1)
      }
      norm_weights <- weights/sum(weights)
      snis <- colSums(H * norm_weights)
      
      if (ess.emp[t] <= min_ess){
        iter.emp[t] <- iter.emp[t] + step
        varIbar <- (t(H - snis)%*%(norm_weights*(H-snis)))/iter.emp[t]
        emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
        ess.emp[t] <- iter.emp[t]*((det(varIbar)/det(emp.var))^(1/p))
        means.emp[t,] <- snis
      }
      
      if (ess.kong[t] <= min_ess){
        iter.kong[t] <- iter.kong[t] + step
        ess.kong[t] <- 1/sum(norm_weights^2)
        means.kong[t,] <- snis
      }
      
      if (ess.true[t] <= min_ess){
        iter.true[t] <- iter.true[t] + step
        ess.true[t] <- iter.true[t]*(det(var.f)/det(true.var))^(1/p)
        means.true[t,] <- snis
      }
      
    }
  }
  return(list("kong" = cbind(iter.kong, ess.kong, means.kong), 
              "emp" = cbind(iter.emp, ess.emp, means.emp), "true" = cbind(iter.true, ess.true, means.true)))
}

fun <- function(x) x