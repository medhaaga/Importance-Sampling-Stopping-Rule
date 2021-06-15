set.seed(1)
rm(list = ls())
library(mvtnorm)
library(mcmcse)

var.snis <- function(Sigma, Lambda, mu=0){
  foo <- 2*solve(Lambda) - solve(Sigma)
  return((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo)))) * (solve(foo)))
}

var.uis <- function(Sigma, Lambda, mu=0)
{
  foo <- 2*solve(Lambda) - solve(Sigma)
  return(((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo)))) * (solve(foo) + outer(mu, mu)))  - outer(mu, mu))
}

obj.det <- function(par, Lambda){
  sigma <- par[1]
  rho <- par[2]
  Sigma <- matrix(c(sigma, rho*sigma, rho*sigma, sigma), 2, 2)
  foo <- 2*solve(Lambda) - solve(Sigma)
  return((sqrt(det(Sigma))/(det(Lambda) * sqrt(det(foo))))^2 / det(foo))
}

is_ESS <- function(min_ESS, step, loop, N_min, Sigma, Lambda, mu, fun, h=0, p){
  
  iter.emp.uis <- rep(N_min, loop)
  iter.emp.snis <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  iter.true <- rep(N_min, loop)
  
  if (h==0){
    means.emp.uis <- matrix(0, nrow = loop, ncol = p)
    means.emp.snis <- matrix(0, nrow = loop, ncol = p)
    #means.kong <- matrix(0, nrow = loop, ncol = p)
    #means.true <- matrix(0, nrow = loop, ncol = p)
  } else{
    means.emp.uis <- matrix(0, nrow = loop, ncol = 1)
    means.emp.snis <- matrix(0, nrow = loop, ncol = 1)
    #means.kong <- matrix(0, nrow = loop, ncol = 1)
    #means.true <- matrix(0, nrow = loop, ncol = 1)
  }
  ess.emp.uis <- rep(0, loop)
  ess.emp.snis <- rep(0, loop)
  #ess.kong <- rep(0, loop)
  #ess.true <- rep(0, loop)
  
  for(t in 1:loop){
    
    print(t)
    is <- rmvnorm(N_min, mean = mu, sigma = Sigma )
    weights <- dmvnorm(is, mean=mu, sigma=Lambda)/dmvnorm(is, mean=mu, sigma = Sigma)   
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
    uis <- colSums(H * weights)/N_min
    if(h==0){
      foo1 <- matrix(uis, nrow = iter.emp.uis[t], ncol = p, byrow = TRUE)
      foo2 <- matrix(snis, nrow = iter.emp.snis[t], ncol = p, byrow = TRUE)
    } else{
      foo1 <- matrix(uis, nrow = iter.emp.uis[t], ncol = 1, byrow = TRUE)
      foo2 <- matrix(snis, nrow = iter.emp.snis[t], ncol = 1, byrow = TRUE)
    }
    
    #true.var.uis <- var.uis(Sigma, Lambda)
    #true.var.snis <- var.snis(Sigma, Lambda)  ## N*true IS variance
    #var.f <- Lambda  
    varIbar.uis <- (t(H - foo1) %*% (weights * (H - foo1)))/(N_min^2)
    varIbar.snis <- (t(H - foo2) %*% (norm_weights * (H - foo2)))/N_min
    
    emp.var.uis <-  (t(weights*H - foo1) %*% (weights*H - foo1))/(iter.emp.uis[t]^2)
    emp.var.snis <- (t(norm_weights*(H - foo2)) %*% (norm_weights * (H - foo2)))
    ess.emp.uis[t] <- iter.emp.uis[t]*((det(varIbar.uis)/det(emp.var.uis))^(1/p))
    ess.emp.snis[t] <- iter.emp.snis[t]*((det(varIbar.snis)/det(emp.var.snis))^(1/p))
    
    #ess.kong[t]  <- 1/sum(norm_weights^2)  ##comment to not calculate kong
    #ess.true[t] <- N_min*((det(var.f)/det(true.var))^(1/p))
    means.emp.uis[t,] <- uis
    means.emp.snis[t,] <- snis
    #means.kong[t,] <- snis   ##comment to not calculate kong
    
    while (ess.emp.uis[t] <= min_ess || ess.emp.snis[t] <= min_ess){
      
      if(ess.emp.uis[t]%%step==0) print(ess.emp[t])
      
      samp <- rmvnorm(step, mean=mu, sigma = Sigma )
      weights.more <- dmvnorm(samp, mean=mu, sigma=Lambda)/dmvnorm(samp, mean=mu, sigma = Sigma)   
      weights <- c(weights, weights.more)
      is.more <- rbind(is.more, fun(samp))
      if (h==0)
        H <- as.matrix(is.more, ncol = p) 
      for (i in 1:p){
        if(h==i)
          H <- as.matrix(is.more[,i], ncol = 1)
      }
      norm_weights <- weights/sum(weights)
      uis <- colSums(H * weights)/length(weights)
      snis <- colSums(H * norm_weights)
      
      
      if (ess.emp.uis[t] <= min_ess){
        iter.emp.uis[t] <- iter.emp.uis[t] + step
        if(h==0)
          foo1 <- matrix(uis, nrow = iter.emp.uis[t], ncol = p, byrow = TRUE) else
            foo1 <- matrix(uis, nrow = iter.emp.uis[t], ncol = 1, byrow = TRUE)
        
        varIbar.uis <- (t(H - foo1)%*%(weights*(H-foo1)))/(iter.emp.uis[t]^2)
        emp.var.uis <-  (t(weights*H - foo1) %*% (weights*H - foo1))/(iter.emp.uis[t]^2)
        ess.emp.uis[t] <- iter.emp.uis[t]*((det(varIbar.uis)/det(emp.var.uis))^(1/p))
        means.emp.uis[t,] <- uis
      }
      
      if (ess.emp.snis[t] <= min_ess){
        iter.emp.snis[t] <- iter.emp.snis[t] + step
        if(h==0)
          foo2 <- matrix(snis, nrow = iter.emp.snis[t], ncol = p, byrow = TRUE) else
            foo2 <- matrix(snis, nrow = iter.emp.snis[t], ncol = 1, byrow = TRUE)
        varIbar.snis <- (t(H - foo2)%*%(norm_weights*(H-foo2)))/iter.emp.snis[t]
        emp.var.snis <- (t(H - foo2) %*% (norm_weights^2 * (H - foo2)))
        ess.emp.snis[t] <- iter.emp.snis[t]*((det(varIbar.snis)/det(emp.var.snis))^(1/p))
        means.emp.snis[t,] <- snis
      }
      
      # if (ess.kong[t] <= min_ess){
      #   iter.kong[t] <- iter.kong[t] + step
      #   ess.kong[t] <- 1/sum(norm_weights^2)
      #   means.kong[t,] <- snis
      # }
      # 
      # if (ess.true[t] <= min_ess){
      #   iter.true[t] <- iter.true[t] + step
      #   ess.true[t] <- iter.true[t]*(det(var.f)/det(true.var))^(1/p)
      #   means.true[t,] <- snis
      # }
      
    }
  }
  return(list("emp_uis" = cbind(iter.emp.uis, ess.emp.uis, means.emp.uis), "emp_snis" = cbind(iter.emp.snis, ess.emp.snis, means.emp.snis)))
}

fun <- function(x) x