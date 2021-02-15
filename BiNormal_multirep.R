set.seed(1)
library(mvtnorm)
library(ellipse)
library(mcmcse)

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

lambda <- .5
rho <- .5
sigma <- (p+1)/p
Lambda <- matrix(c(2, lambda*sqrt(2), lambda*sqrt(2), 1), 2, 2)
Sigma <- matrix(c(sigma, rho*sigma, rho*sigma, sigma), 2, 2)
plot(ellipse(Lambda) , type = "l", xlim = c(-5,5), ylim = c(-5,5))
lines(ellipse(Sigma), col = "blue")

is_ESS <- function(min_ESS, step, loop, N_min, Sigma, Lambda){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  iter.true <- rep(N_min, loop)
  means.emp <- matrix(0, nrow = loop, ncol = 2)
  means.kong <- matrix(0, nrow = loop, ncol = 2)
  means.true <- matrix(0, nrow = loop, ncol = 2)
  ess.emp <- rep(0, loop)
  ess.kong <- rep(0, loop)
  ess.true <- rep(0, loop)
  
  for(t in 1:loop){
    
    print(t)
    is <- rmvnorm(N_min, mean=rep(0,2), sigma = Sigma )
    is.more <- is
    weights <- dmvnorm(is, mean=rep(0,2), sigma=Lambda)/dmvnorm(is, mean=rep(0,2), sigma = Sigma)   
    norm_weights <- weights/sum(weights)
    H <- is
    snis <- colSums(H * norm_weights)
    
    true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
    var.f <- Lambda  
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[t] <- N_min*(det(varIbar)/det(emp.var))
    ess.kong[t]  <- 1/sum(norm_weights^2)
    ess.true[t] <- N_min*(det(var.f)/det(true.var))
    means.emp[t,] <- snis
    means.kong[t,] <- snis
    means.true[t,] <- snis
    
    while (ess.emp[t] <= min_ess ||ess.kong[t] <= min_ess || ess.true[t] <= min_ess){
      
      samp <- rmvnorm(step, mean=rep(0,2), sigma = Sigma)
      is.more <- rbind(is.more, samp)
      H <- is.more
      weights <- dmvnorm(is.more, mean=rep(0,2), sigma=Lambda)/dmvnorm(is.more, mean=rep(0,2), sigma = Sigma)   
      norm_weights <- weights/sum(weights)
      snis <- colSums(H * norm_weights)
      
      if (ess.emp[t] <= min_ess){
        iter.emp[t] <- iter.emp[t] + step
        varIbar <- (t(H - snis)%*%(norm_weights*(H-snis)))/iter.emp[t]
        emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
        ess.emp[t] <- iter.emp[t]*(det(varIbar)/det(emp.var))
        means.emp[t,] <- snis
      }
      
      if (ess.kong[t] <= min_ess){
        iter.kong[t] <- iter.kong[t] + step
        ess.kong[t] <- 1/sum(norm_weights^2)
        means.kong[t,] <- snis
      }
      
      if (ess.true[t] <= min_ess){
        iter.true[t] <- iter.true[t] + step
        ess.true[t] <- iter.true[t]*det(var.f)/det(true.var)
        means.true[t,] <- snis
      }
      
    }
  }
  return(list("kong" = cbind(iter.kong, ess.kong, means.kong), 
              "emp" = cbind(iter.emp, ess.emp, means.emp), "true" = cbind(iter.true, ess.true, means.true)))
}

min_ess <- minESS(2)
step <- 100
N_min <- 1e3
loop <- 100

all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda)

pdf(file = "biNorm-multirep.pdf", height = 5, width = 10)
par(mfrow = c(1,2))

plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = I")
points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)

plot(all_ESS$emp[,1], all_ESS$emp[,4], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,4], all_ESS$emp[,4], all_ESS$true[,4]), xlab = "Iterations", ylab = "Component = II")
points(all_ESS$kong[,1], all_ESS$kong[,4], pch = 19, col = "orange")
points(all_ESS$true[,1], all_ESS$true[,4], pch = 19, col = "green")
legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)

dev.off()