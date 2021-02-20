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

# try another setting with low correlation 
# but inefficient proposal distribution
# and see if the results are not so much in our 
# favor any more. 
# saving this: lambda = .9, rho = .9, sigma = 3
lambda <- .1
rho <- .1
sigma <- 3#(p+1)/p
Lambda <- matrix(c(3, lambda*sqrt(3), lambda*sqrt(3), 1), 2, 2)
Sigma <- matrix(c(sigma, rho*sigma, rho*sigma, sigma), 2, 2)
foo <- 2*solve(Lambda) - solve(Sigma)
det(foo)

pdf(file = "ellipse_highCorr_highCorr.pdf", height = 5, width = 5)
plot(ellipse(Lambda) , type = "l", xlim = c(-8,8), ylim = c(-8,8))
lines(ellipse(Sigma), col = "blue")
lines(ellipse(var.is(Sigma, Lambda)), col = "red")
dev.off()


is_ESS <- function(min_ESS, step, loop, N_min, Sigma, Lambda, h=0, p){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  iter.true <- rep(N_min, loop)
  if (h==0){
    means.emp <- matrix(0, nrow = loop, ncol = 2)
    means.kong <- matrix(0, nrow = loop, ncol = 2)
    means.true <- matrix(0, nrow = loop, ncol = 2)
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
    is <- rmvnorm(N_min, mean=rep(0,2), sigma = Sigma )
    is.more <- is
    weights <- dmvnorm(is, mean=rep(0,2), sigma=Lambda)/dmvnorm(is, mean=rep(0,2), sigma = Sigma)   
    norm_weights <- weights/sum(weights)
    if (h==0)
      H <- as.matrix(is, ncol = 2) else if(h==1)
        H <- as.matrix(is[,1], ncol = 1) else
          H <- as.matrix(is[,2], ncol = 1)
    
    
    snis <- colSums(H * norm_weights)
    
    true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
    var.f <- Lambda  
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[t] <- N_min*((det(varIbar)/det(emp.var))^(1/p))
    ess.kong[t]  <- 1/sum(norm_weights^2)
    ess.true[t] <- N_min*((det(var.f)/det(true.var))^(1/p))
    means.emp[t,] <- snis
    means.kong[t,] <- snis
    means.true[t,] <- snis
    
    while (ess.emp[t] <= min_ess ||ess.kong[t] <= min_ess || ess.true[t] <= min_ess){
      
      samp <- rmvnorm(step, mean=rep(0,2), sigma = Sigma)
      is.more <- rbind(is.more, samp)
      if (h==0)
        H <- as.matrix(is.more, ncol = 2) else if(h==1)
          H <- as.matrix(is.more[,1], ncol = 1) else
            H <- as.matrix(is.more[,2], ncol = 1)
      weights <- dmvnorm(is.more, mean=rep(0,2), sigma=Lambda)/dmvnorm(is.more, mean=rep(0,2), sigma = Sigma)   
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
        ess.true[t] <- iter.true[t]*(det(var.f)/det(true.var)^(1/p))
        means.true[t,] <- snis
      }
      
    }
  }
  return(list("kong" = cbind(iter.kong, ess.kong, means.kong), 
              "emp" = cbind(iter.emp, ess.emp, means.emp), "true" = cbind(iter.true, ess.true, means.true)))
}

min_ess <- minESS(2)
step <- 100
N_min <- 5e3
loop <- 20

all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=0, p=2)

pdf(file = "biNorm-multirep-bivariate_lowCorr.pdf", height = 5, width = 10)
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

min_ess <- minESS(1)
step <- 10
N_min <- 3e3
loop <- 20


pdf(file = "biNorm-multirep-univariate_lowCorr.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=1, p=2)

plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = I")
points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)

all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=2, p=2)

plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = II")
points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)
dev.off()


#####################################################
###### Univarate and bivariate true ESS vs rho ######
#####################################################

lambda <- .1
sigma <- 3#(p+1)/p
rho <- seq(0, .65, .05)

true.ess <- matrix(0, nrow = 3, ncol = length(rho))
for (i in 1:length(rho)){
  corr <- rho[i]
  Lambda <- matrix(c(3, lambda*sqrt(3), lambda*sqrt(3), 1), 2, 2)
  Sigma <- matrix(c(sigma, corr*sigma, corr*sigma, sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  
  true.ess[1,i] <- (det(var.f)/det(true.var))^(1/p)
  true.ess[2,i] <- var.f[1,1]/true.var[1,1]
  true.ess[3,i] <- var.f[2,2]/true.var[2,2]
}

pdf(file = "trueESSvsRho_lowCorr_det.pdf", width = 5, height = 5)
plot(rho, true.ess[1,], type = "l", lwd=2, ylim = range(true.ess), ylab = "True ESS")
lines(rho, true.ess[2,], col = "blue", lwd=2)
lines(rho, true.ess[3,], col = "green", lwd=2)
legend("bottomleft", legend = c("Bivariate", "Comp-1", "Comp2"), col = c("black", "blue", "green"), lwd=2)
dev.off()

lambda <- .9
sigma <- 3#(p+1)/p
rho <- seq(0, .9, .05)

true.ess <- matrix(0, nrow = 3, ncol = length(rho))
for (i in 1:length(rho)){
  corr <- rho[i]
  Lambda <- matrix(c(3, lambda*sqrt(3), lambda*sqrt(3), 1), 2, 2)
  Sigma <- matrix(c(sigma, corr*sigma, corr*sigma, sigma), 2, 2)
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  
  true.ess[1,i] <- det((var.f))^.5/det((true.var))^.5
  true.ess[2,i] <- var.f[1,1]/true.var[1,1]
  true.ess[3,i] <- var.f[2,2]/true.var[2,2]
}

pdf(file = "trueESSvsRho_highCorr_det.pdf", width = 5, height = 5)
plot(rho, true.ess[1,], type = "l", lwd=2, ylim = range(true.ess), ylab = "True ESS")
lines(rho, true.ess[2,], col = "blue", lwd=2)
lines(rho, true.ess[3,], col = "green", lwd=2)
legend("topleft", legend = c("Bivariate", "Comp-1", "Comp2"), col = c("black", "blue", "green"), lwd=2)
dev.off()


