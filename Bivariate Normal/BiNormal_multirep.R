set.seed(1)
library(mvtnorm)
library(ellipse)
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

##### Fixing target covariance matrix #################

p <- 3
sigma <- 2^seq(1:p)
lambda <- .5
Lambda <- diag(sigma)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
    Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
  }
}

########### Fixing proposal covariance matrix ############

tau <- ((p+1)/p)^seq(1,p)
rho <- .7
Sigma <- diag(tau)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Sigma[i,j] <- rho*sqrt(tau[i]*tau[j])
    Sigma[j,i] <- rho*sqrt(tau[i]*tau[j])
  }
}

########### Validity check ##################
foo <- 2*solve(Lambda) - solve(Sigma)
det(foo)
true <- var.is(Sigma, Lambda)
(det(Lambda)/det(true))^(1/p)
############## print ellipse when p=2 ##################
#pdf(file = "ellipse_optimal.pdf", height = 5, width = 5)
#plot(ellipse(Lambda) , type = "l", xlim = c(-4,4), ylim = c(-4,4))
#lines(ellipse(Sigma), col = "blue")
#lines(ellipse(var.is(Sigma, Lambda)), col = "red")
#dev.off()
########################################################

is_ESS <- function(min_ESS, step, loop, N_min, Sigma, Lambda, h=0, p){
  
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
    is.more <- is
    weights <- dmvnorm(is, mean=rep(0,p), sigma=Lambda)/dmvnorm(is, mean=rep(0,p), sigma = Sigma)   
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
      is.more <- rbind(is.more, samp)
      if (h==0)
        H <- as.matrix(is.more, ncol = p) 
      for (i in 1:p){
        if(h==i)
          H <- as.matrix(is.more[,i], ncol = 1)
      }
      weights <- dmvnorm(is.more, mean=rep(0,p), sigma=Lambda)/dmvnorm(is.more, mean=rep(0,p), sigma = Sigma)   
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

step <- 10^p
loop <- 20

### Epsilon = 0.1
min_ess <- minESS(p, eps = 0.1)
N_min <- round(min_ess/2)
all_ESS1 <- list()
for (i in 0:p){
  all_ESS1[[1]] <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS1, file = "ESS_eps1.Rdata")

### Epsilon = 0.07
min_ess <- minESS(p, eps = 0.07)
N_min <- round(min_ess/2)
all_ESS2 <- list()
for (i in 0:p){
  all_ESS2[[1]] <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS1, file = "ESS_eps2.Rdata")

### Epsilon = 0.05
min_ess <- minESS(p, eps = 0.05)
N_min <- round(min_ess/2)
all_ESS3 <- list()
for (i in 0:p){
  all_ESS3[[1]] <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS1, file = "ESS_eps3.Rdata")

### Epsilon = 0.02
min_ess <- minESS(p, eps = 0.02)
N_min <- round(min_ess/2)
all_ESS4 <- list()
for (i in 0:p){
  all_ESS4[[1]] <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS1, file = "ESS_eps4.Rdata")

SE1 <- apply(all_ESS1[[1]]$emp[,3:(p+2)], 1, norm, '2')
SE2 <- apply(all_ESS2[[1]]$emp[,3:(p+2)], 1, norm, '2')
SE3 <- apply(all_ESS3[[1]]$emp[,3:(p+2)], 1, norm, '2')
SE4 <- apply(all_ESS4[[1]]$emp[,3:(p+2)], 1, norm, '2')
plot(all_ESS1[[1]]$emp[,1], SE1, pch=19, col = "blue", xlim = range(all_ESS1[[1]]$emp[,1], all_ESS2[[1]]$emp[,1], all_ESS3[[1]]$emp[,1], all_ESS4[[1]]$emp[,1],), ylim = range(SE1, SE2, SE3, SE4), xlab = "Iterations", ylab = "Component = I")
points(all_ESS2[[1]]$emp[,1], SE2, pch=19, col = "blue")
points(all_ESS3[[1]]$emp[,1], SE3, pch=19, col = "blue")
points(all_ESS4[[1]]$emp[,1], SE4, pch=19, col = "blue")


#pdf(file = "biNorm-multirep-bivariate_setting3.pdf", height = 5, width = 10)
#par(mfrow = c(1,2))

#plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = I")
#points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
#points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
#legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)

#plot(all_ESS$emp[,1], all_ESS$emp[,4], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,4], all_ESS$emp[,4], all_ESS$true[,4]), xlab = "Iterations", ylab = "Component = II")
#points(all_ESS$kong[,1], all_ESS$kong[,4], pch = 19, col = "orange")
#points(all_ESS$true[,1], all_ESS$true[,4], pch = 19, col = "green")
#legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)
#dev.off()


#min_ess <- minESS(1)
#step <- 10
#N_min <- 3e3
#loop <- 20


#pdf(file = "biNorm-multirep-univariate_setting1.pdf", height = 5, width = 10)
#par(mfrow = c(1,2))
#all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=1, p=2)

#plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = I")
#points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
#points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
#legend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)

#all_ESS <- is_ESS(min_ESS, step, loop, N_min, Sigma, Lambda, h=2, p=2)

#plot(all_ESS$emp[,1], all_ESS$emp[,3], pch=19, col = "blue", xlim = range(all_ESS$kong[,1], all_ESS$emp[,1], all_ESS$true[,1]), ylim = range(all_ESS$kong[,3], all_ESS$emp[,3], all_ESS$true[,3]), xlab = "Iterations", ylab = "Component = II")
#points(all_ESS$kong[,1], all_ESS$kong[,3], pch = 19, col = "orange")
#points(all_ESS$true[,1], all_ESS$true[,3], pch = 19, col = "green")
#egend("topright", legend = c("Truth", "Empirical", "Kong"), col = c("green", "blue", "orange"), pch=19)
#dev.off()


#####################################################
###### Univarate and bivariate true ESS vs rho ######
#####################################################

##### Fixing target covariance matrix #################

p <- 3
sigma <- 2^seq(1:p)
lambda <- .5
Lambda <- diag(sigma)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
    Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
  }
}

########### Fixing proposal covariance matrix ############

tau <- ((p+1)/p)^seq(1,p)
rho <- .3
Sigma <- diag(tau)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Sigma[i,j] <- rho*sqrt(tau[i]*tau[j])
    Sigma[j,i] <- rho*sqrt(tau[i]*tau[j])
  }
}

########### Validity check ##################
foo <- 2*solve(Lambda) - solve(Sigma)
det(foo)

rho <- seq(.4, .9, .01)
true.ess <- matrix(0, nrow = p+1, ncol = length(rho))
for (t in 1:length(rho)){
  corr <- rho[t]
  Sigma <- diag(tau)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      Sigma[i,j] <- corr*sqrt(tau[i]*tau[j])
      Sigma[j,i] <- corr*sqrt(tau[i]*tau[j])
    }
  }
  true.var <- var.is(Sigma, Lambda)  ## N*true IS variance
  var.f <- Lambda  
  true.ess[1,i] <- (det(var.f)/det(true.var))^(1/p)
  for (i in 1:p){
    true.ess[i+1,t] <- var.f[i,i]/true.var[i,i]  
  }
}

pdf(file = "trueESSvsRho_det.pdf", width = 5, height = 5)
plot(rho, true.ess[1,], type = "l", lwd=2, ylim = range(true.ess), ylab = "True ESS")
for (i in 1:p)
  lines(rho, true.ess[i+1,], col = adjustcolor("blue", alpha.f = 0.3), lwd=2)
legend("bottomleft", legend = c("Bivariate", "Univariates"), col = c("blue", adjustcolor("blue", alpha.f = 0.3)), lwd=2)
dev.off()
