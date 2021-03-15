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
sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
lambda <- .6
Lambda <- diag(sigma)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
    Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
  }
}

########### Fixing proposal covariance matrix ############

tau <- c(rep(2,ceiling(p/2)), rep(((p+1)/p), floor(p/2)))
rho <- .5
Sigma <- diag(tau)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Sigma[i,j] <- rho*sqrt(tau[i]*tau[j])
    Sigma[j,i] <- rho*sqrt(tau[i]*tau[j])
  }
}

########### Validity check ##################
foo <- 2*solve(Lambda) - solve(Sigma)
min(eigen(foo)$values)
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

step <- 1e2
loop <- 20


###########################################################
############ Don't run, load ##############################
###########################################################

### Epsilon = 0.1
min_ess <- minESS(p, eps = 0.1)
N_min <- round(min_ess/2)
all_ESS1 <- list()
all_ESS1[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.1)
  all_ESS1[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS1, file = paste("ESS_eps1_p", p, ".Rdata", sep = ""))

### Epsilon = 0.07
min_ess <- minESS(p, eps = 0.08)
N_min <- round(min_ess/2)
all_ESS2[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.08)
  all_ESS2[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS2, file = paste("ESS_eps2_p", p, ".Rdata", sep = ""))



### Epsilon = 0.04
min_ess <- minESS(p, eps = 0.06)
N_min <- round(min_ess/2)
all_ESS3[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.06)
  all_ESS3[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS3, file = paste("ESS_eps3_p", p, ".Rdata", sep = ""))

### Epsilon = 0.01
min_ess <- minESS(p, eps = 0.04)
N_min <- round(min_ess/2)
all_ESS4[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.04)
  all_ESS4[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS4, file = paste("ESS_eps4_p", p, ".Rdata", sep = ""))

### Epsilon = 0.01
min_ess <- minESS(p, eps = 0.03)
N_min <- round(min_ess/2)
all_ESS5[[1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.03)
  all_ESS5[[i+1]] <- is_ESS(min_ess, step, loop, N_min, Sigma, Lambda, h=i, p=p)
}
save(all_ESS5, file = paste("ESS_eps5_p", p, ".Rdata", sep = ""))

#################################################################
#################################################################

for (i in 1:5)
  load(file = paste("ESS_eps", i, "_p", p, ".Rdata", sep = ""))

SE1 <- matrix(0, nrow = p+1, ncol = loop)
for (i in 1:(p+1))
  SE1[i,] <- apply(as.matrix(all_ESS1[[i]]$emp[,3], row = loop), 1, norm, '2')
SE2 <- matrix(0, nrow = p+1, ncol = loop)
for (i in 1:(p+1))
  SE2[i,] <- apply(as.matrix(all_ESS2[[i]]$emp[,3], row = loop), 1, norm, '2')
SE3 <- matrix(0, nrow = p+1, ncol = loop)
for (i in 1:(p+1))
  SE3[i,] <- apply(as.matrix(all_ESS3[[i]]$emp[,3], row = loop), 1, norm, '2')
SE4 <- matrix(0, nrow = p+1, ncol = loop)
for (i in 1:(p+1))
  SE4[i,] <- apply(as.matrix(all_ESS4[[i]]$emp[,3], row = loop), 1, norm, '2')
SE5 <- matrix(0, nrow = p+1, ncol = loop)
for (i in 1:(p+1))
  SE5[i,] <- apply(as.matrix(all_ESS5[[i]]$emp[,3], row = loop), 1, norm, '2')


min_ess_eps <- c(minESS(p, eps = .1), minESS(p, eps=.08), minESS(p, eps=.06), minESS(p, eps=.04), minESS(p, eps=.03))
truth <- min_ess_eps/((det(Lambda)/det(var.is(Sigma, Lambda)))^(1/p))
  
  
pdf(file = paste("multivariateNormal-multirep_p", p, ".pdf", sep = ""))

plot(all_ESS1[[1]]$emp[,1], SE1[1,], pch=19, col = "blue", xlim = range(truth[1], truth[5], all_ESS1[[1]]$emp[,1], all_ESS5[[1]]$emp[,1], all_ESS1[[2]]$emp[,1], all_ESS5[[2]]$emp[,1]), ylim = range(SE1, SE2, SE3, SE4), xlab = "Iterations", ylab = "Squared Error")
points(all_ESS2[[1]]$emp[,1], SE2[1,], pch=19, col = "green")
points(all_ESS3[[1]]$emp[,1], SE3[1,], pch=19, col = "orange")
points(all_ESS4[[1]]$emp[,1], SE4[1,], pch=19, col = "pink")
points(all_ESS5[[1]]$emp[,1], SE5[1,], pch=19, col = "red")

points(all_ESS1[[2]]$emp[,1], SE1[2,], pch = 1, col = adjustcolor("blue", alpha.f = 0.3))
points(all_ESS2[[2]]$emp[,1], SE2[2,], pch = 1, col = adjustcolor("green", alpha.f = 0.3))
points(all_ESS3[[2]]$emp[,1], SE3[2,], pch = 1, col = adjustcolor("orange", alpha.f = 0.3))
points(all_ESS4[[2]]$emp[,1], SE4[2,], pch = 1, col = adjustcolor("pink", alpha.f = 0.3))
points(all_ESS5[[2]]$emp[,1], SE5[2,], pch = 1, col = adjustcolor("red", alpha.f = 0.3))

abline(v = truth[1], lty=2, col = "blue")
abline(v = truth[2], lty=2, col = "green")
abline(v = truth[3], lty=2, col = "orange")
abline(v = truth[4], lty=2, col = "pink")
abline(v = truth[5], lty=2, col = "red")
legend("topright", legend = c("epsilon = 0.1", "epsilon = 0.08", "epsilon = 0.06", "epsilon = 0.04", "epsilon = 0.03"), col = c("blue", "green", "orange", "pink", "red"), pch = 19)
dev.off()


pdf(file = paste("multivariateNormal-multirep_estimates_p", p, ".pdf", sep = ""))

plot(all_ESS1[[1]]$emp[,1], all_ESS1[[1]]$emp[,3], pch=19, col = "blue", 
     xlim = range(truth[1], truth[5], all_ESS1[[1]]$emp[,1], all_ESS5[[1]]$emp[,1], all_ESS1[[2]]$emp[,1], all_ESS5[[2]]$emp[,1]), ylim = range(all_ESS1[[1]]$emp[,3], all_ESS5[[1]]$emp[,3], all_ESS1[[2]]$emp[,3], all_ESS1[[2]]$emp[,3]), xlab = "Iterations", ylab = "Squared Error")
points(all_ESS2[[1]]$emp[,1], all_ESS2[[1]]$emp[,3], pch=19, col = "green")
points(all_ESS3[[1]]$emp[,1], all_ESS3[[1]]$emp[,3], pch=19, col = "orange")
points(all_ESS4[[1]]$emp[,1], all_ESS4[[1]]$emp[,3], pch=19, col = "pink")
points(all_ESS5[[1]]$emp[,1], all_ESS5[[1]]$emp[,3], pch=19, col = "red")

points(all_ESS1[[2]]$emp[,1], all_ESS1[[2]]$emp[,3], pch = 1, col = adjustcolor("blue", alpha.f = 0.3))
points(all_ESS2[[2]]$emp[,1], all_ESS2[[2]]$emp[,3], pch = 1, col = adjustcolor("green", alpha.f = 0.3))
points(all_ESS3[[2]]$emp[,1], all_ESS3[[2]]$emp[,3], pch = 1, col = adjustcolor("orange", alpha.f = 0.3))
points(all_ESS4[[2]]$emp[,1], all_ESS4[[2]]$emp[,3], pch = 1, col = adjustcolor("pink", alpha.f = 0.3))
points(all_ESS5[[2]]$emp[,1], all_ESS5[[2]]$emp[,3], pch = 1, col = adjustcolor("red", alpha.f = 0.3))

abline(v = truth[1], lty=2, col = "blue")
abline(v = truth[2], lty=2, col = "green")
abline(v = truth[3], lty=2, col = "orange")
abline(v = truth[4], lty=2, col = "pink")
abline(v = truth[5], lty=2, col = "red")
abline(h = 0, lty=2, col = "black")
legend("topright", legend = c("epsilon = 0.1", "epsilon = 0.08", "epsilon = 0.06", "epsilon = 0.04", "epsilon = 0.03"), col = c("blue", "green", "orange", "pink", "red"), pch = 19)
dev.off()


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

p <- 10
sigma <- c(rep(2,ceiling(p/2)), rep(1, floor(p/2)))
lambda <- .6
Lambda <- diag(sigma)
for (i in 1:(p-1)){
  for (j in (i+1):p){
    Lambda[i,j] <- lambda*sqrt(sigma[i]*sigma[j])
    Lambda[j,i] <- lambda*sqrt(sigma[i]*sigma[j])
  }
}


########### Fixing proposal covariance matrix ############

tau <- c(rep(2,ceiling(p/2)), rep(((p+1)/p), floor(p/2)))
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

rho <- seq(.3, .7, .01)
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
  true.ess[1,t] <- (det(var.f)/det(true.var))^(1/p)
  for (i in 1:p){
    true.ess[i+1,t] <- var.f[i,i]/true.var[i,i]  
  }
}

pdf(file = paste("multiNormal-trueESSvsRho_det_p", p, ".pdf", sep = ""), width = 5, height = 5)

plot(rho, true.ess[1,], type = "l", lwd=2, ylim = range(true.ess), ylab = "True ESS")
for (i in 1:p)
  lines(rho, true.ess[i+1,], col = adjustcolor("blue", alpha.f = 0.3), lwd=2)
legend("topleft", legend = c("Multivariate", "Univariates"), col = c("black", adjustcolor("blue", alpha.f = 0.3)), lwd=2)
dev.off()
