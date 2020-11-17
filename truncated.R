
############ Visualizing target and proposal distribution #############
x <- seq(-10, 10, .1)
plot(x, dnorm(x, mean = -2, sd=1), "l")
lines(x, .5*dnorm(x, mean = -4, sd=1) + .5*dnorm(x, mean=4, sd =1), col = "red")
#######################################################################


N <- 5000       ## no. of importance samples
X <- rep(0, N)
for (i in 1:N){
  if (runif(1) <= .5)
    X[i] = rnorm(1, mean = -3, sd  = 1)
  else X[i] = rnorm(1, mean = 3, sd = 1)
}
W = dnorm(X, mean = -2, sd = 1)/(.5*dnorm(X, mean = -3, sd=1) + .5*dnorm(X, mean=3, sd =1))  ##importance weights
snis <- sum(W*X)/sum(W)    ## Wighted IS estimator
iid.var <- (cov.wt(as.matrix(X, ncol  =1), wt = W/sum(W))$cov)/N    ## varinace of IID Monte Carlo estimator, Lambda


################################################
########## Non-truncated Bootstrap #############
################################################

B <- N          ## no. of bootstrap steps       
boot_means <- rep(0, B)
for (i in 1:B){
  sample_idx <- sample(1:N, size= N, replace = TRUE)
  boot_means[i] <- sum(X[sample_idx] * W[sample_idx])/sum(W[sample_idx])
}
ess.boot <- N * iid.var/(var(boot_means))   ## non-truncated bootstrap estimator for variance 
ess.boot

###############################################
######### Triuncated Bootstrap ################
###############################################

Delta <- rep(0, B)  
tau <- max(.04*abs(snis), .01)
for (i in 1:B){
  if(boot_means[i] - snis > tau)
  {
    print("Truncated")
    Delta[i] = tau
  }
  else if(abs(boot_means[i] - snis) <= tau)
    Delta[i] = boot_means[i] - snis
  else {
    print("Truncated")
    Delta[i] = -tau
  }
    
}
ess.boot.trunc <- N * iid.var/var(Delta)
ess.boot.trunc
