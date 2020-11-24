
############ Visualizing target and proposal distribution #############
x <- seq(-10, 10, .1)
plot(x, dnorm(x, mean = -2, sd=1), "l")
lines(x, .2*dnorm(x, mean = -4, sd=1) + .8*dnorm(x, mean=4, sd =1), col = "red")
#######################################################################

alpha <- seq(.5,3,.1)
ess.boot <- alpha
ess.boot.trunc <- alpha
for (k in 1:length(alpha)){
  N <- 1000       ## no. of importance samples
  X <- rep(0, N)
  for (i in 1:N){
    if (runif(1) <= 1)
      X[i] = rnorm(1, mean = 0, sd  = 1)
    else X[i] = rnorm(1, mean = 3, sd = 1)
  }
  W = dnorm(X, mean = 0, sd = 1)/dnorm(X, mean = 0, sd = 1)  ##importance weights
  H <- rep(0,N)
  for (j in 1:N)
    if(abs(X[j]) > alpha[k]) H[j] <- 1
  
  snis <- sum(W*H)/sum(W)    ## Wighted IS estimator
  iid.var <- (cov.wt(as.matrix(H, ncol  =1), wt = W/sum(W))$cov)/N    ## varinace of IID Monte Carlo estimator, Lambda
  
  
  ################################################
  ########## Non-truncated Bootstrap #############
  ################################################
  
  B <- N          ## no. of bootstrap steps       
  boot_means <- rep(0, B)
  for (i in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    boot_means[i] <- sum(H[sample_idx] * W[sample_idx])/sum(W[sample_idx])
  }
  ess.boot[k] <- N * iid.var/(var(boot_means))   ## non-truncated bootstrap estimator for variance 
  print(ess.boot[k])
  
  ###############################################
  ######### Truncated Bootstrap ################
  ###############################################
  
  Delta <- rep(0, B)  
  tau <- max(.04*abs(snis), .01)
  for (i in 1:B){
    if(boot_means[i] - snis > tau)
    {
      #print("Truncated")
      Delta[i] = tau
    }
    else if(abs(boot_means[i] - snis) <= tau)
      Delta[i] = boot_means[i] - snis
    else {
      #print("Truncated")
      Delta[i] = -tau
    }
    
  }
  ess.boot.trunc[k] <- N * iid.var/var(Delta)
  print(ess.boot.trunc[k])
}

plot(alpha, ess.boot.trunc, "l")
lines(alpha, ess.boot, col = "red")
legend("topright", legend = c("Trunc. Bootstrap", "Bootstrap"), col = c("black", "red"), lty=1)
