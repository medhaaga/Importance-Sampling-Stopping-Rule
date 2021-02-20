library(EnvStats)
library(mcmcse)
proj.hour = c(
  387,182,244,600,627,
  332,418,300,798,584,
  660,39,274,174,50,
  34,1895,158,974,345,
  1755,1752,473,81,954,
  1407,230,464,380,131,
  1205
)

post <- function(lambda, beta, t)
{
  if( lambda < 0 || beta < 0)
  {
    foo <- 0
  } else{
  #foo <- lambda^(32.5)*beta^(31) * (prod(t))^(beta - 1) + exp(-lambda * sum(t^beta)) * exp(-beta) * exp(-2350*lambda)
  foo <- 32.5*log(lambda) + 31*log(beta) + (beta - 1)*sum(log(t)) - lambda*sum(t^beta) - beta - 2350*lambda
  foo <- exp(foo)
  }
  return(foo)
}

joint_draw <- function(N, t, h)
{
  mle.fit <- eweibull(proj.hour)
  beta.fit <- mle.fit$parameters[1]
  lambda.fit <- mle.fit$parameters[2]
  
  beta <- rnorm(N, mean = beta.fit, sd = h)
  
  lambda <- numeric(length = N)
  importance <- numeric(length = N)
  posterior <- numeric(length = N)
  
  for(i in 1:N)
  {
    lambda[i] <- rgamma(1, shape = 33.5, rate = 2350 + sum(t^beta[i]))
    
    importance[i] <- dnorm(beta[i], mean = beta.fit, sd = h)
    importance[i] <- importance[i]*dgamma(lambda[i], shape = 33.5, rate = 2350 + sum(t^beta[i]))
    
    posterior[i] <- post(lambda[i], beta[i], t = proj.hour)
  }
  weights <- posterior/importance
  return(list(weights, cbind(lambda, beta)))
}

mttf <- function(lambda, beta)
{
  lambda^(-1/beta) *gamma(1 + 1/beta)
}

reliability <- function(lambda, beta)
{
  exp(-lambda*(1500^beta))
}

is_ESS <- function(min_ESS, step, loop, N_min, proj.hour, beta_sd, f=0, p=2){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  if (f==0){
    means.emp <- matrix(0, nrow = loop, ncol = 2)
    means.kong <- matrix(0, nrow = loop, ncol = 2)
  } else{
    means.emp <- matrix(0, nrow = loop, ncol = 1)
    means.kong <- matrix(0, nrow = loop, ncol = 1)
  }
  ess.emp <- rep(0, loop)
  ess.kong <- rep(0, loop)
  
  for(t in 1:loop){
    
    print(t)
    is <- joint_draw(N = N_min, t = proj.hour, h=beta_sd)
    is.more <- is
    L <- is[[2]][,1]
    B <- is[[2]][,2]
    weights <- is[[1]]
    norm_weights <- weights/sum(weights)
    if(f==0)
      H <- as.matrix(cbind(mttf(L,B), reliability(L,B)), ncol = 2) else if(f == 1)
        H <- as.matrix(mttf(L,B), ncol  = 1) else
          H <- as.matrix(reliability(L,B), ncol  =1)
    
    snis <- colSums(H * norm_weights)
    
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[t] <- N_min*((det(varIbar)/det(emp.var))^(1/p))
    ess.kong[t]  <- 1/sum(norm_weights^2)
    means.emp[t,] <- snis
    means.kong[t,] <- snis
    
    while (ess.emp[t] <= min_ess ||ess.kong[t] <= min_ess){
      
      samp <- joint_draw(N = step, t = proj.hour, h=beta_sd)
      is.more[[1]] <- c(is.more[[1]], samp[[1]])
      is.more[[2]] <- rbind(is.more[[2]], samp[[2]])
      L <- is.more[[2]][,1]
      B <- is.more[[2]][,2]
      if(f==0)
        H <- as.matrix(cbind(mttf(L,B), reliability(L,B)), ncol = 2) else if(f == 1)
          H <- as.matrix(mttf(L,B), ncol  = 1) else
            H <- as.matrix(reliability(L,B), ncol  =1)
      weights <- is.more[[1]]
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
      
    }
  }
  return(list("kong" = cbind(iter.kong, ess.kong, means.kong), "emp" = cbind(iter.emp, ess.emp, means.emp)))
}

step <- 10
loop <- 10

# beta_sd = .2 and .1
min_ess <- minESS(2)
N_min <- 5e3

all_h20_bi <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .2, f=0, p=2)
all_h10_bi <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .1, f=0, p=2)

min_ess <- minESS(1)
N_min <- 3e3

all_h20_uni1 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .2, f=1, p=2)
all_h20_uni2 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .1, f=1, p=2)
all_h10_uni1 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .2, f=2, p=2)
all_h10_uni2 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .1, f=2, p=2)

save(all_h10_bi, all_h20_bi, all_h20_uni1, all_h20_uni2, all_h10_uni1, all_h10_uni2, file = "weibull-multirep.Rdata")
load(file = "weibull-multirep.Rdata")


## Bivariate, sd = 0.1
pdf(file = "weibull-multirep_sd1_bi.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(all_h10_bi$emp[,1], all_h10_bi$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.1, MTTF vs Iterations", xlab = "Iterations", ylab = "Estimated MTTF", xlim = range(all_h10_bi$emp[,1], all_h10_bi$kong[,1]), ylim = range(all_h10_bi$emp[,3], all_h10_bi$kong[,3]))
points(all_h10_bi$kong[,1], all_h10_bi$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)

plot(all_h10_bi$emp[,1], all_h10_bi$emp[,4], pch=19, col = "blue", main = "Beta SD = 0.1, Reliability vs Iterations", xlab = "Iterations", ylab = "Estimated Reliability", xlim = range(all_h10_bi$emp[,1], all_h10_bi$kong[,1]), ylim = range(all_h10_bi$emp[,4], all_h10_bi$kong[,4]))
points(all_h10_bi$kong[,1], all_h10_bi$kong[,4], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)
dev.off()

## Univariate, sd = 0.1
pdf(file = "weibull-multirep_sd1_uni.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(all_h10_uni1$emp[,1], all_h10_uni1$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.1, MTTF vs Iterations", xlab = "Iterations", ylab = "Estimated MTTF", xlim = range(all_h10_uni1$emp[,1], all_h10_uni1$kong[,1]), ylim = range(all_h10_uni1$emp[,3], all_h10_uni1$kong[,3]))
points(all_h10_uni1$kong[,1], all_h10_uni1$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)

plot(all_h10_uni2$emp[,1], all_h10_uni2$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.1, Reliability vs Iterations", xlab = "Iterations", ylab = "Estimated Reliability", xlim = range(all_h10_uni2$emp[,1], all_h10_uni2$kong[,1]), ylim = range(all_h10_uni2$emp[,3], all_h10_uni2$kong[,3]))
points(all_h10_uni2$kong[,1], all_h10_uni2$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)
dev.off()

## Bivariate, sd = 0.2
pdf(file = "weibull-multirep_sd2_bi.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(all_h20_bi$emp[,1], all_h20_bi$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.2, MTTF vs Iterations", xlab = "Iterations", ylab = "Estimated MTTF", xlim = range(all_h20_bi$emp[,1], all_h20_bi$kong[,1]), ylim = range(all_h20_bi$emp[,3], all_h20_bi$kong[,3]))
points(all_h20_bi$kong[,1], all_h20_bi$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)

plot(all_h20_bi$emp[,1], all_h20_bi$emp[,4], pch=19, col = "blue", main = "Beta SD = 0.2, Reliability vs Iterations", xlab = "Iterations", ylab = "Estimated Reliability", xlim = range(all_h20_bi$emp[,1], all_h20_bi$kong[,1]), ylim = range(all_h20_bi$emp[,4], all_h20_bi$kong[,4]))
points(all_h20_bi$kong[,1], all_h20_bi$kong[,4], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)
dev.off()

## Univariate, sd = 0.2
pdf(file = "weibull-multirep_sd2_uni.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(all_h20_uni1$emp[,1], all_h20_uni1$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.2, MTTF vs Iterations", xlab = "Iterations", ylab = "Estimated MTTF", xlim = range(all_h20_uni1$emp[,1], all_h20_uni1$kong[,1]), ylim = range(all_h20_uni1$emp[,3], all_h20_uni1$kong[,3]))
points(all_h20_uni1$kong[,1], all_h20_uni1$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)

plot(all_h20_uni2$emp[,1], all_h20_uni2$emp[,3], pch=19, col = "blue", main = "Beta SD = 0.2, Reliability vs Iterations", xlab = "Iterations", ylab = "Estimated Reliability", xlim = range(all_h20_uni2$emp[,1], all_h20_uni2$kong[,1]), ylim = range(all_h20_uni2$emp[,3], all_h20_uni2$kong[,3]))
points(all_h20_uni2$kong[,1], all_h20_uni2$kong[,3], pch=19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1, pch=19)
dev.off()


########################################################
########################################################
########### ESS vs Beta standard deviation #############
########################################################
########################################################


N <- 1e4
SD <- seq(.1, .2, .01)
reps <- 10
ess.kong.mttf <- matrix(0, nrow = reps, ncol = length(SD))
ess.emp.mttf <- matrix(0, nrow = reps, ncol = length(SD))
ess.emp.bi <- matrix(0, nrow = reps, ncol = length(SD))
ess.kong.rel <- matrix(0, nrow = reps, ncol = length(SD))
ess.emp.rel <- matrix(0, nrow = reps, ncol = length(SD))

for (y in 1:length(SD))
{
  print(SD[y])
  for(i in 1:reps)
  {
    is <- joint_draw(N = N, t = proj.hour, h=SD[y])
    L <- is[[2]][,1]
    B <- is[[2]][,2]
    weights <- is[[1]]
    norm_weights <- weights/sum(weights)
    H1 <- mttf(L,B)
    H2 <- reliability(L,B)
    H3 <- as.matrix(cbind(H1, H2), ncol = 2)
    snis1 <- sum(H1 * norm_weights)
    snis2 <- sum(H2 * norm_weights)
    snis3 <- colSums(H3 * norm_weights)
    ## Kong's ESS
    ess.kong.mttf[i,y] <- 1/sum(norm_weights^2)
    ess.kong.rel[i,y] <- 1/sum(norm_weights^2)
    
    ## Empirical ESS
    varIbar1 <- sum((H1 - snis1)^2 * norm_weights)/N
    emp.var1 <- sum(norm_weights^2 * (H1 - snis1)^2)
    ess.emp.mttf[i,y] <- N*varIbar1/emp.var1
    varIbar2 <- sum((H2 - snis2)^2 * norm_weights)/N
    emp.var2 <- sum(norm_weights^2 * (H2 - snis2)^2)
    ess.emp.rel[i,y] <- N*varIbar2/emp.var2
    varIbar3 <- (t(H3 - snis3) %*% (norm_weights * (H3 - snis3)))/N
    emp.var3 <- (t(H3 - snis3) %*% (norm_weights^2 * (H3 - snis3)))
    ess.emp.bi[i,y] <- N*((det(varIbar3)/det(emp.var3))^(1/p))
  }
}

pdf(file = "weibull-ESSvsBetaSD.pdf", height = 5, width = 5)
plot(SD, colMeans(ess.emp.mttf/N), lwd=2, type = "l", main = "ESS/N vs Beta SD", col = "blue", ylab = "ESS/N", xlab = "Beta SD", ylim = range(colMeans(ess.emp.mttf)/N, colMeans(ess.emp.rel)/N, colMeans(ess.emp.bi/N), colMeans(ess.kong.mttf)/N))
lines(SD, colMeans(ess.emp.rel/N), col = "green", lwd=2)
lines(SD, colMeans(ess.emp.bi/N), col = "black", lwd=2)
lines(SD, colMeans(ess.kong.mttf/N), col = "orange", lwd=2)
legend("topright", legend = c("Estimated Bivariate", "Estimated MTTF", "Estimated Rel", "Kong"), lwd = 2, col = c("black","blue", "green", "orange"), lty=1)
dev.off()
