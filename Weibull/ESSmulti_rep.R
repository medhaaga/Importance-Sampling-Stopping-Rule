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

h <- function(lambda, beta)
{
  lambda^(-1/beta) *gamma(1 + 1/beta)
}

is_ESS <- function(min_ESS, step, loop, N_min, proj.hour, beta_sd){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  means.emp <- rep(0, loop)
  means.kong <- rep(0, loop)
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
    H <- h(L,B)
    snis <- sum(H * norm_weights)
    
    varIbar <- sum((H - snis)^2 * norm_weights)/N_min
    emp.var <- sum(norm_weights^2 * (H - snis)^2)
    
    ess.emp[t] <- N_min*varIbar/emp.var
    ess.kong[t]  = 1/sum(norm_weights^2)
    means.emp[t] <- snis
    means.kong[t] <- snis
    
    while (ess.emp[t] <= min_ess ||ess.kong[t] <= min_ess){
      
      samp <- joint_draw(N = step, t = proj.hour, h=beta_sd)
      is.more[[1]] <- c(is.more[[1]], samp[[1]])
      is.more[[2]] <- rbind(is.more[[2]], samp[[2]])
      L <- is.more[[2]][,1]
      B <- is.more[[2]][,2]
      H <- h(L, B)
      weights <- is.more[[1]]
      norm_weights <- weights/sum(weights)
      snis <- sum(H * norm_weights)
      
      if (ess.emp[t] <= min_ess){
        iter.emp[t] <- iter.emp[t] + step
        varIbar <- sum((H - snis)^2 * norm_weights)/iter.emp[t]
        emp.var <- sum(norm_weights^2 * (H - snis)^2)
        ess.emp[t] <- iter.emp[t]*varIbar/emp.var
        means.emp[t] <- snis
      }
      
      if (ess.kong[t] <= min_ess){
        iter.kong[t] <- iter.kong[t] + step
        ess.kong[t] <- 1/sum(norm_weights^2)
        means.kong[t] <- snis
      }
      
    }
  }
  return(list("kong" = cbind(iter.kong, means.kong, ess.kong), "emp" = cbind(iter.emp, means.emp, ess.emp)))
}

min_ess <- minESS(1)
step <- 10
N_min <- 5e3
loop <- 100


# beta_sd = .2 and .1
all_h20 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .2)
all_h10 <- is_ESS(min_ESS, step, loop, N_min, proj.hour, beta_sd = .1)
save(all_h10, all_h20, file = "weibull-multirep.Rdata")


pdf(file = "weibull-multirep.pdf", height = 5, width = 10)
par (mfrow = c(1,2))
plot(all_h10$emp[,1], all_h10$emp[,2], pch=19, col = "blue", main = "Beta SD = 0.1", xlab = "Iterations", ylab = "Target Estimate", xlim = range(all_h10$emp[,1], all_h10$kong[,1]), ylim = range(all_h20$emp[,2], all_h20$kong[,2], all_h10$emp[,2], all_h10$kong[,2]))
points(all_h10$kong[,1], all_h10$kong[,2], pch = 19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1.2, pch=19)

plot(all_h20$emp[,1], all_h20$emp[,2], pch=19, col = "blue", main = "Beta SD = 0.2", xlab = "Iterations", ylab = "Target Estimate", xlim = range(all_h20$emp[,1], all_h20$kong[,1]), ylim = range(all_h20$emp[,2], all_h20$kong[,2], all_h10$emp[,2], all_h10$kong[,2]))
points(all_h20$kong[,1], all_h20$kong[,2], pch = 19, col = "orange")
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4"), col = c("blue", "orange"), cex=1.2, pch=19)
dev.off()

########################################################
########################################################
########### ESS vs Bets standard deviation #############
########################################################
########################################################


N <- 1e4
SD <- seq(.1, .25, .01)
reps <- 10
ess.kong <- matrix(0, nrow = reps, ncol = length(SD))
ess.emp <- matrix(0, nrow = reps, ncol = length(SD))

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
    H <- h(L,B)
    snis <- sum(H * norm_weights)
    
    ## Kong's ESS
    ess.kong[i,y] <- 1/sum(norm_weights^2)
    
    ## Empirical ESS
    varIbar <- sum((H - snis)^2 * norm_weights)/N
    emp.var <- sum(norm_weights^2 * (H - snis)^2)
    ess.emp[i,y] <- N*varIbar/emp.var
  }
}

pdf(file = "weibull-ESSvsBetaSD.pdf", height = 5, width = 5)
plot(SD, colMeans(ess.emp/N), type = "l", col = "blue", ylab = "ESS/N", xlab = "h parameter", ylim = range(colMeans(ess.emp), colMeans(ess.kong))/N)
lines(SD, colMeans(ess.kong/N), col = "orange")
legend("topright", legend = c("Estimated", "Kong"), col = c("blue", "orange"), lty=1)
dev.off

