library(mcmcse)


h <- function(x,m) x^m
m <- 3
shape = 5
prop_shape = 4
var.snis <- function (a1, a2, m){
  mu <- gamma(a1 + m)/gamma(a1)
  var.snis <- (gamma(a2)/(gamma(a1))^2) * (gamma(2*a1 - a2 + 2*m) - (2*mu*gamma(2*a1 - a2 + m)) + (gamma(2*a1 - a2)*mu^2)) 
  return (var.snis)
}

min_ess <- minESS(1)
step <- 10
N_min <- 5e3
loop <- 100
iter.emp <- rep(N_min, loop)
iter.truth <- rep(N_min, loop)
iter.kong <- rep(N_min, loop)
means.emp <- rep(0, loop)
means.truth <- rep(0, loop)
means.kong <- rep(0, loop)
ess.emp <- rep(0, loop)
ess.truth <- rep(0, loop)
ess.kong <- rep(0, loop)
for(t in 1:loop){
  
  print(t)
  is <- rgamma(n=N_min, shape=prop_shape, rate=1)
  is.more <- is
  weights <- dgamma(is, shape=shape, rate=1)/dgamma(is, shape=prop_shape, rate=1)
  norm_weights <- weights/sum(weights)
  snis <- sum(h(is,m)*weights)/sum(weights)
  varIbar <- sum((h(is,m) - snis)^2 * norm_weights)/N_min
  emp.var <- sum(norm_weights^2 * (h(is,m) - snis)^2)
  ess.emp[t] <- N_min*varIbar/emp.var
  true.var.is <- var.snis(shape, prop_shape, m)
  var.f <- (gamma(shape + 2*m)/gamma(shape)) - (gamma(shape + m)/gamma(shape))^2
  ess.truth[t] <- N_min*(var.f)*(1/true.var.is)
  ess.kong[t]  = 1/sum(norm_weights^2)
  means.emp[t] <- snis
  means.kong[t] <- snis
  means.truth[t] <- snis
  
  while (ess.emp[t] <= min_ess || ess.truth[t] <= min_ess || ess.kong[t] <= min_ess){
    
    is.more <- c(is.more, rgamma(n=step, shape = prop_shape, rate=1))
    weights <- dgamma(is.more, shape=shape, rate=1)/dgamma(is.more, shape=prop_shape, rate=1)
    norm_weights <- weights/sum(weights)
    snis <- sum(h(is.more,m)*weights)/sum(weights)
    
    if (ess.emp[t] <= min_ess){
      iter.emp[t] <- iter.emp[t] + step
      varIbar <- sum((h(is.more,m) - snis)^2 * norm_weights)/iter.emp[t]
      emp.var <- sum(norm_weights^2 * (h(is.more,m) - snis)^2)
      ess.emp[t] <- iter.emp[t]*varIbar/emp.var
      means.emp[t] <- snis
    }
    
    if (ess.truth[t] <= min_ess){
      iter.truth[t] <- iter.truth[t] + step
      ess.truth[t] <- iter.truth[t]*(var.f)*(1/true.var.is)
      means.truth[t] <- snis
    }
    
    if (ess.kong[t] <= min_ess){
      iter.kong[t] <- iter.kong[t] + step
      ess.kong[t] <- 1/sum(norm_weights^2)
      means.kong[t] <- snis
    }
    
  }
}
  
save(ess.emp, ess.truth, ess.kong, iter.emp, iter.truth, iter.kong, means.emp, means.truth, means.kong, file = "gamma_multi_rep_m3_alpha4.Rdata")



## plots for alpha = 4 and alpha = 7 combined
pdf(file = "gamma-multirep_m3.pdf", height = 5, width = 10)
load(file = "gamma_multi_rep_m3_alpha4.Rdata")

iter.emp1 <-  iter.emp
iter.kong1 <- iter.kong
iter.truth1 <- iter.truth[1]
means.emp1 <- means.emp
means.kong1 <- means.kong

load(file = "gamma_multi_rep_m3_alpha7.Rdata")
iter.emp2 <-  iter.emp
iter.kong2 <- iter.kong
iter.truth2 <- iter.truth[1]
means.emp2 <- means.emp
means.kong2 <- means.kong

plot(iter.emp1, means.emp1, pch=19, col = "blue", xlab = "Iterations", ylab = "Target Estimate", xlim = range(iter.emp1, iter.kong1, iter.emp2, iter.kong2), ylim = range(means.emp1, means.kong1, means.emp2, means.kong2))
points(iter.kong1, means.kong1, pch = 19, col = "orange")
abline(v = iter.truth1, lty=2)

points(iter.emp2, means.emp2, pch = 19, col = "lightblue")
points(iter.kong2, means.kong2, pch = 19, col = "pink")
abline(v = iter.truth2, lty=2)
legend("topright", legend = c("Estimated, a = 4", "Kong, a = 4", "Estimated, a = 7", "Kong, a = 7"), cex=1.5, col = c("blue", "orange", "lightblue", "pink"), pch  =19)
dev.off()
