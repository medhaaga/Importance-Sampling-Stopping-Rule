set.seed(1)
shape <- 5  #target alpha
prop_shape <- seq(4, 9, .5)  #proposal alphas
N <- c(1e3, 1e4)    #Number of samples
ess <- matrix(0, nrow = length(N), ncol = length(prop_shape))
ess.snis <- matrix(0, nrow = length(N), ncol = length(prop_shape))
ess.boot <- matrix(0, nrow = length(N), ncol = length(prop_shape))
ess.emp <- matrix(0, nrow = length(N), ncol = length(prop_shape))

m <- 3  # moment of Gamma distr. to be estimated using IS
h <- function(x,m) x^m  # x^m

######### function calculates true N*Variance of IS estiamtor ###########
var.is <- function(a1, a2, m){
  var.is <- (((gamma(a2)/(gamma(a1))^2) * gamma(2*a1 - a2 + 2*m)) - (gamma(a1 + m)/gamma(a1))^2) 
  return (var.is)
}

######### function calculates true N*Variance of SNIS estiamtor ###########
var.snis <- function (a1, a2, m){
  mu <- gamma(a1 + m)/gamma(a1)
  var.snis <- (gamma(a2)/(gamma(a1))^2) * (gamma(2*a1 - a2 + 2*m) - (2*mu*gamma(2*a1 - a2 + m)) + (gamma(2*a1 - a2)*mu^2)) 
  return (var.snis)
}

######### function calculates true N*MSE of SNIS estiamtor ###########

mse.snis <- function (a1, a2, m, N){
  mu <- gamma(a1 + m)/gamma(a1)
  var.snis <- (gamma(a2)/(gamma(a1))^2) * (gamma(2*a1 - a2 + 2*m) - (2*mu*gamma(2*a1 - a2 + m)) + (gamma(2*a1 - a2)*mu^2)) 
  bias <- (1/N) * (gamma(a2)/(gamma(a1))^2) * ((mu * gamma(2*a1 - a2)) - gamma(2*a1 - a2 + m))
  mse.snis <- var.snis + bias^2
  return (mse.snis)
}


##################################################
########## Weighted Importance Sampling# #########
##################################################

for (x in 1:length(N)){
  for (y in 1:length(prop_shape)){
    
    true.var.is <- var.snis(shape, prop_shape[y], m)
    is <- rgamma(n=N[x], shape=prop_shape[y], rate=1)
    weights <- dgamma(is, shape=shape, rate=1)/dgamma(is, shape=prop_shape[y], rate=1)
    norm_weights <- weights/sum(weights)
    snis <- sum(h(is,m)*weights)/sum(weights)
    ess.snis[x,y] <- 1/sum(norm_weights^2)
    var.f <- (gamma(shape + 2*m)/gamma(shape)) - (gamma(shape + m)/gamma(shape))^2
    ess[x,y] <- N[x]*(var.f)*(1/true.var.is)
    B <- 500
    boot_means <- rep(0, B)
    
    varIbar <- sum((h(is,m) - snis)^2 * weights)/(sum(weights)*N[x])
    emp.var <- sum(norm_weights^2 * (h(is,m) - snis)^2)
    ess.emp[x,y] <- N[x]*(varIbar/emp.var)
    
    for (i in 1:B){
      sample_idx <- sample(1:N[x], size= N[x], replace = TRUE)
      boot_means[i] <- sum(h(is,m)[sample_idx] * weights[sample_idx])/sum(weights[sample_idx])
    }
    ess.boot[x,y] <- N[x] * varIbar/( var(boot_means))
    print(paste(ess[x,y], ess.snis[x,y], ess.boot[x,y]))
  }
  ess.boot[x,ess.boot[x,]==Inf] <- 5e29
  ess[x,ess[x,]==Inf] <- 2*max(ess.boot[x,])
  
}

# plots for different sample size. x axis is the shape parameter of proposal.

## N=1e3
par(mfrow = c(1,2))
plot(prop_shape, ess.snis[1,]/N[1], type = "l", main = "No of samples = 1e2"
     , ylim = range(ess.snis[1,]/N[1], ess.boot[1,]/N[1], ess[1,]/N[1]), xlab = "alpha_2", ylab = "ESS/N")
lines(prop_shape, ess[1,]/N[1], lty=2, col = "red")
lines(prop_shape, ess.boot[1,]/N[1], lty=3, col = "blue")
lines(prop_shape, ess.emp[1,]/N[1], lty=4, col = "green3")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth", "Empirical"), lty = c(1,3,2,4), col = c("black", "blue", "red", "green3"))
abline(v=shape, col = "pink")
abline(v=prop_shape[ess.boot[1,]==max(ess.boot[1,])], col = "pink")


## N=1e4
plot(prop_shape, ess.snis[2,]/N[2], type = "l", main = "No of samples = 1e3"
     , ylim = range(ess.snis[2,]/N[2], ess.boot[2,]/N[2], ess[2,]/N[2]), xlab = "alpha_2", ylab = "ESS/N")
lines(prop_shape, ess[2,]/N[2], lty=2, col = "red")
lines(prop_shape, ess.boot[2,]/N[2], lty=3, col = "blue")
lines(prop_shape, ess.emp[2,]/N[2], lty=4, col = "green3")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth", "Empirical"), lty = c(1,3,2,4), col = c("black", "blue", "red", "green3"))
abline(v = shape, col = "pink")
abline(v=prop_shape[ess.boot[2,]==max(ess.boot[2,])], col = "pink")

