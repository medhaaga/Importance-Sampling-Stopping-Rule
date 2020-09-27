set.seed(1)
prop_sigma <- seq((1/sqrt(2))+.1, 2.5, .1)  #proposal sigmas
N <- c(1e2, 1e3, 1e4)    #Number of samples
ess <- matrix(0, nrow = 3, ncol = length(prop_sigma))
ess.snis <- matrix(0, nrow = 3, ncol = length(prop_sigma))
ess.boot <- matrix(0, nrow = 3, ncol = length(prop_sigma))

m <- 1  # moment of Normal distr. to be estimated using IS
h <- function(x,m) x^m  # x^m

######### function calculates true N*Variance of IS estiamtor ###########
var.is <- function(sigma){
  v <- sigma^2
  var.is <- (v^2)/((2*v - 1)^1.5)
  return (var.is)
}

######### function calculates true N*Variance of SNIS estiamtor ###########
var.snis <- function(sigma){
  v <- sigma^2
  var.is <- (v^2)/((2*v - 1)^1.5)
  return (var.is)
}

##############################################
###### Simple Importance Sampling ############
##############################################

for (x in 1:length(N)){
  for (y in 1:length(prop_sigma)){
    
    true.var.is <- var.is(prop_sigma[y])  ## N*true IS variance
    var.f <- 1  # true Var_p(h(x))
    is <- rnorm(n=N[x], mean=0, sd = (prop_sigma[y]))   ## Importance samples
    weights <- dnorm(is, mean=0, sd=1)/dnorm(is, mean=0, sd = (prop_sigma[y]))   
    norm_weights <- weights/sum(weights)
    snis <- mean(h(is,m)*weights)
    
    ess.snis[x,y] <- 1/sum(norm_weights^2)      ## estimated SNIS ESS 
    ess[x,y] <- N[x]*(var.f)*(1/true.var.is)
    
    B <- 500
    boot_means <- rep(0, B)
    
    varIbar <- mean((h(is,m) - snis)^2 * weights)/N[x]   ## empirical Var_p(h(x))
    for (i in 1:B){
      sample_idx <- sample(1:N[x], size= N[x], replace = TRUE)
      boot_means[i] <- mean(h(is,m)[sample_idx] * weights[sample_idx])
    }
    ess.boot[x,y] <- N[x] * varIbar/( var(boot_means))
    print(paste(ess[x,y], ess.snis[x,y], ess.boot[x,y]))
  }
}
save(ess.snis, ess, ess.boot, file = paste("Out_gaussian/gaussian-IS_m", m, ".Rdata", sep = ""))

# plots for different sample size

## N=1e2
pdf(file = paste("Out_gaussian/gaussian-IS_m", m, ".pdf", sep = ""), height = 4, width = 10)

par(mfrow = c(1,3))
plot(prop_sigma, ess.snis[1,]/N[1], type = "l", main = "No of samples = 1e2"
     , ylim = range(ess.snis[1,]/N[1], ess.boot[1,]/N[1]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[1,]/N[1], lty=2, col = "red")
lines(prop_sigma, ess.boot[1,]/N[1], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=sqrt(2), col = "pink")

## N=1e3
plot(prop_sigma, ess.snis[2,]/N[2], type = "l", main = "No of samples = 1e3"
     , ylim = range(ess.snis[2,]/N[2], ess.boot[2,]/N[2]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[2,]/N[2], lty=2, col = "red")
lines(prop_sigma, ess.boot[2,]/N[2], lty=3, col = "blue")
abline(v = sqrt(2), col = "pink")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))

## N=1e4
plot(prop_sigma, ess.snis[3,]/N[3], type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.snis[3,]/N[3], ess.boot[3,]/N[3]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[3,]/N[3], lty=2, col = "red")
lines(prop_sigma, ess.boot[3,]/N[3], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = sqrt(2), col = "pink")
dev.off()

##################################################
########## Weighted Importance Sampling# #########
##################################################

for (x in 1:length(N)){
  for (y in 1:length(prop_sigma)){
    
    true.var.is <- var.snis(prop_sigma[y])
    is <- rnorm(n=N[x], mean=0, sd=prop_sigma[y])
    weights <- dnorm(is, mean=0, sd=1)/dnorm(is, mean=0, sd = (prop_sigma[y]))   
    norm_weights <- weights/sum(weights)
    snis <- sum(h(is,m)*weights)/sum(weights)
    ess.snis[x,y] <- 1/sum(norm_weights^2)
    var.f <- 1
    ess[x,y] <- N[x]*(var.f)*(1/true.var.is)
    B <- 500
    boot_means <- rep(0, B)
    
    varIbar <- sum((h(is,m) - snis)^2 * weights)/(sum(weights)*N[x])
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
save(ess.snis, ess, ess.boot, file = paste("Out_gaussian/gaussian-SNIS_m", m, ".Rdata", sep = ""))


# plots for different sample size

## N=1e2
pdf(file = paste("Out_gaussian/gaussian-SNIS_m", m, ".pdf", sep = ""), height = 4, width = 10)
par(mfrow = c(1,3))
plot(prop_sigma, ess.snis[1,]/N[1], type = "l", main = "No of samples = 1e2"
     , ylim = range(ess.snis[1,]/N[1], ess.boot[1,]/N[1], ess[1,]/N[1]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[1,]/N[1], lty=2, col = "red")
lines(prop_sigma, ess.boot[1,]/N[1], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=sqrt(2), col = "pink")


## N=1e3
plot(prop_sigma, ess.snis[2,]/N[2], type = "l", main = "No of samples = 1e3"
     , ylim = range(ess.snis[2,]/N[2], ess.boot[2,]/N[2], ess[2,]/N[2]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[2,]/N[2], lty=2, col = "red")
lines(prop_sigma, ess.boot[2,]/N[2], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v=sqrt(2), col = "pink")


## N=1e4
plot(prop_sigma, ess.snis[3,]/N[3], type = "l", main = "No of samples = 1e4"
     , ylim = range(ess.snis[3,]/N[3], ess.boot[3,]/N[3], ess[3,]/N[3]), xlab = "sigma", ylab = "ESS/N")
lines(prop_sigma, ess[3,]/N[3], lty=2, col = "red")
lines(prop_sigma, ess.boot[3,]/N[3], lty=3, col = "blue")
legend("topright", legend = c("SNIS", "Bootstrap", "Truth"), lty = c(1,3,2), col = c("black", "blue", "red"))
abline(v = sqrt(2), col = "pink")
dev.off()

 
