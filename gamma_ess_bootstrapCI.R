library(mcmcse)

shape <- 5  #target alpha
prop_shape <- seq(4, 9, .5)  #proposal alphas
N <- 5e3    #Number of samples
ess.truth <- rep(0, length(prop_shape))
ess.kong <- rep(0, length(prop_shape))
ess.emp <- rep(0, length(prop_shape))
ess.emp.ci <- matrix(0, nrow = 2, ncol = length(prop_shape))
ess.kong.ci <- matrix(0, nrow = 2, ncol = length(prop_shape))
m <- 3  # moment of Gamma distr. to be estimated using IS
h <- function(x,m) x^m  # x^m

######### function calculates true N*Variance of SNIS estiamtor ###########
var.snis <- function (a1, a2, m){
  mu <- gamma(a1 + m)/gamma(a1)
  var.snis <- (gamma(a2)/(gamma(a1))^2) * (gamma(2*a1 - a2 + 2*m) - (2*mu*gamma(2*a1 - a2 + m)) + (gamma(2*a1 - a2)*mu^2)) 
  return (var.snis)
}


for (y in 1:length(prop_shape)){
  
  is <- rgamma(n=N, shape=prop_shape[y], rate=1)
  weights <- dgamma(is, shape=shape, rate=1)/dgamma(is, shape=prop_shape[y], rate=1)
  norm_weights <- weights/sum(weights)
  snis <- sum(h(is,m)*weights)/sum(weights)
  
  ## Kong's ESS
  ess.kong[y] <- 1/sum(norm_weights^2)
  
  ## True ESS
  var.f <- (gamma(shape + 2*m)/gamma(shape)) - (gamma(shape + m)/gamma(shape))^2
  true.var.is <- var.snis(shape, prop_shape[y], m)
  ess.truth[y] <- N*(var.f)*(1/true.var.is)
  
  ## Empirical ESS
  varIbar <- sum((h(is,m) - snis)^2 * norm_weights)/N
  emp.var <- sum(norm_weights^2 * (h(is,m) - snis)^2)
  ess.emp[y] <- N*varIbar/emp.var
  
  print(paste(ess.truth[y], ess.emp[y], ess.kong[y]))
  
  B <- N/2
  boot.emp.ess <- rep(0, B)
  boot.kong.ess <- rep(0, B)
  for (j in 1:B){
    sample_idx <- sample(1:N, size = N, replace = TRUE)
    is.boot <- is[sample_idx]
    norm.weights.boot <- norm_weights[sample_idx]
    snis.boot <- sum(h(is.boot,m)*norm.weights.boot)
    varIbar <- sum((h(is.boot,m) - snis.boot)^2 * norm.weights.boot)/N
    emp.var <- sum(norm.weights.boot^2 * (h(is.boot,m) - snis.boot)^2)
    boot.kong.ess[j] <- 1/(sum(norm.weights.boot^2))
    boot.emp.ess[j] <- N*varIbar/emp.var
  }
  ess.emp.ci[,y] <- quantile(boot.emp.ess, probs = c(.05, .95)) #var(boot.emp.ess)
  ess.kong.ci[,y] <- quantile(boot.kong.ess, probs = c(.05, .95)) #var(boot.kong.ess)
  
}

pdf(file = "gamma_m3_proposals_CI.pdf")
plot(prop_shape, ess.emp/N, type = "l", col = "blue", ylim = range(ess.emp, ess.truth, ess.kong)/N)
lines(prop_shape, ess.kong/N, col = "orange")
lines(prop_shape, ess.truth/N, col = "green")
segments(x0 = prop_shape, y0 = ess.emp.ci[1,]/N, x1 = prop_shape, y1 = ess.emp.ci[2,]/N, col = adjustcolor("blue", alpha.f=.4))
segments(x0 = prop_shape, y0 = ess.kong.ci[1,]/N, x1 = prop_shape, y1 = ess.kong.ci[2,]/N, col = adjustcolor("orange", alpha.f=.4))
legend("topright", legend = c("Empirical", "Kong", "Truth"), col = c("blue", "orange", "green"), lwd=2)
dev.off()

