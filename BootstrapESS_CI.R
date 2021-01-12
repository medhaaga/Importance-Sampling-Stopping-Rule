x <- seq(-10,10,.1)
tar1 <- .5*dnorm(x, mean = 5, sd=1) + .5*dnorm(x, mean = -5, sd=1)
prop <- dcauchy(x, scale = 2)
plot(x, tar1, "l", ylim = range(tar1))
lines(x, prop, col = "red")
###############################################################

h <- function(x) x

N = 1e4
sep <- seq(0, 6, .1)
ess.emp <- sep
ess.emp.ci <- matrix(0, nrow = 2, ncol = length(sep))
ess.kong <- sep
ess.kong.ci <- matrix(0, nrow = 2, ncol = length(sep))

for (i in 1:length(sep))
{
  print(i)
  is <- rcauchy(n = N, location = 0, scale = 2)
  weights <- (.5*dnorm(is, mean = (-1*sep[i]), sd=1) + .5*dnorm(is, mean = sep[i], sd=1))/dcauchy(is, location = 0, scale = 2)
  norm.weights <- weights/sum(weights)
  snis <- sum(h(is)*norm.weights)
  ess.kong[i] <- 1/(sum(norm.weights^2))
  
  varIbar <- sum((h(is) - snis)^2 * norm.weights)/N
  emp.var <- sum(norm.weights^2 * (h(is) - snis)^2)
  ess.emp[i] <- N*varIbar/emp.var
  
  
  B <- N/2
  boot.emp.ess <- rep(0, B)
  boot.kong.ess <- rep(0, B)
  for (j in 1:B){
    sample_idx <- sample(1:N, size= N, replace = TRUE)
    is.boot <- is[sample_idx]
    norm.weights.boot <- norm.weights[sample_idx]
    snis.boot <- sum(h(is.boot)*norm.weights.boot)
    varIbar <- sum((h(is.boot) - snis.boot)^2 * norm.weights.boot)/N
    emp.var <- sum(norm.weights.boot^2 * (h(is.boot) - snis.boot)^2)
    boot.kong.ess[j] <- 1/(sum(norm.weights.boot^2))
    boot.emp.ess[j] <- N*varIbar/emp.var
  }
  ess.emp.ci[,i] <- quantile(boot.emp.ess, probs = c(.05, .95)) #var(boot.emp.ess)
  ess.kong.ci[,i] <- quantile(boot.kong.ess, probs = c(.05, .95)) #var(boot.kong.ess)
}

pdf(file = "sepVSess.pdf", height = 5, width = 5)
plot(sep, ess.emp/N, "l", ylim = range(ess.emp.ci,ess.kong.ci)/N)
lines(sep, ess.kong/N, col = "red")
segments(x0 = sep, y0 = ess.emp.ci[1,]/N, x1 = sep, y1 = ess.emp.ci[2,]/N, col = adjustcolor("black", alpha.f=.4))
segments(x0 = sep, y0 = ess.kong.ci[1,]/N, x1 = sep, y1 = ess.kong.ci[2,]/N, col = adjustcolor("red", alpha.f=.4))
legend("topright", legend = c("Empirical", "Kong"), col = c("black", "red"), lwd=2)
dev.off()
