x <- seq(-10,10,.1)
tar1 <- .5*dnorm(x, mean = 5, sd=1) + .5*dnorm(x, mean = -5, sd=1)
prop <- dcauchy(x, scale = 2)
plot(x, tar1, "l", ylim = range(tar1, tar2))
lines(x, prop, col = "red")
###############################################################

h <- function(x) x
B = 100
N = 1e4
sep <- seq(0, 6, .1)
ess.emp <- matrix(0, nrow = B, ncol = length(sep))
ess.kong <- matrix(0, nrow = B, ncol = length(sep))

for (b in 1:B){
  print(b)
  for (i in 1:length(sep))
  {
    print(i)
    is <- rcauchy(n = N, location = 0, scale = 2)
    weights <- (.5*dnorm(is, mean = (-1*sep[i]), sd=1) + .5*dnorm(is, mean = sep[i], sd=1))/dcauchy(is, location = 0, scale = 2)
    norm.weights <- weights/sum(weights)
    snis <- sum(h(is)*norm.weights)
    ess.kong[b,i] <- 1/(sum(norm.weights^2))
    
    varIbar <- sum((h(is) - snis)^2 * norm.weights)/N
    emp.var <- sum(norm.weights^2 * (h(is) - snis)^2)
    ess.emp[b,i] <- N*varIbar/emp.var
    
  }
}

emp.var <- apply(ess.emp, 2, var)
kong.var <- apply(ess.kong, 2, var)


pdf(file = "sepVSess_noBootstrap.pdf", height = 5, width = 5)
plot(sep, colMeans(ess.emp), "l", ylim = range(ess.kong+ess.kong.var, ess.emp+ess.emp.var, ess.kong-ess.kong.var, ess.emp-ess.emp.var))
lines(sep, colMeans(ess.kong), col = "red")
segments(x0 = sep, y0 = colMeans(ess.emp)-.5*emp.var, x1 = sep, y1 = colMeans(ess.emp)+.5*emp.var, col = adjustcolor("black", alpha.f=.4))
segments(x0 = sep, y0 = colMeans(ess.kong)-.5*kong.var, x1 = sep, y1 = colMeans(ess.kong)+.5*kong.var, col = adjustcolor("red", alpha.f=.4))
legend("topright", legend = c("Empirical", "Kong"), col = c("black", "red"), lwd=2)
dev.off()
