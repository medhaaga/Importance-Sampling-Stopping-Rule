objective <- function(par){
  sig <- par[1]
  lambda <- par[2]
  num <- (sig^4) * (1 - lambda^2)
  denom <- 4*(sig^4) - 4*(sig^4)*(lambda^2) - 4*(sig^2) + 1
  return(num^3/denom^2)
}


corr <- seq(-1,1,.01)
optim_corr <- corr
for(i in 1:length(corr)) optim_corr[i] <- optim(c(sqrt(2), corr[i]), objective, method = "BFGS")$par[2]
plot(corr, optim_corr, "l")

obj_corr <- corr
for(i in 1:length(corr)) obj_corr[i] <- objective(c(sqrt(1.5), corr[i]))
pdf(file = "detVScorr.pdf")
plot(corr, obj_corr, "l", ylim = c(0,10), ylab = "Determinant of variance matrix", xlab = "Correlation")
abline(h = obj_corr[101], lty=2, col = "red")
dev.off()
vari <- seq(sqrt(3)^-1, 10, .1)
obj_var <- vari
for(i in 1:length(vari)) obj_var[i] <- objective(c(vari[i], 0))
plot(vari, obj_var, "l")
