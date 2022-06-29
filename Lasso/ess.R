library(ggplot2)
library(scales)
library(ggpubr)
library(ISLR)

load("Out/lasso_is.Rdata")



data(Hitters)
Hitters <- na.omit(Hitters)

is_samp <- is_mod$eta
is_weights <- is_mod$weight
N <- length(is_weights)

n_samp <- seq(1e3, N, (N/10))
mESS <- rep(0, length(n_samp))
ESS <- matrix(0, nrow = 5, ncol = length(n_samp))
kESS <- rep(0, length(n_samp))

for (i in (1:length(n_samp))){
  n <- n_samp[i]
  H <- is_samp[1:n,]
  W <- is_weights[1:n]/sum(is_weights[1:n])
  snis <- colSums(H * W)
  ESS_num <- (t(H - snis) %*% (W * (H - snis)))/n
  ESS_denom <- (t(W*(H - snis)) %*% (W * (H - snis)))
  mESS[i] <- n*(det(ESS_num)/det(ESS_denom))^(1/5)
  ESS[,i] <- n*(diag(ESS_num)/diag(ESS_denom))
  kESS[i] <- 1/(sum(W^2))
}

pdf(file = "Out/is_mESS.pdf", height = 5, width  = 5)
plot(n_samp, mESS, type = "l", lwd = 2, xlab = "n", main = "IS-INLA ESS vs n", ylab = "Effective Sample Size", col = "blue", ylim = range(mESS, ESS, kESS))
lines(n_samp, ESS[1,], lty=2)
lines(n_samp, ESS[2,], lty=2)
lines(n_samp, ESS[3,], lty=2)
lines(n_samp, ESS[4,], lty=2)
lines(n_samp, ESS[5,], lty=2)
lines(n_samp, kESS, col = "red")
legend("topleft", legend = c("m-ESS", "ESS", "k-ESS"), lty = c(1,2,1), lwd = c(2,1,1), col = c("blue", "black", "red"), cex = 0.8)
dev.off()
