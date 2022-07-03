library(ggplot2)
library(scales)
library(ggpubr)


load("Out/lasso_is_init.Rdata")



data(Hitters)
Hitters <- na.omit(Hitters)


is_samp <- is_mod$eta
is_weights <- is_mod$weight
N <- nrow(is_samp)
n_samp <- seq(1e3, N, (N/10))
mESS <- rep(0, length(n_samp))
ESS <- matrix(0, nrow = 19, ncol = length(n_samp))
kESS <- rep(0, length(n_samp))

for (i in (1:length(n_samp))){
  n <- n_samp[i]
  H <- is_samp[1:n,]
  W <- is_weights[1:n]/sum(is_weights[1:n])
  snis <- colSums(H * W)
  ESS_num <- (t(H - snis) %*% (W * (H - snis)))/n
  ESS_denom <- (t(W*(H - snis)) %*% (W * (H - snis)))
  mESS[i] <- n*(det(ESS_num)/det(ESS_denom))^(1/19)
  ESS[,i] <- n*(diag(ESS_num)/diag(ESS_denom))
  kESS[i] <- 1/(sum(W^2))
}

pdf(file = "Out/is_mESS.pdf", height = 5, width  = 5)
plot(n_samp, mESS/N, type = "l", lwd = 2, xlab = "n", main = "IS-INLA ESS vs n", ylab = "Effective Sample Size", col = "blue", ylim = range(mESS, ESS, kESS)/N)
for(i in 1:19)
{
	lines(n_samp, ESS[i,]/N)
}

lines(n_samp, kESS/N, col = "red")
legend("topright", legend = c("m-ESS", "ESS", "k-ESS"), lwd = c(2,1,1), col = c("blue", "black", "red"), cex = 0.8)
dev.off()