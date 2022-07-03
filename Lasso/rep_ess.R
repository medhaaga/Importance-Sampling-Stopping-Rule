
load("Out/lasso_reps.Rdata")

d <- dim(all_out[[1]]$eta)[2]
N <- length(all_out[[1]]$weight)
reps <- 10


n_samp <- seq(1e3, 1e4, (N/10))
mESS <- matrix(0, nrow = length(n_samp), ncol = reps)
ESS <- array(0, dim = c(length(n_samp), d, reps))
oESS <- array(0, dim = c(length(n_samp), d, reps))
kESS <- matrix(0, nrow = length(n_samp), ncol = reps)


#############################################
###### Calculate ESS for each rep and n #####
#############################################

for (k in 1:reps){
  is_samp <- all_out[[k]]$eta
  is_weights <- all_out[[k]]$weight
  
  for (i in (1:length(n_samp))){
    n <- n_samp[i]
    H <- is_samp[1:n,]
    W <- is_weights[1:n]/sum(is_weights[1:n])
    snis <- colSums(H * W)
    ESS_num <- (t(H - snis) %*% (W * (H - snis)))/n
    ESS_denom <- (t(W*(H - snis)) %*% (W * (H - snis)))
    mESS[i,k] <- ((det(ESS_num)/det(ESS_denom))^(1/d))
    ESS[i,,k] <- (diag(ESS_num)/diag(ESS_denom))
    kESS[i,k] <- 1/(n*sum(W^2))

    Wh <- abs(H)*is_weights[1:n]/colSums(abs(H)*is_weights[1:n])
    oESS[i,,k] <- 1/(n*colSums(Wh^2))
  }
}

#############################################
### Calculate mean and SE of ESS ############
#############################################

mean_kESS <- rowMeans(kESS)
mean_mESS <- rowMeans(mESS)
mean_oESS <- apply(oESS, c(1,2), mean)
mean_ESS <- apply(ESS, c(1,2), mean)
se_kESS <- 2*apply(kESS, 1, sd)/sqrt(reps)
se_mESS <- 2*apply(mESS, 1, sd)/sqrt(reps)
#se_ESS1 <- 2*apply(ESS[,1,], 1, sd)/sqrt(reps)
#se_ESS2 <- 2*apply(ESS[,2,], 1, sd)/sqrt(reps)
#se_ESS3 <- 2*apply(ESS[,3,], 1, sd)/sqrt(reps)
#se_ESS4 <- 2*apply(ESS[,4,], 1, sd)/sqrt(reps)
#se_ESS5 <- 2*apply(ESS[,5,], 1, sd)/sqrt(reps)

### Plot

pdf(file = "Out/reps_ESS.pdf", height = 5, width  = 5)
plot(n_samp, mean_mESS, type = "l", lwd = 2, xlab = "Sample Size", ylab = "ESS/n", col = "blue", ylim = c(0,.31))
segments(x0 = n_samp, y0 = mean_mESS - se_mESS, y1 = mean_mESS + se_mESS, col = adjustcolor("blue", alpha = 0.5))
for(i in 1:d)
{
  lines(n_samp, mean_oESS[,i], lty=2, col = adjustcolor("black", alpha = 0.7))
}
lines(n_samp, mean_kESS, col = "red")
segments(x0 = n_samp, y0 = mean_kESS - se_kESS, y1 = mean_kESS + se_kESS, col = adjustcolor("red", alpha = 0.5))
legend("topright", legend = c("m-ESS", "Owen", "Kong"), lty = c(1,2,1), lwd = c(2,1,1), col = c("blue", "black", "red"), cex = 0.8)
dev.off()

