set.seed(1)
rm(list= ls())
source("functions.R")
library(mvtnorm)
library(ellipse)
library(mcmcse)


############################################################
########### Covariance ellipse for diff settings ###########
########### Traget, proposal, and SNIS variance ############
############################################################

p <- 2
lambdas <- c(0.1, 0.9, 0.5)
rhos <- c(0.1, .7, .5)
Lambda <- diag(c(3,1))
Upsilon <- diag(c(3,3))


for(i in 1:3)
{
  Lambda[1,2] <- Lambda[2,1] <- lambdas[i]*sqrt(prod(diag(Lambda)))
  Upsilon[1,2] <- Upsilon[2,1] <- rhos[i]*sqrt(prod(diag(Upsilon)))
  true.var <- var.is(Upsilon, Lambda)
  pdf(file = paste("biNorm-ellipse_setting", i, ".pdf", sep = ""))
  plot(ellipse(Lambda), type = 'l', xlab = "Component - 1", ylab = "Component - 2", xlim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.var)), ylim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.var)))
  lines(ellipse(Upsilon), col = "blue")
  lines(ellipse(true.var), col = "red")
  legend("topleft", legend = c("Target", "proposal", "SNIS variance"), col = c("black", "blue", "red"), lty=1)
  dev.off()
}


############################################################
######### True ESS/n vs rho for diff settings @#############
################# in bivariate Normal case##################
############################################################

p <- 2
sel_rhos <- c(.1, .7, .5)
rho_max <- c(.6, .8, .7)

for (t in 1:3)
{
  load(file = paste("trueESSvsRho_p2_setting", t, ".Rdata", sep=""))
  rhos <- seq(0, rho_max[t], .1)
  pdf(file = paste("biNormal-trueESSvsRho_p2_setting", t, ".pdf", sep = ""), width = 5, height = 5)
  par(xpd = TRUE)
  plot(rhos, true.ess[1,], type = "l", lwd=2, ylim = range(true.ess, 1.2), ylab = "True ESS/n", xlab = expression(rho))
  for (i in 1:p)
    lines(rhos, true.ess[(i+1),], col = adjustcolor(col="blue", alpha.f = 1), lwd=2)
  abline(v = sel_rhos[t], lty=2)
  legend("bottomleft", inset= c(.5,1), legend = c("Multivariate", "Univariate"), col = c("black", (adjustcolor("blue", alpha.f = 1))), lwd=2, cex = 1.2)
  dev.off()
}



#####################################################
###### Univarate and multivate true ESS vs rho  ######
########### for dimensions p = 2, 3, 10 ##############
#####################################################

dims <- c(2,3,10)
rhos <- seq(.3, .7, .01)
for (i in 1:length(dims))
{
  p <- dims[i]
  load(file = paste("trueESSvsRho_p", p, ".Rdata", sep=""))
  
  pdf(file = paste("trueESSvsRho_p", p, ".pdf", sep=""))
  plot(rhos, true.ess[1,], type = 'l', lwd=2, ylim = range(true.ess, 1.2), ylab = "True ESS/n", xlab = expression(rho))
  for (t in 1:p)
    lines(rhos, true.ess[t+1,], col = "blue")
  legend("topleft", legend = c("Multivariate", "Univariate"), col = c("black", "blue"), lwd=c(2,1), cex = 1.5)
  dev.off()
}

##########################################################
##### Termination for diff epsilon in settings-1,2,3 #####
##########################################################

loop <- 20
dims <- c(2,3,10)

for (j in 1:3)
{
  p <- dims[j]
  load(file = paste("ESS_eps_p", p, ".Rdata", sep = ""))
  SE1 <- matrix(0, nrow = p+1, ncol = loop)
  for (i in 1:(p+1))
    SE1[i,] <- apply(as.matrix(all_ESS1[[i]]$emp[,3], row = loop), 1, norm, '2')
  SE2 <- matrix(0, nrow = p+1, ncol = loop)
  for (i in 1:(p+1))
    SE2[i,] <- apply(as.matrix(all_ESS2[[i]]$emp[,3], row = loop), 1, norm, '2')
  SE3 <- matrix(0, nrow = p+1, ncol = loop)
  for (i in 1:(p+1))
    SE3[i,] <- apply(as.matrix(all_ESS3[[i]]$emp[,3], row = loop), 1, norm, '2')
  SE4 <- matrix(0, nrow = p+1, ncol = loop)
  for (i in 1:(p+1))
    SE4[i,] <- apply(as.matrix(all_ESS4[[i]]$emp[,3], row = loop), 1, norm, '2')
  SE5 <- matrix(0, nrow = p+1, ncol = loop)
  for (i in 1:(p+1))
    SE5[i,] <- apply(as.matrix(all_ESS5[[i]]$emp[,3], row = loop), 1, norm, '2')
  
  pdf(file = paste("multivariateNormal-sqeVSiter_p", p, ".pdf", sep = ""), width = 5, height = 5)
  
  plot(all_ESS1[[1]]$emp[,1], SE1[1,], pch=19, col = "blue", main = paste("p = ", p), xlim = range(truth[1], truth[5], all_ESS1[[1]]$emp[,1], all_ESS5[[1]]$emp[,1], all_ESS1[[2]]$emp[,1], all_ESS5[[2]]$emp[,1]), ylim = range(SE1, SE2, SE3, SE4, 0.12), xlab = "Iterations", ylab = "Squared Error")
  points(all_ESS2[[1]]$emp[,1], SE2[1,], pch=19, col = "green")
  points(all_ESS3[[1]]$emp[,1], SE3[1,], pch=19, col = "orange")
  points(all_ESS4[[1]]$emp[,1], SE4[1,], pch=19, col = "pink")
  points(all_ESS5[[1]]$emp[,1], SE5[1,], pch=19, col = "red")
   
  abline(v = truth[1], lty=2, col = "blue")
  abline(v = truth[2], lty=2, col = "green")
  abline(v = truth[3], lty=2, col = "orange")
  abline(v = truth[4], lty=2, col = "pink")
  abline(v = truth[5], lty=2, col = "red")
  legend("topright", inset = c(.1,0), legend = c(expression(paste(epsilon, "=0.1")), expression(paste(epsilon, "=0.08")), expression(paste(epsilon, "=0.06")), expression(paste(epsilon, "=0.04")), expression(paste(epsilon, "=0.03"))), col = c("blue", "green", "orange", "pink", "red"), pch = 19, cex = 1.2)
  dev.off()
}




load(file = "uisVSsnis.Rdata")
samp_size <- seq(1e3, 1e4, 1e3)
reps <- 1e2
samps <- length(samp_size)
kong.avg <- colMeans(ESS.kong)
uis.avg <- colMeans(ESS.uis)
snis.avg <- colMeans(ESS.snis)
uis.avg.num <- colMeans(ESS.uis.num)
snis.avg.num <- colMeans(ESS.snis.num)
uis.avg.denom <- colMeans(ESS.uis.denom)
snis.avg.denom <- colMeans(ESS.snis.denom)
uis.se <- apply(ESS.uis, MARGIN=2, FUN=sd)/sqrt(reps)
snis.se <- apply(ESS.snis, MARGIN=2, FUN=sd)/sqrt(reps)
kong.se <- apply(ESS.kong, MARGIN=2, FUN=sd)/sqrt(reps)

pdf(file = "uisVSsnis.pdf", width = 5, height = 5)
plot(samp_size, uis.avg, type = 'l', lwd=2, col = "blue", ylim = range(kong.avg-kong.se, uis.avg, snis.avg+snis.se))
lines(samp_size, snis.avg, col = "green", lwd=2)
lines(samp_size, kong.avg, lwd=2)
segments(x0 = samp_size, y0=uis.avg-uis.se, y1 = uis.avg+uis.se, col = "blue")
segments(x0 = samp_size, y0=snis.avg-snis.se, y1 = snis.avg+snis.se, col = "green")
segments(x0 = samp_size, y0=kong.avg-kong.se, y1 = kong.avg+kong.se, col = "black")
abline(h = trueESS.uis, lty=2, col = "blue")
abline(h = trueESS.snis, lty=2, col = "green")
dev.off()