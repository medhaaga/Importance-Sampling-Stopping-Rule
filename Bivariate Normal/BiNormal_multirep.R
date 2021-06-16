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
lambdas <- c(0.1, 0.8, 0.5)
rhos <- c(0.1, .7, .5)
Lambda <- diag(c(2,1))
Upsilon <- diag(c(2,2))
mu <- rep(2, p)

for(i in 1:3)
{
  Lambda[1,2] <- Lambda[2,1] <- lambdas[i]*sqrt(prod(diag(Lambda)))
  Upsilon[1,2] <- Upsilon[2,1] <- rhos[i]*sqrt(prod(diag(Upsilon)))
  true.uis.var <- var.uis(Upsilon, Lambda, mu)
  true.snis.var <- var.snis(Upsilon, Lambda, mu)
  print(det(Lambda)/det(true.uis.var))
  print(det(Lambda)/det(true.snis.var))
  pdf(file = paste("biNorm-ellipse_setting", i, ".pdf", sep = ""), height=5, width=5)
  plot(ellipse(Lambda), type = 'l', lwd=2, lty=2, xlab = "Component - 1", ylab = "Component - 2", xlim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.uis.var), ellipse(true.snis.var), -6, 6), ylim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.uis.var), ellipse(true.snis.var), -6, 6))
  lines(ellipse(Upsilon), col = "blue", lwd=2, lty=2)
  lines(ellipse(true.uis.var), col = "red", lwd=2, lty=1)
  lines(ellipse(true.snis.var), col = "orange", lwd=2, lty=1)
  legend("topleft", legend = c("Target variance", "Proposal variance", "UIS variance", "SNIS variance"), bty = 'n', cex=.8, col = c("black", "blue", "red", "orange"), lty=1, lwd=2)
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
  plot(rhos, true.ess.uis[1,], type = "l", lwd=2, ylim = range(true.ess.uis, true.ess.snis, 1.2), ylab = "True ESS/n", xlab = expression(rho))
  lines(rhos, true.ess.snis[1,], type = "l", lwd=2, col = "blue")
  for (i in 1:p)
  {
    lines(rhos, true.ess.uis[(i+1),], lty=2)
    lines(rhos, true.ess.snis[(i+1),], col = adjustcolor(col="blue", alpha.f = 1), lty=2)
  }
    
  
  abline(v = sel_rhos[t], lty=3, col = "red")
  legend("topright", legend = c("UIS", "SNIS"), bty = 'n', cex=.8, col = c("black", (adjustcolor("blue", alpha.f = 1))), lwd=2)
  dev.off()
}


##########################################################
##### Termination for diff epsilon in settings-1,2,3 #####
##########################################################

loop <- 100
dims <- c(2,10)

for (j in 1:2)
{
  for (a in 1:3)
  {
    p <- dims[j]
    mu <- rep(2, p)
    load(file = paste("ESS_eps_p", p, "_setting", a, ".Rdata", sep = ""))
    SE1_uis <- rowSums(abs(all_ESS1[[1]]$emp_uis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS1[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE2_uis <- rowSums(abs(all_ESS2[[1]]$emp_uis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS2[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE3_uis <- rowSums(abs(all_ESS3[[1]]$emp_uis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS3[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE4_uis <- rowSums(abs(all_ESS4[[1]]$emp_uis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS4[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    
    SE1_snis <- rowSums(abs(all_ESS1[[1]]$emp_snis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS1[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE2_snis <- rowSums(abs(all_ESS2[[1]]$emp_snis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS2[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE3_snis <- rowSums(abs(all_ESS3[[1]]$emp_snis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS3[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    SE4_snis <- rowSums(abs(all_ESS4[[1]]$emp_snis[,-c(1,2)]-mu))#apply(as.matrix(all_ESS4[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '1')
    
    pdf(file = paste("multinorm-l1VSiter_p", p, "_setting", a, ".pdf", sep = ""), width = 5, height = 5)
    
    plot(log(all_ESS1[[1]]$emp_uis[,1]), SE1_uis, pch=1, col = "blue", main = paste("p = ", p), xlim = log(range(truth.snis, truth.uis, all_ESS1[[1]]$emp_uis[,1],  all_ESS1[[1]]$emp_snis[,1], all_ESS4[[1]]$emp_uis[,1],  all_ESS4[[1]]$emp_snis[,1])), ylim = range(SE1_uis, SE2_uis, SE3_uis, SE4_uis, SE1_snis, SE2_snis, SE3_snis, SE4_snis), xlab = "Log of Iterations", ylab = "L1 norm of error")
    points(log(all_ESS2[[1]]$emp_uis[,1]), SE2_uis, pch=1, col = "green")
    points(log(all_ESS3[[1]]$emp_uis[,1]), SE3_uis, pch=1, col = "orange")
    points(log(all_ESS4[[1]]$emp_uis[,1]), SE4_uis, pch=1, col = "pink")
    
    points(log(all_ESS1[[1]]$emp_snis[,1]), SE1_snis, pch=19, col = "blue")
    points(log(all_ESS2[[1]]$emp_snis[,1]), SE2_snis, pch=19, col = "green")
    points(log(all_ESS3[[1]]$emp_snis[,1]), SE3_snis, pch=19, col = "orange")
    points(log(all_ESS4[[1]]$emp_snis[,1]), SE4_snis, pch=19, col = "pink")
    
    abline(v = log(truth.uis[1]), lty=2, col = "blue")
    abline(v = log(truth.uis[2]), lty=2, col = "green")
    abline(v = log(truth.uis[3]), lty=2, col = "orange")
    abline(v = log(truth.uis[4]), lty=2, col = "pink")
    
    abline(v = log(truth.snis[1]), lty=3, col = "blue")
    abline(v = log(truth.snis[2]), lty=3, col = "green")
    abline(v = log(truth.snis[3]), lty=3, col = "orange")
    abline(v = log(truth.snis[4]), lty=3, col = "pink")
    
    legend("topright", border = NA, bty = 'n', legend = c(expression(paste(epsilon, "=0.10")), expression(paste(epsilon, "=0.06")), expression(paste(epsilon, "=0.04")), expression(paste(epsilon, "=0.02"))), col = c("blue", "green", "orange", "pink"), pch = 19, cex = 1.2)
    legend("topright", inset = c(.25, 0), bty = 'n', legend = c("UIS", "SNIS"), pch = c(1,19))
    dev.off()
  }
  
}



