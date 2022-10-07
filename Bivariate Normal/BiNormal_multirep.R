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
lambdas <- c(0.1, 0.5, .8)
rhos <- c(0.1, .5, .7)
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

  pdf(file = paste("Out/biNorm-ellipse_setting", i, ".pdf", sep = ""), height=6, width=6)
  plot(ellipse(Lambda), type = 'l', lwd=2, lty=2, xlab = "Component 1", 
    ylab = "Component 2", cex.lab = 1.3,
    xlim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.uis.var), ellipse(true.snis.var), -6, 6), 
    ylim = range(ellipse(Lambda), ellipse(Upsilon), ellipse(true.uis.var), ellipse(true.snis.var), -4.5, 6))
  lines(ellipse(Upsilon), col = "blue", lwd=2, lty=2)
  lines(ellipse(true.uis.var), col = "red", lwd=2, lty=1)
  lines(ellipse(true.snis.var), col = "orange", lwd=2, lty=1)
  legend("topleft", legend = c("Target variance", "Proposal variance", "UIS variance", "SNIS variance"), 
    bty = 'n', cex=1.3, col = c("black", "blue", "red", "orange"), 
    lty=c(2,2,1,1), lwd=2)
  dev.off()
}


############################################################
######### True ESS/n vs rho for diff settings @#############
################# in bivariate Normal case##################
############################################################

p <- 2
sel_rhos <- c(.1, .5, .7)
rho_max <- c(.6, .7, .8)

for (t in 1:3)
{
  load(file = paste("Out/trueESSvsRho_p2_setting", t, ".Rdata", sep=""))
  rhos <- seq(0, rho_max[t], .1)

  pdf(file = paste("Out/biNormal-trueESSvsRho_p2_setting", t, ".pdf", sep = ""), width = 6, height = 6)
  plot(rhos, true.ess.uis[1,], type = "l", col = "red", lwd=2, 
    ylim = range(true.ess.uis, true.ess.snis, 1.2), ylab = "True ESS/n", 
    xlab = expression(rho), cex.lab = 1.3)
  lines(rhos, true.ess.snis[1,], type = "l", lwd=2, col = "orange")
  for (i in 1:p)
  {
    lines(rhos, true.ess.uis[(i+1),], lty=2, col = "red")
    lines(rhos, true.ess.snis[(i+1),], col = adjustcolor(col="orange", alpha.f = 1), lty=2)
  }


  abline(v = sel_rhos[t], lty=3)
  if(t==1)
    legend("topright", legend = c("UIS", "SNIS"), bty = 'n', cex=1.3, col = c("red", "orange"), lwd=2) else
      legend("topleft", legend = c("UIS", "SNIS"), bty = 'n', cex=1.3, col = c("red", "orange"), lwd=2)

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
    load(file = paste("Out/ESS_eps_p", p, "_setting", a, ".Rdata", sep = ""))
    SE1_uis <- apply(as.matrix(all_ESS1[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '2')
    SE2_uis <- apply(as.matrix(all_ESS2[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '2')
    SE3_uis <- apply(as.matrix(all_ESS3[[1]]$emp_uis[,-c(1,2)]-mu, row = loop), 1, norm, '2')

    SE1_snis <- apply(as.matrix(all_ESS1[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '2')
    SE2_snis <- apply(as.matrix(all_ESS2[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '2')
    SE3_snis <- apply(as.matrix(all_ESS3[[1]]$emp_snis[,-c(1,2)]-mu, row = loop), 1, norm, '2')

    if(p==2)
      ymax <- 0.08 else
        ymax <- 0.3

    pdf(file = paste("Out/multinorm-sqeVSiter_p", p, "_setting", a, ".pdf", sep = ""), width = 6, height = 6)

    plot(log(all_ESS1[[1]]$emp_uis[,1]), SE1_uis, pch=1, col = "blue", main = paste("p = ", p), 
      xlim = log(range(truth.snis, truth.uis, all_ESS1[[1]]$emp_uis[,1],  all_ESS1[[1]]$emp_snis[,1], all_ESS3[[1]]$emp_uis[,1],  all_ESS3[[1]]$emp_snis[,1])), 
      ylim = range(ymax, SE1_uis, SE2_uis, SE3_uis, SE1_snis, SE2_snis, SE3_snis), xlab = "Number of Samples", 
      ylab = "Squared error", xaxt = "n", cex.lab = 1.3)

    x_range <- log(range(truth.snis, truth.uis, all_ESS1[[1]]$emp_uis[,1],  all_ESS1[[1]]$emp_snis[,1], all_ESS3[[1]]$emp_uis[,1],  all_ESS3[[1]]$emp_snis[,1]))
    axis(1, at=seq(x_range[1], x_range[2], length.out=6), labels=round(exp(seq(x_range[1], x_range[2], length.out=6))))
    points(log(all_ESS2[[1]]$emp_uis[,1]), SE2_uis, pch=1, col = "green")
    points(log(all_ESS3[[1]]$emp_uis[,1]), SE3_uis, pch=1, col = "orange")

    points(log(all_ESS1[[1]]$emp_snis[,1]), SE1_snis, pch=2, col = "blue")
    points(log(all_ESS2[[1]]$emp_snis[,1]), SE2_snis, pch=2, col = "green")
    points(log(all_ESS3[[1]]$emp_snis[,1]), SE3_snis, pch=2, col = "orange")

    abline(v = log(truth.uis[1]), lty=2, col = "blue")
    abline(v = log(truth.uis[2]), lty=2, col = "green")
    abline(v = log(truth.uis[3]), lty=2, col = "orange")

    abline(v = log(truth.snis[1]), lty=3, col = "blue")
    abline(v = log(truth.snis[2]), lty=3, col = "green")
    abline(v = log(truth.snis[3]), lty=3, col = "orange")

    legend("topright", border = NA, bty = 'n', 
      legend = c(expression(paste(epsilon, "=0.06")), 
        expression(paste(epsilon, "=0.04")),
         expression(paste(epsilon, "=0.02"))), 
      col = c("blue", "green", "orange"), pch = 19, cex = 1.3)
    legend("topright", inset = c(.25, 0.02), bty = 'n', 
      legend = c("UIS", "SNIS"), pch = c(1,2), cex = 1.3)
    dev.off()
  }

}



