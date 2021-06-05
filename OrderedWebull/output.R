set.seed(1)
rm(list = ls())
source("functions.R")
library(corrplot)

########## Fish dataset1 #################
t <- c(83.5, 91.0, 91.0, 97.0, 107.0, 109.5, 114.0, 115.41, 128.61, 133.53, 138.58, 140.0, 152.08, 155.1)
t <- (t-80)/100
failures <- c(0,6,3,3,2)
stress <- c(0, .3, .5, .7, 1)
m <- 3
n <- 14
r <- 14
p <- m+2

####### hyperparams for priors ###########
val <- .8
shape.a <- val 
shape.b <- val 
ord.a <- val 
ord.b <- val 
ord.A <- rep(1, m+1)
p_H1 <- p
p_H2 <- p-1

n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])


########################################################
########### Correlation plot for IS estimator ########
########################################################


sim <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, 1e5)
H1 <- fun1(sim$samp)
H2 <- fun2(sim$samp)
colnames(H1) <- c("a", "l1", "l2", "l3", "l4")
colnames(H2) <- c("mu1", "mu2", "mu3", "mu4")
run_weights <- sim$weights
norm_weights <- run_weights/sum(run_weights)
snis1 <- colSums(H1*norm_weights)
snis2 <- colSums(H2*norm_weights)
emp.var1 <- (t(H1 - snis1) %*% (norm_weights^2 * (H1 - snis1)))
emp.var2 <- (t(H2 - snis2) %*% (norm_weights^2 * (H2 - snis2)))

pdf(file = "ISestimator_corrplot_H1.pdf", height=5, width =5)
corrplot(cov2cor(emp.var1), method = "circle")
dev.off()
pdf(file = "ISestimator_corrplot_H2.pdf", height=5, width =5)
corrplot(cov2cor(emp.var2), method = "circle")
dev.off()

###################################################################
########## ESS vs Sample size for kong, uni, and multi ESS ########
###################################################################


samp_size <- seq(1e3, 1e4, 1e3)
samps <- length(samp_size)
reps <- 5e2


load(file = "ESSvsSampSize_objects_kong_uni_mult.Rdata")
se.kong <- 2*apply(kong/samp_size, 2, sd)/sqrt(reps)
se.H1 <- 2*apply(multESS1/samp_size, 2, sd)/sqrt(reps)
se.H2 <- 2*apply(multESS2/samp_size, 2, sd)/sqrt(reps)
kong_means <- rowMeans(kong/samp_size)
multESS1_means <- rowMeans(multESS1/samp_size)
num.multESS1_means <- rowMeans(num.multESS1)
denom.multESS1_means <- rowMeans(denom.multESS1)
multESS2_means <- rowMeans(multESS2/samp_size)
num.multESS2_means <- rowMeans(num.multESS2)
denom.multESS2_means <- rowMeans(denom.multESS2)

pdf(file = "ESSvsSampSizeH1.pdf", height = 5, width = 5)
plot(samp_size, kong_means, type = 'l', lwd=2, col = "red", xlab = "Sample Size", ylab = "ESS/N", ylim = range(c(multESS2_means-se.H2, multESS1_means+se.H1, apply(apply(uniESS1, 3, rowMeans), 1, '/', samp_size))))
lines(samp_size, multESS1_means, col = "blue", lwd=2)
for(k in 1:5)
{
  lines(samp_size, apply(uniESS1, 3, rowMeans)[k,]/samp_size, col = "green")
}
legend("bottomleft", legend = c("Kong", "multiESS", "uniESS"), cex = 0.6, col = c("red", "blue", "green"), lty=1)
dev.off()


pdf(file = "ESSvsSampSizeH1H2.pdf", height = 5, width = 5)
par(xpd=TRUE)
plot(samp_size, kong_means, type = 'l',lwd=2, col = "red",  xlab = "Sample Size", ylab = "ESS/N", ylim = range(c(multESS2_means-se.H2, multESS1_means+se.H1, apply(apply(uniESS1, 3, rowMeans), 1, '/', samp_size))))
segments(x0 = samp_size, y0 = kong_means - se.kong, y1 = kong_means+se.kong, col = adjustcolor("red", alpha=.5))
lines(samp_size, multESS1_means, col = "blue", lwd=2)
segments(x0 = samp_size, y0 = multESS1_means-se.H1, y1 = multESS1_means+se.H1, col = adjustcolor("blue", alpha = 0.5))
lines(samp_size, num.multESS1_means, col = "lightblue", lwd=1)
lines(samp_size, denom.multESS1_means, col = "lightblue", lwd=1)
lines(samp_size, multESS2_means, col = "black", lwd=2)
segments(x0 = samp_size, y0 = multESS2_means-se.H2, y1 = multESS2_means+se.H2, col = adjustcolor("black", alpha = 0.5))
lines(samp_size, num.multESS2_means, col = "grey", lwd=1)
lines(samp_size, denom.multESS2_means, col = "grey", lwd=1)
legend("bottomleft", cex = 0.6, legend = c("Kong", "multiESS_H1", "multiESS_H2"), col = c("red", "blue", "black"), lty=1, lwd=2)
dev.off()


pdf(file = "ESSvsSampSizeH1H2_num_denom.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
plot(samp_size, num.multESS1_means, col = "blue", lwd=2, type='l', xlab = "Sample Size", ylab = "Num and Denom", ylim = range(num.multESS1_means, denom.multESS1_means))
lines(samp_size, denom.multESS1_means, col = "lightblue", lwd=2, type='l')
legend("topright", legend = c("Numerator", "Denominator"), col = c("blue", "lightblue"), lwd=2)

plot(samp_size, num.multESS2_means, col = "blue", lwd=2, type='l', xlab = "Sample Size", ylab = "Num and Denom", ylim = range(num.multESS2_means, denom.multESS2_means))
lines(samp_size, denom.multESS2_means, col = "lightblue", lwd=2, type='l')
legend("topright", legend = c("Numerator", "Denominator"), col = c("blue", "lightblue"), lwd=2)
dev.off()
###########################################################
############ Termination point vs epsilon ##############################
###########################################################

load(file = "ESS_eps_H1_H2.Rdata")

for (i in 1:p_H1)
  print(colMeans(all_ESS_H1[[i]]$emp))

for (i in 1:p_H2)
  print(colMeans(all_ESS_H2[[i]]$emp))

