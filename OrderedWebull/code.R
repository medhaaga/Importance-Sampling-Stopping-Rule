set.seed(1)
library(mcmcse)
library(corrplot)
################ posterior density (unnormalised) ###############
posterior <- function(alpha, Lambda, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, n_bar, stress, t)
{
  lambda <- sum(Lambda)
  D <- rep(0, m+1)
  post1 <- 1
  for (i in 1:m+1){
    D[i] <- ((n - n_bar[i+1])*(stress[i+1]^alpha)) + ((n - n_bar[i])*(stress[i]^alpha)) +  sum(t[(n_bar[i]+1):n_bar[i+1]]^alpha)
    post1 <- post1 * (Lambda[i]^failures[i]) * exp(-Lambda[i]*(ord.b + D[i]))
  }
  
  perms <- permutations(n=m+1, r=m+1, v=1:m+1, repeats.allowed = FALSE)
  post2 <- 0
  for(k in 1:factorial(m+1)){
    post2 <- post2 + Lambda[perms[k,]]^(ord.A - 1)
  }
  
  post <- (lambda^(ord.a - sum(ord.A))) * post1 * post2 * (alpha^(r+shape.a-1)) * (exp(-shape.b*alpha)) * prod(t^alpha)
}


############## \pi_1(lambda | alpha, data) (refer paper) #############
pi1_samp <- function(a, b, A, m){
  
  X <- rgamma(n=1, shape=a, rate=b)
  Z <- rep(0, (m+1))
  for (i in 1:(m+1))
    Z[i] <- rgamma(n=1, shape=A[i], rate=1)
  Y <- Z/sum(Z)
  Lambda <- Y*X
  Lambda <- sort(Lambda, decreasing = FALSE)
  return(Lambda)
}


############## g function (refer paper) ########################
g <- function(lambda, failures, J, D, S, ord.a, ord.b, m){
  foo <- 1
  for (i in 1:m+1){
    foo <- foo * ((lambda[i]^(failures[i+1] - J)) * exp(-lambda[i] * (D[i] - S)))
  }
  foo <- foo/((ord.b + S)^(ord.a + (m+1)*J))
}


############## function returns a list of #######################
############## IS samples and corresponsding weights ############

is <- function(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, M){
  J <- min(failures[-1])
  D <- rep(0, m+1)
  samp <- matrix(0, nrow = M, ncol = m+2)
  weights <- rep(0, M)
  
  for (k in 1:M){
    alpha <- rgamma(n=1, shape = r+shape.a, rate = shape.b - sum(log(t)))
    for (i in 1:(m+1)){
      D[i] <- ((n - n_bar[i+1])*(stress[i+1]^alpha)) - ((n - n_bar[i])*(stress[i]^alpha)) + sum(t[(n_bar[i]+1):n_bar[i+1]]^alpha)
    }
    S <- min(D)
    lambda <- pi1_samp((ord.a + (m+1)*J), (ord.b + S), (ord.A + J), m)
    samp[k,] <- c(alpha, lambda)
    weights[k] <- g(lambda, failures, J, D, S, ord.a, ord.b, m)
  }
  weights <- weights #/sum(weights)
  return (list("samp" = samp, "weights" = weights))
}

#########################################################
######## function calculates emp ESS and kong ESS #######
######## and terminates acc to minESS criterion #########
#########################################################

is_ESS <- function(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, h=0, p){
  
  iter.emp <- rep(N_min, loop)
  iter.kong <- rep(N_min, loop)
  
  if (h==0){
    means.emp <- matrix(0, nrow = loop, ncol = p)
    means.kong <- matrix(0, nrow = loop, ncol = p)
  } else{
    means.emp <- matrix(0, nrow = loop, ncol = 1)
    means.kong <- matrix(0, nrow = loop, ncol = 1)
  }
  ess.emp <- rep(0, loop)
  ess.kong <- rep(0, loop)
  
  for(j in 1:loop){
    
    print(j)
    is.samp_weights <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, N_min)
    is.samp <- is.samp_weights$samp
    is.more <- is.samp
    run_weights <- is.samp_weights$weights
    norm_weights <- run_weights / sum(run_weights )
    
    if (h==0)
      H <- as.matrix(is.samp, ncol = p) 
    for (i in 1:p){
      if(h==i)
        H <- as.matrix(is.samp[,i], ncol = 1)
    }
    if(h==0)
      eff_p = p
    else
      eff_p = 1
    
    snis <- colSums(H * norm_weights)
    
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[j] <- N_min*((det(varIbar)/det(emp.var))^(1/eff_p))
    ess.kong[j]  <- 1/sum(norm_weights^2)  ##comment to not calculate kong
    means.emp[j,] <- snis
    means.kong[j,] <- snis   ##comment to not calculate kong
    
    while (ess.emp[j] <= min_ess ||ess.kong[j] <= min_ess){
      #print(paste("Emp Iterations: ", iter.emp[j], "Emp ESS: ", ess.emp[j], "Kong Iterations: ", iter.kong[j], "Kong ESS: ",  ess.kong[j]))
      more.samp <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, step)
      is.more <- rbind(is.more, more.samp$samp)
      if (h==0)
        H <- as.matrix(is.more, ncol = p) 
      for (i in 1:p){
        if(h==i)
          H <- as.matrix(is.more[,i], ncol = 1)
      }
      run_weights <- c(run_weights, more.samp$weights)
      norm_weights <- run_weights / sum(run_weights )
      snis <- colSums(H * norm_weights)
      
      if (ess.emp[j] <= min_ess){
        iter.emp[j] <- iter.emp[j] + step
        varIbar <- (t(H - snis)%*%(norm_weights*(H-snis)))/iter.emp[j]
        emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
        ess.emp[j] <- iter.emp[j]*((det(varIbar)/det(emp.var))^(1/eff_p))
        means.emp[j,] <- snis
      }
      
      if (ess.kong[j] <= min_ess){
        iter.kong[j] <- iter.kong[j] + step
        ess.kong[j] <- 1/sum(norm_weights^2)
        means.kong[j,] <- snis
      }
      
      
    }
  }
  return(list("kong" = cbind(iter.kong, ess.kong, means.kong), 
              "emp" = cbind(iter.emp, ess.emp, means.emp)))
}

########## Fish dataset1 #################
t <- c(83.5, 91.0, 91.0, 97.0, 107.0, 109.5, 114.0, 115.41, 128.61, 133.53, 138.58, 140.0, 152.08, 155.1)
t <- (t-80)/100
failures <- c(0,6,3,3,2)
stress <- c(0, .3, .5, .7, 1)
m <- 3
n <- 14
r <- 14
p <- m+2

####### hyperparams for non-informative priors ###########
val <- .5
shape.a <- val 
shape.b <- val 
ord.a <- val 
ord.b <- val 
ord.A <- rep(1, m+1)


n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])


########################################################
########### Correlation plot for IS estimator ########
########################################################


sim <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, 1e5)
H1 <- sim$samp
H2 <- (H1[,2:5]^(-1/H1[,1]))*gamma(1 + (1/H1[,1]))
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

reps <- 5e2
samp_size <- seq(1e3, 1e5, 1e4)
samps <- length(samp_size)
kong <- matrix(0, nrow = samps, ncol = reps)
multESS1 <- matrix(0, nrow = samps, ncol = reps)
uniESS1 <- array(0, dim = c(p, reps, samps))
multESS2 <- matrix(0, nrow = samps, ncol = reps)
uniESS2 <- array(0, dim = c(4, reps, samps))

for (l in 1:reps)
{
  print(l)
  sim <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, samp_size[length(samp_size)])
  is.samp <- sim$samp
  run_weights <- sim$weights
  
  for (z in 1:samps)
  {
    M <- samp_size[z]
    H1 <- is.samp[1:M,]
    H2 <- (is.samp[1:M,2:5]^(-1/is.samp[1:M,1]))*gamma(1 + (1/is.samp[1:M,1]))
    norm_weights <- run_weights[1:M]/sum(run_weights[1:M])
    
    kong[z,l] <- 1/sum(norm_weights^2)
    
    snis1 <- colSums(H1 * norm_weights)
    snis2 <- colSums(H2 * norm_weights)
    
    varIbar1 <- (t(H1 - snis1) %*% (norm_weights * (H1 - snis1)))/M
    varIbar2 <- (t(H2 - snis2) %*% (norm_weights * (H2 - snis2)))/M
    
    emp.var1 <- (t(H1 - snis1) %*% (norm_weights^2 * (H1 - snis1)))
    emp.var2 <- (t(H2 - snis2) %*% (norm_weights^2 * (H2 - snis2)))
    
    p1 <- length(snis1)
    p2 <- length(snis2)
    
    multESS1[z,l] <- M*((det(varIbar1)/det(emp.var1))^(1/p1))
    uniESS1[,l,z] <- M* diag(varIbar1)/ diag(emp.var1)
    
    multESS2[z,l] <- M*((det(varIbar2)/det(emp.var2))^(1/p2))
    uniESS2[,l,z] <- M* diag(varIbar2)/ diag(emp.var2)
    
  }
}
save(kong, uniESS1, uniESS2, multESS1, multESS2, file = "ESSvsSampSize_objects_kong_uni_mult.Rdata")



############### Run from here after loading objects #################
load(file = "ESSvsSampSize_objects_kong_uni_mult.Rdata")
se.kong <- 2*apply(kong/samp_size, 2, sd)
se.H1 <- 2*apply(multESS1/samp_size, 2, sd)/sqrt(reps)
se.H2 <- 2*apply(multESS2/samp_size, 2, sd)/sqrt(reps)
kong_means <- rowMeans(kong/samp_size)
multESS1_means <- rowMeans(multESS1/samp_size)
multESS2_means <- rowMeans(multESS2/samp_size)

pdf(file = "ESSvsSampSizeH1.pdf", height = 5, width = 5)
plot(samp_size, kong_means, type = 'l', lwd=2, col = "red", xlab = "Sample Size", ylab = "ESS/N", ylim = range(c(kong/samp_size, multESS1/samp_size, apply(apply(uniESS1, 3, rowMeans), 1, '/', samp_size))))
segments(x0 = samp_size, y0 = kong_means - se.kong, y1 = kong_means+se.kong, col = adjustcolor("red", alpha=.5))
lines(samp_size, multESS1_means, col = "blue", lwd=2)
segments(x0 = samp_size, y0 = multESS1_means-se.H1, y1 = multESS1_means+se.H1, col = adjustcolor("blue", alpha = 0.5))
for(k in 1:5)
{
  lines(samp_size, apply(uniESS1, 3, rowMeans)[k,]/samp_size, col = "green")
}
legend("topleft", legend = c("Kong", "multiESS", "uniESS"), cex = 0.8, col = c("red", "blue", "green"), lty=1)
dev.off()


pdf(file = "ESSvsSampSizeH1H2.pdf", height = 5, width = 5)
plot(samp_size, kong_means, type = 'l',lwd=2, col = "red",  xlab = "Sample Size", ylab = "ESS/N", ylim = range(c(multESS2_means-se.H2, multESS1_means+se.H1)))
segments(x0 = samp_size, y0 = kong_means - se.kong, y1 = kong_means+se.kong, col = adjustcolor("red", alpha=.5))
lines(samp_size, multESS1_means, col = "blue", lwd=2)
segments(x0 = samp_size, y0 = multESS1_means-se.H1, y1 = multESS1_means+se.H1, col = adjustcolor("blue", alpha = 0.5))
lines(samp_size, multESS2_means, col = "black", lwd=2)
segments(x0 = samp_size, y0 = multESS2_means-se.H2, y1 = multESS2_means+se.H2, col = adjustcolor("black", alpha = 0.5))
legend("bottomleft", cex = 0.8, legend = c("Kong", "multiESS_H1", "multiESS_H2"), col = c("red", "blue", "black"), lty=1, lwd=2)
dev.off()

###########################################################
############ Termination point vs epsilon ##############################
###########################################################

step <- 1e2
loop <- 10


### Epsilon = 0.1
min_ess <- minESS(p, eps = 0.1)
N_min <- round(min_ess/2)
all_ESS1 <- list()
all_ESS1[[1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, h=0, p=p)
for (i in 1:5){
  min_ess <- minESS(1, eps = 0.1)
  N_min <- round(min_ess/2)
  all_ESS1[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, h=i, p)
}
save(all_ESS1, file = paste("ESS_eps1_p", p, ".Rdata", sep = ""))

print(all_ESS1[[1]]$emp)
print(all_ESS1[[2]]$emp)
print(all_ESS1[[3]]$emp)
print(all_ESS1[[4]]$emp)
print(all_ESS1[[5]]$emp)
print(all_ESS1[[6]]$emp)
