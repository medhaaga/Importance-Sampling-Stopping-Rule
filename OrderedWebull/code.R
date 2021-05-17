set.seed(1)
library(mcmcse)

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
  weights <- weights/sum(weights)
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
    norm_weights <- is.samp_weights$weights
    
    if (h==0)
      H <- as.matrix(is.samp, ncol = p) 
    for (i in 1:p){
      if(h==i)
        H <- as.matrix(is.samp[,i], ncol = 1)
    }
    
    snis <- colSums(H * norm_weights)
    
    varIbar <- (t(H - snis) %*% (norm_weights * (H - snis)))/N_min
    emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
    
    ess.emp[j] <- N_min*((det(varIbar)/det(emp.var))^(1/p))
    ess.kong[j]  <- 1/sum(norm_weights^2)  ##comment to not calculate kong
    means.emp[j,] <- snis
    means.kong[j,] <- snis   ##comment to not calculate kong
    
    while (ess.emp[j] <= min_ess ||ess.kong[j] <= min_ess){
      print(paste("Emp Iterations: ", iter.emp[j], "Emp ESS: ", ess.emp[j], "Kong Iterations: ", iter.kong[j], "Kong ESS: ",  ess.kong[j]))
      more.samp <- is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, step)
      is.more <- rbind(is.more, more.samp$samp)
      if (h==0)
        H <- as.matrix(is.more, ncol = p) 
      for (i in 1:p){
        if(h==i)
          H <- as.matrix(is.more[,i], ncol = 1)
      }
      norm_weights <- c(norm_weights, more.samp$weights)
      snis <- colSums(H * norm_weights)
      
      if (ess.emp[j] <= min_ess){
        iter.emp[j] <- iter.emp[j] + step
        varIbar <- (t(H - snis)%*%(norm_weights*(H-snis)))/iter.emp[j]
        emp.var <- (t(H - snis) %*% (norm_weights^2 * (H - snis)))
        ess.emp[j] <- iter.emp[j]*((det(varIbar)/det(emp.var))^(1/p))
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
shape.a <- 1
shape.b <- 1
ord.a <- 1
ord.b <- 1
ord.A <- rep(1, m+1)


n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])


step <- 1e3
loop <- 20


###########################################################
############ Termination point vs epsilon ##############################
###########################################################

### Epsilon = 0.1
min_ess <- minESS(p, eps = 0.1)
N_min <- round(min_ess/2)
all_ESS1 <- list()
all_ESS1[[1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, h=0, p=p)
for (i in 1:p){
  min_ess <- minESS(1, eps = 0.1)
  all_ESS1[[i+1]] <- is_ESS(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, h=i, p=1)
}
save(all_ESS1, file = paste("ESS_eps1_p", p, ".Rdata", sep = ""))


