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

is_ESS <- function(min_ESS, step, loop, N_min, m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, fun, h=0, p){
  
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
    is.samp <- fun(is.samp)
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
      more.is.samp <- fun(more.samp$samp)
      is.more <- rbind(is.more, more.is.samp)
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

fun1 <- function(x) x

fun2 <- function(is.samp)
{
  p <- dim(is.samp)[2]
  return(is.samp[,2:p]^(-1/is.samp[,1]))*gamma(1 + (1/is.samp[,1]))
}
