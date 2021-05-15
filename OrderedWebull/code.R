
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
g <- function(lambda, failures, J, D, S, org.a, org.b, m){
  foo <- 1
  for (i in 1:m+1){
    foo <- foo * ((lambda[i]^(failures[i+1] - J)) * exp(-lambda[i] * (D[i] - S)))
  }
  foo <- foo/((org.b + S)^(org.a + (m+1)*J))
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
      D[i] <- ((n - n_bar[i+1])*(stress[i+1]^alpha)) + ((n - n_bar[i])*(stress[i]^alpha)) + sum(t[(n_bar[i]+1):n_bar[i+1]]^alpha)
    }
    S <- min(D)
    lambda <- pi1_samp((org.a + (m+1)*J), (org.b + S), (org.A + J), m)
    samp[k,] <- c(alpha, lambda)
    weights[k] <- g(lambda, failures, J, D, S, org.a, org.b, m)
  }
  weights <- weights/sum(weights)
  return (list("samp" = samp, "weights" = weights))
}


########## Fish dataset1 #################
t <- c(83.5, 91.0, 91.0, 97.0, 107.0, 109.5, 114.0, 115.41, 128.61, 133.53, 138.58, 140.0, 152.08, 155.1)
t <- (t-80)/100
failures <- c(0,6,3,3,2)
stress <- c(0, .3, .5, .7, 1)
m <- 3
n <- 14
r <- 14


####### hyperparams for non-informative priors ###########
shape.a <- 1e-3
shape.b <- 1e-3
org.a <- 1e-3
org.b <- 1e-3
org.A <- rep(1, m+1)


n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])

M <- 1000
is.samp_weights = is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, M)
colSums(is.samp_weights$samp*is.samp_weights$weights)




t <- c(91, 93, 94, 98.2, 115.81, 116, 116.5, 117.25, 126.75, 127.5, 154.33, 159.5, 164, 184.14, 188.33)
t <- (t-80)/150
failures <- c(0, 4, 6, 0, 3, 2)
stress <- c(0, .2, .33, .46, .6, 1)
m <- 4
n <- 15
r <- 15


####### hyperparams for non-informative priors ###########
shape.a <- 1e-3
shape.b <- 1e-3
org.a <- 1e-3
org.b <- 1e-3
org.A <- rep(1, m+1)


n_bar <- rep(0, m+2)
for (i in 1:(m+2)) 
  n_bar[i] <- sum(failures[1:i])

M <- 1000
is.samp_weights = is(m, n, r, shape.a, shape.b, ord.a, ord.b, ord.A, failures, stress, t, n_bar, M)
colSums(is.samp_weights$samp*is.samp_weights$weights)

