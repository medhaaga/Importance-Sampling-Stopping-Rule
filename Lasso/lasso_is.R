# loading required packages
library(INLA)
library(ISLR)
library(glmnet)
library(smoothmest)
library(mvtnorm)

# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("inlaIS.R")
source("genFuncs.R")

data(Hitters)

#Check NA's and fix

Hitters <- na.omit(Hitters)

#
# The Lasso
#

#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)
n.beta <- ncol(df$x)

# ml estimates
ml = summary(lm(y~-1 + x, data = df))$coefficients[,1:2]

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)

#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid,intercept = F)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1,intercept=F)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)


#Fitted values
lasso.fitted <- predict(out, s = bestlam, newx = x)
# importing dataset

# finding inverse of the precision
stdev.samp <- .25 * solve(t(x)%*%x)


fit.inla <- function(data, eta) {
  data$oset = data$x %*% eta
  res = inla(y ~ -1 + offset(oset), data = data)
  res = inla.rerun(res)
  return(list(mlik = res$mlik[[1]],
              dists = list(tau = res$marginals.hyperpar[[1]]),
              stats = list(tau = as.numeric(res$summary.hyperpar[1]))))
}



prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))
  
  if(!log) { res <- exp(res) }
  
  return(res)
}

# initial parameters of the proposal distribution
init = list(mu = rep(0,n.beta), cov = 4*stdev.samp)
# proposal distribution
## evaluate
dq.beta <- function(y, theta = init, log =TRUE) {
  #dmvnorm(y,mean = x, sigma = sigma,log = log)
  dmvt(y,delta=theta[[1]],sigma=theta[[2]],df=3,log=log,type = "shifted")
}
## sample
rq.beta <- function(theta) {
  #rmvnorm(1,mean=x,sigma = sigma)
  as.vector(rmvt(1,sigma = theta[[2]], df=3, delta = theta[[1]], type = "shifted"))
}


### IS
set.seed(1)
start.time <- Sys.time()
is_mod = inlaIS(data = df, init = init,
                prior.beta, dq.beta, rq.beta, fit.inla, N_0 = 2000, N = 10000,ncores = 5, N0iter = 10)
end.time <- Sys.time()
print(end.time - start.time)
save(is_mod, file = "Out/lasso_is_init.Rdata")
