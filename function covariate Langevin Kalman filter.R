# Load packages
# library(Rhabit)
library(here)
source(here("functions/Rhabit_functions.R"))
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(optimParallel)

##############################
# Making Covariate Functions #
##############################

# Covariate maps parameters -------------------------------------------------

IP1 <- 6;           IP2 <- 6
DP1 <- c(40, 40);   DP2 <- c(40, 40)
CP1 <- c(0, 0);     CP2 <- c(-2, pi/2)
FP1 <- c(0.6, 0.2); FP2 <- c(0.1, 0.5)

alpha1 <- alpha2 <- 6
as1 <- rep(0, 2)
as2 <- c(-2, pi / 2)
sigs1 <- sigs2 <- rep(40, 2)
oms1 <- c(0.6, 0.2)
oms2 <- c(0.1, 0.5)




# Global Functions -----------------------------------------------------------

# Main function (1d)
exp_sin_function <- function(x, a, sigma, omega){
  exp(-(x - a)^2 / sigma) * sin(omega * (x -a))
}

# Main derivative (1d)
exp_sin_derivative <- function(x, a, sigma, omega){
  (-2 * (x - a) / sigma * exp_sin_function(x, a, sigma, omega)
   + omega * exp(-(x - a)^2 / sigma) * cos(omega * (x -a)))
}

# main second derivative (1d)
exp_sin_second_derivative <- function(x, a, sigma, omega){
  -2 / sigma * exp_sin_function(x, a, sigma, omega) - 
  2 * (x - a) / sigma * exp_sin_derivative(x, a, sigma, omega) - 
  2 * omega * (x - a) / sigma * exp(-(x - a)^2 / sigma) * cos(omega * (x -a))- 
  omega^2 * exp_sin_function(x, a, sigma, omega)
}


#HEssian of covariate field

hessian <- function(x, beta){
  beta[1]*matrix(c(alpha1^2 * exp_sin_second_derivative(x[1], a = as1[1], sigma = sigs1[1], omega = oms1[1]) * exp_sin_function(x[2], as1[2], sigs1[2], oms1[2]),
           alpha1^2 * exp_sin_derivative(x[1], a = as1[1], sigma = sigs1[1], omega = oms1[1]) * exp_sin_derivative(x[2], as1[2], sigs1[2], oms1[2]),
           alpha1^2 * exp_sin_derivative(x[1], a = as1[1], sigma = sigs1[1], omega = oms1[1]) * exp_sin_derivative(x[2], as1[2], sigs1[2], oms1[2]),
           alpha1^2 * exp_sin_function(x[1], a = as1[1], sigma = sigs1[1], omega = oms1[1]) * exp_sin_second_derivative(x[2], as1[2], sigs1[2], oms1[2])),
         ncol = 2) +
    beta[2]*matrix(c(alpha2^2 * exp_sin_second_derivative(x[1], a = as2[1], sigma = sigs2[1], omega = oms2[1]) * exp_sin_function(x[2], as2[2], sigs2[2], oms2[2]),
             alpha2^2 * exp_sin_derivative(x[1], a = as2[1], sigma = sigs2[1], omega = oms2[1]) * exp_sin_derivative(x[2], as2[2], sigs2[2], oms2[2]),
             alpha2^2 * exp_sin_derivative(x[1], a = as2[1], sigma = sigs2[1], omega = oms2[1]) * exp_sin_derivative(x[2], as2[2], sigs2[2], oms2[2]),
             alpha2^2 * exp_sin_function(x[1], a = as2[1], sigma = sigs2[1], omega = oms2[1]) * exp_sin_second_derivative(x[2], as2[2], sigs2[2], oms2[2])),
           ncol = 2) -
    beta[3]*diag(2,2,2)
  
}



# Main gradient (2d)
exp_sin_gradient <- function(x, as, sigmas, omegas, alpha){
  if(!(all(c(length(x), length(as), length(sigmas), length(omegas)) == 2))){
    stop("x, as, sigmas, omegas must all be of length 2")
  }
  alpha * exp_sin_function(x, as, sigmas, omegas)[2:1] * 
    exp_sin_derivative(x, as, sigmas, omegas)
}



# Functions for covariates -----------------------------------------

fun_cov1 <- function(x) 
  alpha1 * prod(exp_sin_function(x, a = as1, sigma = sigs1, omega = oms1))

fun_cov2 <- function(x) 
  alpha2 * prod(exp_sin_function(x, a = as2, sigma = sigs2, omega = oms2))

fun_cov3 <- function(x) - sum(x^2)

fun_cov_list <- list(fun_cov1, fun_cov2, fun_cov3)


# Functions for covariates gradients -------------------------------

gradient_cov1 <- function(x)
  exp_sin_gradient(x, as = as1, sigmas = sigs1, omegas = oms1, alpha = alpha1)

gradient_cov2 <- function(x)
  exp_sin_gradient(x, as = as2, sigmas = sigs2, omegas = oms2, alpha = alpha2)

gradient_cov3 <- function(x) - 2 * x # derivative of - ||x||^2

gradient_list = list(gradient_cov1, gradient_cov2, gradient_cov3)


grad <- function(x, beta){
  beta[1]*gradient_cov1(x) +  beta[2]*gradient_cov2(x) +  beta[3]*gradient_cov3(x)
}


# Movement parameters

# parameters should be in the same order as the gradient list
beta_true <- c(beta1 = -1, beta2 = 0.5, beta3 = 0.05)




#####################
# Generating Tracks #
#####################

#max time for track
Tmax <- 5000
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 1


# Time grids
alltimes <- list()
for(i in 1:ntrack)
  alltimes[[i]] <- time

# Generate tracks -- This is very computational and may take a few hours.
alldat <- lapply(alltimes, function(times) {
  simLangevinMM(beta = beta_true, gamma2 = 1,
                times = times, loc0 = c(0, 0),
                grad_fun = gradient_list)
})

#changed the start position from runif(2, -20, 20) to c(0,0)



# Add ID column
for(zoo in 1:ntrack)
  alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])



############
# Thinning #
############

#parameters for thinning
thin = 250
#divided by six because of 5 extra points

X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]




##############
# Likelihood #
##############



#likelihood using extended kalman filter
#assuming R = 0

lik <- function(par, delta, X, grad){
  
  func <- function(x){
    return(par[1]*fun_cov1(x) + par[2]*fun_cov2(x) + par[3]*fun_cov3(x))
  }

  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #initial covariance guess
    P = 10*Q
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) {
      #control vector
      u = grad(z, par[1:3]) 
      
      
      F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, par[1:3])
      
      #predicted state estimate
      z_p = z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation covariance
      S = P 
      
      #updated state estimate
      z =  X[i, ]
      
      #updated estimate covariance
      P = diag(0,2,2)
      
      #adding likelihood contribution of i-th state
      l = l + dmvnorm(c(X[i, ] - z_p), mean = c(0,0), sigma = S, log = T)
      
      
      for (k in 1:m) {
        #control vector
        u = grad(z, par[1:3]) 
        
        
        F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, par[1:3])
        
        #predicted state estimate
        z_p = z + B %*% u
        
        #predicted estimate covariance 
        P = F_k %*% P %*% t(F_k) + Q
        
        #innovation covariance
        S = P 
        
        #updated state estimate
        z = z_p 
        
        #updated estimate covariance
        P = P
        
      }
      
    }
  }
  min(-l, 1e7)
  
}


##############
# Estimation #
##############

m = 1
delta = dt*thin/(m+1)

t1 = Sys.time()
cl <- makeCluster(4)     # set the number of processor cores
clusterExport(cl, varlist = c("lik", "grad", "hessian", "dmvnorm", "gradient_cov1", "gradient_cov1", "gradient_cov2", "gradient_cov3", "exp_sin_gradient", "exp_sin_derivative", "exp_sin_function", "exp_sin_second_derivative",
                             "IP1", "DP1", "CP1", "FP1", "IP2", "DP2", "CP2", "FP2", "alpha1", "as1", "sigs1", "oms1", "alpha2", "as2", "sigs2", "oms2", "m"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o1 = optimParallel(par=c(0,0,0,1), fn=lik, delta = delta, X = X, grad = grad, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1


m = 10
delta = dt*thin/(m+1)

t1 = Sys.time()
cl <- makeCluster(4)     # set the number of processor cores
clusterExport(cl, varlist = c("lik", "grad", "hessian", "dmvnorm", "gradient_cov1", "gradient_cov1", "gradient_cov2", "gradient_cov3", "exp_sin_gradient", "exp_sin_derivative", "exp_sin_function", "exp_sin_second_derivative",
                              "IP1", "DP1", "CP1", "FP1", "IP2", "DP2", "CP2", "FP2", "alpha1", "as1", "sigs1", "oms1", "alpha2", "as2", "sigs2", "oms2", "m"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o2 = optimParallel(par=c(0,0,0,1), fn=lik, delta = delta, X = X, grad = grad, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1


m = 500
delta = dt*thin/(m+1)

t1 = Sys.time()
cl <- makeCluster(12)     # set the number of processor cores
clusterExport(cl, varlist = c("lik", "grad", "hessian", "dmvnorm", "gradient_cov1", "gradient_cov1", "gradient_cov2", "gradient_cov3", "exp_sin_gradient", "exp_sin_derivative", "exp_sin_function", "exp_sin_second_derivative",
                              "IP1", "DP1", "CP1", "FP1", "IP2", "DP2", "CP2", "FP2", "alpha1", "as1", "sigs1", "oms1", "alpha2", "as2", "sigs2", "oms2", "m", "fun_cov1", "fun_cov2", "fun_cov3"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o3 = optimParallel(par=c(0,0,0,1), fn=lik, delta = delta, X = X, grad = grad, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o3



#testing hessian function
func <- function(x){
  return(fun_cov1(x)+fun_cov2(x)+fun_cov3(x))
}
par = c(1,1,1)
numDeriv::hessian(func, x)

hessian(x, c(1,1,1))



#estimate using thinned data
thin = 250
X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
N = nrow(X)
X = X[(0:(N%/%thin -1))*thin +1, ]
gradArray <- array(rep(0, dim(X)[1]*2*3), c(dim(X)[1], 2, 3))

for (i in 1:(dim(X)[1])) {
  gradArray[i,1:2, 1] = gradient_cov1(X[i, ])
  gradArray[i,1:2, 2] = gradient_cov2(X[i, ])
  gradArray[i,1:2, 3] = gradient_cov3(X[i, ])
}
locs = X
times = (alldat[[1]]$t)[(0:(N%/%thin -1))*thin +1]
ID = (alldat[[1]]$ID)[(0:(N%/%thin -1))*thin +1]
fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat




#estimate using full data
X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
gradArray <- array(rep(0, dim(X)[1]*2*3), c(dim(X)[1], 2, 3))

for (i in 1:(dim(X)[1])) {
  gradArray[i,1:2, 1] = gradient_cov1(X[i, ])
  gradArray[i,1:2, 2] = gradient_cov2(X[i, ])
  gradArray[i,1:2, 3] = gradient_cov3(X[i, ])
}

locs = X
times = alldat[[1]]$t
ID = alldat[[1]]$ID
fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat



















