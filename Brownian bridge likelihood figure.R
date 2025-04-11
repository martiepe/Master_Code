library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(mvtnorm)
library(ambient)

set.seed(123)

#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 1
covlist <- list()

#perlin covariates
covlist <- list()
xgrid <- seq(lim[1], lim[2], by = resol)
ygrid <- seq(lim[3], lim[4], by = resol)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}
# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[ncov+1]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))




###################
## Simulate data ##
###################
#max time for track
Tmax <- 5000
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 1
#speed parameter for Langevin model
speed <- 5


# Generate tracks
beta = c(4,-0.1)

X = simLangevinMM(beta = beta, gamma2 = speed, times = time, loc0 = c(0, 0), cov_list = covlist)


############
# Thinning #
############

#parameters for thinning
thin = 10
#divided by six because of 5 extra points

X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]


##############
# Likelihood #
##############
delta = dt*thin

# Euler-Maruyama likelihood
EM_lik <- function(par){
  grad = bilinearGradArray(X, covlist) 
  u = par[1]*grad[,,1] + par[2]*grad[,,2]
  N = nrow(X)
  l = 0
  for (i in 1:(N-1)) {
    l = l + dmvnorm(X[i+1, ] - X[i, ] - par[4]*delta*u[i, ]/2, mean = c(0,0), sigma = diag(delta*par[4],2,2), log = TRUE)
  }
  return(-l)
}



#importance sampling using brownian bridge
lik <- function(par, N, M){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 2:(nrow(X)-1)) {
    L = 0
    #generating M brownian bridges
    for (j in 1:M) {
      L_k = 1
      #generating nodes
      sigma = matrix(nrow = N, ncol = N)
      
      mu_x = c()
      mu_y = c()
      #sigma is the same as long as delta is the same, so it can be computed once
      for (k in 1:N) {
        for (m in 1:k) {
          sigma[k,m] = delta*(1 - k/(N+1))*(m/(N+1))
          sigma[m,k] = delta*(1 - k/(N+1))*(m/(N+1))
        }
        mu_x = c(mu_x, X[i, 1] + k*(X[i+1, 1] - X[i, 1])/(N+1))
        mu_y = c(mu_y, X[i, 2] + k*(X[i+1, 2] - X[i, 2])/(N+1))
      }
      
      x = rmvnorm(1, mean = mu_x, sigma = par[4]*sigma)
      y = rmvnorm(1, mean = mu_y, sigma = par[4]*sigma)
      
      L_k = L_k/dmvnorm(x, mean = mu_x, sigma = par[4]*sigma)
      L_k = L_k/dmvnorm(y, mean = mu_y, sigma = par[4]*sigma)
      
      
      #
      grad = bilinearGrad(X[i, ], covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2])/(2*(N+1))
      L_k = L_k*dmvnorm(X[i, ] , mean = c(x[1], y[1]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      
      for (k in 1:(N-1)) {
        grad = bilinearGrad(c(x[k], y[k]), covlist)
        u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2])/(2*(N+1))
        L_k = L_k*dmvnorm(c(x[k+1], y[k+1]) , mean = c(x[k], y[k]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      }
      grad = bilinearGrad(c(x[N], y[N]), covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2])/(2*(N+1))
      L_k = L_k*dmvnorm(X[i+1, ] , mean = c(x[N], y[N]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      
      L = L + L_k/M
    }
    
    l = l + log(L)
    
  }
  
  return(-l)
  
}

X

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "delta", "covlist", "rmvnorm", "dmvnorm", "bilinearGrad", "EM_lik", "lik"))

#M = 50
#N=9 intermediate steps
compute_likelihood <- function(beta1){
  return(lik(c(beta1,-0.1, 5), 9, 50))
}
results1 <- parLapply(cl, seq(-8, 8, 0.5), compute_likelihood)
#N = 4intermediate steps
compute_likelihood <- function(beta1){
  return(lik(c(beta1, -0.1, 5), 4, 50))
}
results2 <- parLapply(cl, seq(-8, 8, 0.5), compute_likelihood)


#M = 20
#N=9 intermediate steps
compute_likelihood <- function(beta1){
  return(lik(c(beta1, -0.1, 5), 9, 20))
}
results3 <- parLapply(cl, seq(-8, 8, 0.5), compute_likelihood)
#N = 4intermediate steps
compute_likelihood <- function(beta1){
  return(lik(c(beta1, -0.1, 5), 4, 20))
}
results4 <- parLapply(cl, seq(-8, 8, 0.5), compute_likelihood)


#Euler-Maruyama likelihood
compute_likelihood <- function(beta1){
  return(lik(c(beta1, -0.1, 5), 4, 50))
}
resultsEM <- parLapply(cl, seq(-8, 8, 0.5), compute_likelihood)

stopCluster(cl)


ggplot()+
  geom_line(aes(seq(-8, 8, 0.5), results1)) +
  geom_line(aes(seq(-8, 8, 0.5), results2)) +
  geom_line(aes(seq(-8, 8, 0.5), results3)) +
  geom_line(aes(seq(-8, 8, 0.5), results4)) +
  geom_line(aes(seq(-8, 8, 0.5), resultsEM))


