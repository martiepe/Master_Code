
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




#perlin covariates
library(ambient)
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
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))








###################
## Simulate data ##
###################
beta <- c(4,2,-0.1)
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
#covariate coefficients
beta <- c(4,2,-0.1)
# Time grids
alltimes <- list()
for(i in 1:ntrack)
  alltimes[[i]] <- time

# Generate tracks -- This is very computational and may take a few hours.
alldat <- lapply(alltimes, function(times) {
  simLangevinMM(beta = beta, gamma2 = speed, times = times,
                loc0 = c(0, 0), cov_list = covlist)
})

#changed the start position from runif(2, -20, 20) to c(0,0)



# Add ID column
for(zoo in 1:ntrack)
  alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])





#estimate using thinned data
thin = 10
X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
N = nrow(X)
X = X[(0:(N%/%thin -1))*thin +1, ]
#X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)

gradArray = bilinearGradArray(X, covlist)

locs = X
times = alldat[[1]]$t[(0:(N%/%thin -1))*thin +1]
ID = alldat[[1]]$ID[(0:(N%/%thin -1))*thin +1]

fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat









###### one brownian bridge importance sampling estimate #####
N = thin-1
M = 200
delta = dt*thin
#find the diffusion constant to be used
gradArray = bilinearGradArray(X, covlist)
locs = X
gammasq = langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)$gamma2Hat

#simulate brownian bridges
sigma = matrix(nrow = N, ncol = N)
for (k in 1:N) {
  for (m in 1:k) {
    sigma[k,m] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
    sigma[m,k] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
  }
}

B <- array(data = NA, c(2, nrow(X)-1, M, N))
bridge <- array(data = NA, c(2, M, N))
#probability of each bridge
P <- c()
for (j in 1:M) {
  bridge[1, j, 1:N] = rmvnorm(1, mean = rep(0, N), sigma = sigma)
  bridge[2, j, 1:N] = rmvnorm(1, mean = rep(0, N), sigma = sigma)
  
  P = c(P, dmvnorm(bridge[1, j, 1:N], rep(0, N), sigma)*dmvnorm(bridge[2, j, 1:N], rep(0, N), sigma))
}


for (k in 1:(nrow(X)-1)) {
  mu_x = c()
  mu_y = c()
  for (i in 1:N) {
    mu_x = c(mu_x, X[k, 1] + i*(X[k+1, 1] - X[k, 1])/(N+1))
    mu_y = c(mu_y, X[k, 2] + i*(X[k+1, 2] - X[k, 2])/(N+1))
  }
  for (j in 1:M) {
    B[1, k, j, 1:N] = mu_x + bridge[1, j, 1:N]
    B[2, k, j, 1:N] = mu_y + bridge[2, j, 1:N]
  }
}


#find the gradient at the bridge nodes
Grad <- array(data = NA, c(ncov +1,nrow(X)-1, M, N, 2))
t1 <- Sys.time()
for (k in 1:(nrow(X)-1)) {
  for(i in 1:M){
    grad = bilinearGradArray(t(B[1:2, k, i, 1:N]), covlist)
    for (j in 1:(ncov+1)) {
      Grad[j, k, i, 1:N, 1] = grad[,1,j]
      Grad[j, k, i, 1:N, 2] = grad[,2,j]
    }
  }
}
Sys.time()- t1






#vectorized likelihood
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 1:(nrow(X)-1)) {
    L_k = 1/P
    
    # Add endpoints to all samples (M x (N+2) matrices)
    # calc initial gradient
    x_samples = B[1, i, , ]
    y_samples = B[2, i, , ]
    
    grad_0 <- array(data = t(gradArray[i, , ]), c(3,1,2))
    
    u_0 <- (delta*par[4]/((N+1)*2)) * 
      (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
    
    full_x <- cbind(X[i,1], x_samples, X[i+1,1])
    full_y <- cbind(X[i,2], y_samples, X[i+1,2])
    # likelihood of all locations
    L_k <- sapply(seq(M), function(j) {
      grads <- Grad[ , i, j, , ]
      us <- (delta*par[4]/((N+1)*2)) * 
        (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
      us <- rbind(u_0, us)
      prod(dmvn((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                   cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us, 
                matrix(c(0,0)),
                diag(delta*par[4]/(N+1), 2, 2)))
    })*L_k
    
    
    l = l + log(sum(L_k/M))
    
  }
  
  return(-l)
  
}


t1 = Sys.time()
lik(c(4,2,-0.1,5))
Sys.time() - t1

library(mvnfast)
library(optimParallel)

t1 = Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "P", "B", "Grad", "gradArray", "dmvn"))
setDefaultCluster(cl=cl)
o = optimParallel(c(0,0,0,1), lik)
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o
