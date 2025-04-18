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
set.seed(123)

#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 2
covlist <- list()
#simulate spatial covariates wuing grf with matern covariance function
for(i in 1:ncov) {
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 1.5, rho = 50, sigma2 = 1, 
                                resol = resol, raster_like = TRUE)
}

# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))



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






# Compute utilisation distribution
beta <- c(4,2,-0.1)
UD <- getUD(covariates=covlist, beta=beta)

# Plot covariates
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
c3plot <- plotRaster(rhabitToRaster(covlist[[3]]), scale.name = expression(c[2])) + ggtheme
UD <- getUD(covariates=covlist, beta=beta)
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

UDplot

c2plot



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

############
# Thinning #
############

#parameters for thinning
thin = 5
#divided by six because of 5 extra points

X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
times = alldat[[1]]$t[(0:(nrow(X)%/%thin -1))*thin +1]
ID = alldat[[1]]$ID[(0:(nrow(X)%/%thin -1))*thin +1]
X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]



############################################
## checking if covariates and thin is ok ###
############################################

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




##############
# Likelihood #
##############

# Euler-Maruyama likelihood
EM_lik <- function(par){
  grad = bilinearGradArray(X, covlist) 
  u = par[1]*grad[,,1] + par[2]*grad[,,2] + par[3]*grad[,,3]
  N = nrow(X)
  l = 0
  for (i in 1:(N-1)) {
    l = l + dmvnorm(X[i+1, ] - X[i, ] - par[4]*delta*u[i, ]/2, mean = c(0,0), sigma = diag(delta*par[4],2,2), log = TRUE)
  }
  return(-l)
}


#importance sampling using brownian bridge
delta = dt*thin
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 2:(nrow(X)-1)) {
    L = 0
    
    #covaance and mean of brownian bridge
    sigma = matrix(nrow = N, ncol = N)
      
    mu_x = c()
    mu_y = c()
    for (k in 1:N) {
      for (m in 1:k) {
        sigma[k,m] = delta*(1 - k/(N+1))*(m/(N+1))
        sigma[m,k] = delta*(1 - k/(N+1))*(m/(N+1))
      }
      mu_x = c(mu_x, X[i, 1] + k*(X[i+1, 1] - X[i, 1])/(N+1))
      mu_y = c(mu_y, X[i, 2] + k*(X[i+1, 2] - X[i, 2])/(N+1))
    }
    
    for (j in 1:M) {
      L_k = 1
      #simulatig brownian bridge
      x = rmvnorm(1, mean = mu_x, sigma = par[4]*sigma)
      y = rmvnorm(1, mean = mu_y, sigma = par[4]*sigma)
      #
      L_k = L_k/dmvnorm(x, mean = mu_x, sigma = par[4]*sigma)
      L_k = L_k/dmvnorm(y, mean = mu_y, sigma = par[4]*sigma)

      
      #adding Langevin process transition probability of proposed brownian bridge
      grad = bilinearGrad(X[i, ], covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
      L_k = L_k*dmvnorm(X[i, ] , mean = c(x[1], y[1]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      
      for (k in 1:(N-1)) {
        grad = bilinearGrad(c(x[k], y[k]), covlist)
        u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
        L_k = L_k*dmvnorm(c(x[k+1], y[k+1]) , mean = c(x[k], y[k]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      }
      grad = bilinearGrad(c(x[N], y[N]), covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
      L_k = L_k*dmvnorm(X[i+1, ] , mean = c(x[N], y[N]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      
      L = L + L_k/M
    }
    
    l = l + log(L)
    
  }
  
  return(-l)
  
}

#parralellized importance sampling using brownian bridge
delta = dt*thin
lik <- function(par, cl) {
  
  
  #compute likelihood of one observation
  compute_L_i <- function(i) {
    L = 0
    
    #making mean and covariance for brownian bridge
    sigma = matrix(nrow = N, ncol = N)
    mu_x = c()
    mu_y = c()
    for (k in 1:N) {
      for (m in 1:k) {
        sigma[k,m] = delta * (1 - k/(N+1)) * (m/(N+1))
        sigma[m,k] = sigma[k,m]
      }
      mu_x = c(mu_x, X[i, 1] + k * (X[i+1, 1] - X[i, 1]) / (N+1))
      mu_y = c(mu_y, X[i, 2] + k * (X[i+1, 2] - X[i, 2]) / (N+1))
      
    }
    
    for (j in 1:M) {
      #simulating brownian bridge
      x = rmvnorm(1, mean = mu_x, sigma = par[4] * sigma)
      y = rmvnorm(1, mean = mu_y, sigma = par[4] * sigma)
      
      #likelihood of brownian bridge
      L_k = 1
      L_k = L_k / dmvnorm(x, mean = mu_x, sigma = par[4] * sigma)
      L_k = L_k / dmvnorm(y, mean = mu_y, sigma = par[4] * sigma)
      
      #adding Langevin process transition probability of proposed brownian bridge
      grad = bilinearGrad(X[i, ], covlist)
      u = delta * par[4] * (grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3]) / (2 * (N+1))
      L_k = L_k * dmvnorm(X[i, ], mean = c(x[1], y[1]) + u, sigma = diag(delta * par[4] / (N+1), 2, 2))
      
      for (k in 1:(N-1)) {
        grad = bilinearGrad(c(x[k], y[k]), covlist)
        u = delta * par[4] * (grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3]) / (2 * (N+1))
        L_k = L_k * dmvnorm(c(x[k+1], y[k+1]), mean = c(x[k], y[k]) + u, sigma = diag(delta * par[4] / (N+1), 2, 2))
      }
      
      grad = bilinearGrad(c(x[N], y[N]), covlist)
      u = delta * par[4] * (grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3]) / (2 * (N+1))
      L_k = L_k * dmvnorm(X[i+1, ], mean = c(x[N], y[N]) + u, sigma = diag(delta * par[4] / (N+1), 2, 2))
      
      L = L + L_k / M
    }
    
    return(log(L))
  }
  
  results <- parLapply(cl, 2:(nrow(X)-1), compute_L_i)
  L<- sum(unlist(results))
  
  return(-L)
}
 

#monte carlo likelihood
delta = dt*thin
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  for (i in 2:(nrow(X)-1)) {
    #for each observation in the track
    L = 0
      for (j in 1:M) {
      L_k = 1
      x = X[i, ]
      
      for (k in 1:N) {
        grad = bilinearGrad(x, covlist)
        u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
        
        x_p = rmvnorm(1, mean = x + u, sigma = diag(par[4]*delta/(N+1), 2, 2))
        
        L_k = L_k*dmvnorm(x_p ,  mean = x + u, sigma = diag(par[4]*delta/(N+1), 2, 2))
        
        x = x_p
      }
      grad = bilinearGrad(x, covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
      
      L_k = L_k*dmvnorm(X[i+1, ] , mean = x + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      L = L + L_k/M
      
      }
    
    l = l + log(L)
    
  }
  
  return(-l)
}



#parralellized monte carlo
delta = dt*thin
lik <- function(par, cl) {

  
  #likelihood of one observation
  compute_L_i <- function(i) {
    L = 0
    for (j in 1:M) {
      L_k = 1
      x = X[i, ]
      
      for (k in 1:N) {
        grad = bilinearGrad(x, covlist)
        u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
        
        x_p = rmvnorm(1, mean = x + u, sigma = diag(par[4]*delta/(N+1), 2, 2))
        
        L_k = L_k*dmvnorm(x_p ,  mean = x + u, sigma = diag(par[4]*delta/(N+1), 2, 2))
        
        x = x_p
      }
      grad = bilinearGrad(x, covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/(2*(N+1))
      
      L_k = L_k*dmvnorm(X[i+1, ] , mean = x + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      L = L + L_k/M
      
    }
    return(log(L))
  }
  
  # Run in parallel
  results <- parLapply(cl, 2:(nrow(X)-1), compute_L_i)
  L <- sum(unlist(results))
  
  return(-L)
}


###### one brownian bridge importance sampling estimate #####
N = thin-1
M = 100
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



#not vectorized
delta = dt*thin
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 1:(nrow(X)-1)) {
    L = 0
    for (j in 1:M) {
      L_k = 1
      #brownian bridge contribution to likelihood
      L_k = L_k/P[j]
      
      #adding Langevin process transition probability of proposed brownian bridge
      
      u = delta*par[4]*(gradArray[i]*par[1] + gradArray[i]*par[2] + gradArray[i]*par[3])/(2*(N+1))
      L_k = L_k*dmvnorm(c(B[1, i, j, 1], B[2, i, j, 1]) , mean = X[i, ] + u, sigma = diag(delta*par[4]/(N+1), 2, 2))

      
      for (k in 1:(N-1)) {
        u = delta*par[4]*(Grad[1, i, j, k, ]*par[1] + Grad[2, i, j, k, ]*par[2] + Grad[3, i, j, k, ]*par[3])/(2*(N+1))
        L_k = L_k*dmvnorm(c(B[1, i, j, k+1], B[2, i, j, k+1]) , mean = c(B[1, i, j, k], B[2, i, j, k]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      }
      
      u = delta*par[4]*(Grad[1, i, j, N, ]*par[1] + Grad[2, i, j, N, ]*par[2] + Grad[3, i, j, N, ]*par[3])/(2*(N+1))
      L_k = L_k*dmvnorm(X[i+1, ] , mean = c(B[1, i, j, N], B[2, i, j, N]) + u, sigma = diag(delta*par[4]/(N+1), 2, 2))
      
      
      L = L + L_k/M
    }
    
    l = l + log(L)
    
  }
  
  return(-l)
  
}

#vectorized
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









################

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "covlist", "rmvnorm", "dmvnorm", "bilinearGrad"))
lik(c(4,2,-0.1,5), cl)
stopCluster(cl)


results = c()
M = 50
N = thin-1

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "covlist", "rmvnorm", "dmvnorm", "bilinearGrad"))

results = c(results, lik(par = c(0, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(2, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(4, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(6, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(8, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(10, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(12, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(14, 2, -0.1, 5), cl = cl))
# Cleanup
stopCluster(cl)

EM = c()
EM = c(EM, EM_lik(par = c(0, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(2, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(4, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(6, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(8, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(10, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(12, 2, -0.1, 5)))
EM = c(EM, EM_lik(par = c(14, 2, -0.1, 5)))

ggplot()+
  geom_line(aes(c(0,2,4,6,8,10,12,14), results)) +
  geom_line(aes(c(0,2,4,6,8,10,12,14), EM), color = "red")



##########################
## Parameter Estimation ##
##########################


#number of simulations
M = 50
#number of nodes
N = thin-1

t1 = Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "covlist", "rmvnorm", "dmvnorm", "bilinearGrad"))
o = optim(c(0,0,0,1), lik, cl = cl, method = "Nelder-Mead", control = list(maxit = 200))
stopCluster(cl)
Sys.time() - t1
o




######################################
## testing lik for a grid of values ##
######################################

M = 40
N = 4
N=thin-1

lik1 <- function(beta) {
  lik(c(beta, 2, -0.5, 5))
}

# Grid of values
beta_grid <- seq(1,7,0.1)

# Set up a cluster
cl <- makeCluster(detectCores()-1)

# Export functions and variables
clusterExport(cl, varlist = c("lik1", "N", "M", "rmvnorm", "dmvnorm", "bilinearGrad", "delta", "lik", "X", "covlist"))

# Run in parallel
results <- parLapply(cl, beta_grid, lik1)

# Shut down cluster
stopCluster(cl)

# Convert list to vector
results <- unlist(results)

print(results)


ggplot() +
  geom_line(aes(beta_grid, results)) +
  labs(x="beta", y = "-l", title = "brownian bridge likelihood") +
  theme_bw()









