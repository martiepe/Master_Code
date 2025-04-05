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
  vals = 2*noise_perlin(c(length(xgrid), length(ygrid)))
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
  # Compute the list of i indices to parallelize
  i_list <- 2:(nrow(X)-1)
  
  # Function to compute log-likelihood contribution for a single i
  compute_L_i <- function(i) {
    L = 0
    sigma = matrix(nrow = N, ncol = N)
    for (k in 1:N) {
      for (m in 1:k) {
        sigma[k,m] = delta * (1 - k/(N+1)) * (m/(N+1))
        sigma[m,k] = sigma[k,m]
      }
    }
    
    for (j in 1:M) {
      mu_x = mu_y = numeric(N)
      for (k in 1:N) {
        mu_x[k] = X[i, 1] + k * (X[i+1, 1] - X[i, 1]) / (N+1)
        mu_y[k] = X[i, 2] + k * (X[i+1, 2] - X[i, 2]) / (N+1)
      }
      
      x = rmvnorm(1, mean = mu_x, sigma = par[4] * sigma)
      y = rmvnorm(1, mean = mu_y, sigma = par[4] * sigma)
      
      L_k = 1
      L_k = L_k / dmvnorm(x, mean = mu_x, sigma = par[4] * sigma)
      L_k = L_k / dmvnorm(y, mean = mu_y, sigma = par[4] * sigma)
      
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
  
  # Run in parallel
  results <- parLapply(cl, i_list, compute_L_i)
  total_log_likelihood <- sum(unlist(results))
  
  return(-total_log_likelihood)
}

results = c()
M = 60
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











#without importance sampling
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
lik(c(4,2,-0.1,5))


#using modified Brownian bridge


##########################
## Parameter Estimation ##
##########################


#number of simulations
M = 40
#number of nodes
N = 4

t1 = Sys.time()
lik(c(0,0,0,1))
Sys.time() - t1

t1 = Sys.time()
o = optim(c(0,0,0,1), lik, method = "Nelder-Mead")
Sys.time() - t1
o


lik(c(4,2,-0.1,5))

######################################
## testing lik for a grid of values ##
######################################

M = 40
<<<<<<< HEAD
N = 4
=======
N=4
>>>>>>> 0ecf6646ff8a42f2fbbdb2c9cd784bd1bc9dc83c
# Your function
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









