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
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 2, rho = 50, sigma2 = 0.1, 
                                resol = resol, raster_like = TRUE)
}

# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))



# Compute utilisation distribution
beta <- c(4,2,-0.1)
UD <- getUD(covariates=covlist, beta=beta)




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



##############
# Likelihood #
##############

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


M = 40
N = 4

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "covlist", "rmvnorm", "dmvnorm", "bilinearGrad"))

results = c()
results = c(results, lik(par = c(4, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(2, 2, -0.1, 5), cl = cl))
results = c(results, lik(par = c(6, 2, -0.1, 5), cl = cl))
# Cleanup
stopCluster(cl)


ggplot()+
  geom_line(aes(c(2,4,6), results))









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









