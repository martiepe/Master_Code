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
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = 0.1, 
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



UD_value <- function(s, UD){
  x = s[1]
  y = s[2]
  
  
  x1 = floor(s[1])
  x2 = ceiling(s[1])
  y1 = floor(s[2])
  y2 = ceiling(s[2])
  
  if(x1 == x2 && y1 == y2){
    UD[x1, y1]
  }
  
  
  f11 = UD[x1+lim[2]+1, y1+lim[2]+1]
  f12 = UD[x1+lim[2]+1, y2+lim[2]+1]
  f21 = UD[x2+lim[2]+1, y1+lim[2]+1]
  f22 = UD[x2+lim[2]+1, y2+lim[2]+1]
  
  #value of UD at location
  ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
  
}





##############
# Likelihood #
##############

#using brownian bridge
delta = dt*thin
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 2:(nrow(X)-1)) {
    L = 0
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

#without importance sampling

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
N=4
# Your function
lik1 <- function(beta) {
  lik(c(beta, 2, -0.5, 5))
}

# Grid of values
beta_grid <- seq(1,7,0.1)

# Set up a cluster
cl <- makeCluster(12)

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









