# Load packages
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
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



#likelihood using extended kalman filter
#assuming R = 0

grad = bilinearGradArray(X, covlist)

lik <- function(par){
  #log-likelihood
  l = 0
  
  U = Delta*par[4]*(grad[,,1]*par[1] + grad[,,2]*par[2] + grad[,,3]*par[3])/2
  UD = getUD(covariates=covlist, beta=par[1:3])$z
  
  for (j in 1:1) {
    for (i in 2:(nrow(X)-1)) {
      #control vector
      u1 = U[i, ]
      u2 = U[i+1, ]

      a = min(1, UD_value(X[i+1, ], UD)*dmvnorm(X[i, ], mean = (X[i+1, ] + u2), sigma = diag(par[4]*Delta, 2, 2))/
                (UD_value(X[i, ], UD)*dmvnorm(X[i+1, ], mean = (X[i, ] + u1), sigma = diag(par[4]*Delta, 2, 2))))
      
      l = l + log(a) + dmvnorm(X[i, ], mean = (X[i+1, ] + u2), sigma = diag(par[4]*Delta, 2, 2), log = T)
      #if(i < 2){
      #  print(dmvnorm(X[i+1, ], mean = (X[i, ] + u1), sigma = diag(par[4]*Delta, 2, 2)))
      #  print((UD_value(X[i, ], UD)))
      #}
      
    }
  }
  return(-l)
  
}
  

Delta = dt*thin



t1 = Sys.time()
lik(c(0,0,1,1), delta, X, grad)
Sys.time() - t1


t1 = Sys.time()
optim(c(1,1,1,1), lik)
Sys.time() - t1

library(optimParallel)
t1 = Sys.time()
cl <- makeCluster(5)     # set the number of processor cores
clusterExport(cl, varlist = c("lik", "grad", "dmvnorm", "X", "UD_value", "getUD", "Delta", "covlist", "lim"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o2 = optimParallel(par=c(0,0,0,1), fn=lik, delta = delta, X = X, grad = grad, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o2


beta
lik(c(beta, 5))
lik(c(0.003425801,  0.001479025, -0.034209244,  5.008598236))
















