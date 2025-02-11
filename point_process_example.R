

# Load packages
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
set.seed(1)



#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*150
resol <- 1
ncov <- 2
covlist <- list()
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




##############################
## Simulate occurrence data ##
##############################

lambda = 0.01*exp(beta[1]*covlist[[1]]$z + beta[2]*covlist[[2]]$z + beta[3]*covlist[[3]]$z)



#intensity of study area
lambda_max = max(lambda)



#simulating number of points
n = rpois(1, 4*lim[2]*lim[2]*lambda_max)

#simulating points
s = matrix(runif(2*n, min = -100, max = 100), ncol = 2)

#proper number of points
N = 0
#proper locations
Y = c()


#trimming number of points to match inhomogeneous poisson process
for (i in 1:n) {
  x = s[i, 1]
  y = s[i, 2]
  
  
  x1 = floor(s[i, 1])
  x2 = ceiling(s[i, 1])
  y1 = floor(s[i, 2])
  y2 = ceiling(s[i, 2])
  
  f11 = lambda[x1+101, y1+101]
  f12 = lambda[x1+101, y2+101]
  f21 = lambda[x2+101, y1+101]
  f22 = lambda[x2+101, y2+101]
  
  #intensity at location
  lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
  
  #acceptance rate
  p = lambda_s/lambda_max
  
  if(p >= runif(1)){
    N = N + 1
    Y = c(Y, s[i, ])
  }
}

Y = matrix(Y, ncol = 2)

N
n


#plotting simulations
ggplot() +
  geom_point(aes(S[,1],S[,2]), colour = "grey80")+
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))









################
## estimation ##
################

#likelihood function
lik <- function(X, Y, B_1, B_2, B_3, kappa, ntrack, Delta){
  #log-likelihood
  l = 0

  #covariate field used to calculate intensity
  cov_field = exp(B_1*covlist[[1]]$z + B_2*covlist[[2]]$z + B_3*covlist[[3]]$z)
  
  #finding intensity of study area
  lambda_sum = 0
  for (i in -lim[2]:(lim[2]-1)) {
    for (j in -lim[2]:(lim[2]-1)) {
      x1 = i
      x2 = i+1
      y1 = j
      y2 = j+1
      
      f11 = kappa*(cov_field[x1+lim[2]+1, y1+lim[2]+1])
      f12 = kappa*(cov_field[x1+lim[2]+1, y2+lim[2]+1])
      f21 = kappa*(cov_field[x2+lim[2]+1, y1+lim[2]+1])
      f22 = kappa*(cov_field[x2+lim[2]+1, y2+lim[2]+1])
      
      
      lambda_sum = lambda_sum + ((y2+y1)*((x2+x1)*f11 + (x2-3*x1)*f21) + (y2-3*y1)*((x2+x1)*f12 + (x2-3*x1)*f22))/4
    }
  }
  l = l - lambda_sum
  
  #intensity of locations
  for (i in 1:nrow(Y)) {
    x = Y[i, 1]
    y = Y[i, 2]
    
    
    x1 = floor(Y[i, 1])
    x2 = ceiling(Y[i, 1])
    y1 = floor(Y[i, 2])
    y2 = ceiling(Y[i, 2])
    
    f11 = kappa*cov_field[x1+lim[2]+1, y1+lim[2]+1]
    f12 = kappa*cov_field[x1+lim[2]+1, y2+lim[2]+1]
    f21 = kappa*cov_field[x2+lim[2]+1, y1+lim[2]+1]
    f22 = kappa*cov_field[x2+lim[2]+1, y2+lim[2]+1]
    
    #intensity at location
    lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
    
    l = l + log(lambda_s)
  }
  
  
  return(l)
}


#par: B1, B2, B3, kappa
par = c(0,0,0,1)

o = optim(lik, par)

