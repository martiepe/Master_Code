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


lik <- function(par){
  #log-likelihood
  l = 0

  for (j in 1:1) {
    for (i in 2:(nrow(X)-1)) {
      #control vector
      z = X[i, ]
      for (k in 1:N) {
        grad = bilinearGrad(z, covlist)
        u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/2
        z_p = z + u*delta*par[4]/2 + rmvnorm(1, mean = c(0,0), sigma = diag(delta*par[4], 2, 2))
        l = l + dmvnorm(z , mean = z_p, sigma = diag(delta*par[4], 2, 2))
        z = z_p
      }
      grad = bilinearGrad(z, covlist)
      u = delta*par[4]*(grad[,1]*par[1] + grad[,2]*par[2] + grad[,3]*par[3])/2
      l = l + dmvnorm(X[i, ] , mean = z + u*delta*par[4]/2, sigma = diag(delta*par[4], 2, 2))
    }
  }
  
  
  return(-l)
  
}



N = 5
delta = dt*thin/(N+1)

t1 = Sys.time()
optim(c(0,0,0,1), lik, method = "Nelder-Mead", lower = c(-Inf, -Inf, -Inf, 0))
Sys.time() - t1


l = c()


t1 = Sys.time()
while (length(l) < 1000) {
  l = c(l, lik(c(4,2,-0.1, 5)))
}
Sys.time() - t1


var(l)
mean(l)

max(l)-min(l)







