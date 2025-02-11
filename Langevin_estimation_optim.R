library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(stats)
set.seed(1)





#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
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
UD <- getUD(covariates=covlist, beta=beta)






#########################
## Simulate track data ##
#########################
Tmax <- 5000
dt <- 0.01
time <- seq(0, Tmax, by = dt)
ntrack <- 10
speed <- 5

# Time grids
alltimes <- list()
for(i in 1:ntrack)
  alltimes[[i]] <- time

# Generate tracks -- This is very computational and may take a few hours.
alldat <- lapply(alltimes, function(times) {
  simLangevinMM(beta = beta, gamma2 = speed, times = times,
                loc0 = c(0,0), cov_list = covlist)
})

# Add ID column
for(zoo in 1:ntrack)
  alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])



#plotting trails

alldat[[1]]

length(alldat[[1]]$x)

ggplot()+
  geom_path(aes(alldat[[1]]$x[(0:2500)*200+1], alldat[[1]]$y[(0:2500)*200+1]),color = "grey") +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))


A = data.frame(alldat[[1]])



tail(A)

for (i in 2:100) {
  A = rbind(A, alldat[[i]])
}
length(A$ID)
ggplot()+
  geom_path(aes(A$x[(0:100000)*500+1], A$y[(0:100000)*500+1]), colour = A$ID[(0:100000)*500+1]) +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))





load("C:/Users/marti/OneDrive/Skrivebord/master/langevin/track simulation.RData")


#number of tracks used in estimation
ntrack = 10
#thinning factor used
thinning = 2
#step length
Delta = 0.01*thinning

#using the first n tracks
X = alldat[1:ntrack]
#thinning the tracks
for (i in 1:ntrack) {
  X[[i]] = matrix(c(alldat[[i]]$x[(0:(500000/thinning))*thinning + 1], alldat[[i]]$y[(0:(500000/thinning))*thinning + 1]), ncol = 2)
}






#finding gradients of tracks
grad1 = array(dim= c(dim(X[[1]])[1], 2, ntrack))
grad2 = array(dim= c(dim(X[[1]])[1], 2, ntrack))
grad3 = array(dim= c(dim(X[[1]])[1], 2, ntrack))

for (i in 1:ntrack) {
  locs = X[[i]]
  
  grad1[,,i] = bilinearGradArray(locs = locs, cov_list = list(covlist[[1]]))
  grad2[,,i] = bilinearGradArray(locs = locs, cov_list = list(covlist[[2]]))
  grad3[,,i] = bilinearGradArray(locs = locs, cov_list = list(covlist[[3]]))
  
}


#likelihood function
lik <- function(par){
  #log-likelihood
  l = 0
  #adding likelihood from track data
  for (i in 1:10) {
    locs = X[[i]]
    len = dim(locs)[1]
    
    l = l + sum(dmvnorm(locs[2:len,] - locs[1:(len-1),] - 0.5*par[4]*Delta*(par[1]*grad1[1:(len-1),,i] + par[2]*grad2[1:(len-1),,i] + par[3]*grad3[1:(len-1),,i]),mean = c(0,0), sigma = par[4]*Delta*diag(1, 2, 2), log = T))
  }
  return(-l)
}









pars = c(0,0,0,1)

optim(pars, lik)

