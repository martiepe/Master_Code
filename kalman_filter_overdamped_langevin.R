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



# Plot covariates
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

pdf("sim2cov.pdf", width=12, height=3)
grid.arrange(c1plot, c2plot, UDplot, nrow=1)
dev.off()



c1plot
c2plot
UDplot

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
ntrack <- 10
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



save.image(file='10Langevin_tracks.RData')

load('10Langevin_tracks.RData')





#################
## likelihood ###
#################

#par: B_1, B_2, B_3, gamma^2

lik <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    #observation matrix
    H = diag(1,2,2)
    #transition matrix
    F_k = diag(1,2,2)
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    #observation covariance
    R = diag(0,2,2)
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    
    #initial covariance guess
    P = 10*Q
    
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) {
      
      #control vector
      u = gradArray[i,,] %*% par[1:3]
      
      #predicted state estimate
      z_p = F_k %*% z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation
      y = X[i, ] - c(H %*% z_p)
      
      #innovation covariance
      S = H %*% P %*% t(H) 
      
      #optimal Kalman gain
      K = P %*% t(H) %*% solve(S)
      
      #updated state estimate
      z = z_p + K %*% y
      
      #updated estimate covariance
      P = (diag(1, 2, 2) - K %*% H) %*% P
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - H %*% z_p), mean = c(0,0), sigma = S, log = T)
      
      
    }
  }
  return(l)
  
}


#The likelihood is giving the using the negative value of Beta. Why?


par = c(-3.72752931, -1.81252269,  0.08142107,  5.01331337)

par2 = c(3.72752931, 1.81252269,  -0.08142107,  5.01331337)


lik(par2)





#parameters for thinning
thin = 2
delta = dt*thin


X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]
gradArray = bilinearGradArray(X, covlist)
pars = c(0,0,0,1)

o = optim(pars, lik)
o

n = dim(alldat[[1]])[1]

locs = X
times = alldat[[1]]$t[(0:(n%/%thin -1))*thin +1]
ID = alldat[[1]]$ID[(0:(n%/%thin -1))*thin +1]
fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat




















