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
library(ggplot2)
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
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 1, rho = 50, sigma2 = 0.1, 
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
Tmax <- 0.1
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 1
#speed parameter for Langevin model
speed <- 5

# Time grids
times = seq(0, 0.1, dt)

#simulation
X = simLangevinMM(beta = beta, gamma2 = speed, times = times, loc0 = c(0, 0), cov_list = covlist)


ggplot() +
  geom_path(aes(X[,1], X[,2]), color = "red") +
  labs(x = "x", y = "y")



#generate 10 tracks
tracks = matrix(data = NA, nrow = 11*10, ncol = 3)
for (i in 1:10) {
  x = simLangevinMM(beta = beta, gamma2 = speed, times = seq(0, 0.09, 0.01), loc0 = c(0, 0), cov_list = covlist)
  tracks[(11*(i-1)+1):(i*11), 1] = seq(0, 0.1, 0.01)
  tracks[(11*(i-1)+1):(i*11), 2] = c(x[,2], X[11, 2])
  tracks[(11*(i-1)+1):(i*11), 3] = i
}
p1 <- ggplot() +
  geom_path(aes(x=tracks[, 1], y = tracks[, 2], group = tracks[, 3])) +
  geom_path(aes(seq(0, 0.1, 0.01), X[,2]), color = "red") +
  labs(x = "t", y = "x", title = "Langevin Process")



#generating brownian bridges
tracks = matrix(data = NA, nrow = 11*10, ncol = 3)

#generating nodes
sigma = matrix(nrow = 9, ncol = 9)
N = 9
delta = 0.1
mu_x = c()
mu_y = c()
#sigma is the same as long as delta is the same, so it can be computed once
for (k in 1:9) {
  for (m in 1:k) {
    sigma[k,m] = delta*(1 - k/(N+1))*(m/(N+1))
    sigma[m,k] = delta*(1 - k/(N+1))*(m/(N+1))
  }
  mu_x = c(mu_x, X[1, 1] + k*(X[11, 1] - X[1, 1])/(N+1))
  mu_y = c(mu_y, X[1, 2] + k*(X[11, 2] - X[1, 2])/(N+1))
}

for (i in 1:10) {
  x = rmvnorm(1, mean = mu_x, sigma = 5*sigma)
  y = rmvnorm(1, mean = mu_y, sigma = 5*sigma)
  tracks[(11*(i-1)+1):(i*11), 1] = seq(0, 0.1, 0.01)
  tracks[(11*(i-1)+1):(i*11), 2] = c(0, y, X[11, 2])
  tracks[(11*(i-1)+1):(i*11), 3] = i
}

p2 <- ggplot() +
  geom_path(aes(x=tracks[, 1], y = tracks[, 2], group = tracks[, 3])) +
  geom_path(aes(seq(0, 0.1, 0.01), X[,2]), color = "red") +
  labs(x = "t", y = "x", title = "Brownian Bridge")

grid.arrange(p1, p2, nrow = 1)






