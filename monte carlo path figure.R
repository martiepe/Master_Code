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
times = seq(0, 0.1, 0.01)

#simulation
simLangevinMM(beta = beta, gamma2 = speed, times = times, loc0 = c(0, 0), cov_list = covlist)


r = rhabitToRaster(covlist[[1]])
r$val
r[1, 2,]
r[1,1, ]



df = expand.grid(x = covlist[[2]]$x, y = covlist[[2]]$y)




z = 
length(z)


z = covlist[[2]]$z[1:201, 1:201]

d = data.frame(x = df$x, y = df$y, z = as.array(covlist[[2]]$z))

d$

df$

ggplot(df) +
  geom_tile(aes(x = x, y = y, fill = z))


p <- plotRaster(rhabitToRaster(covlist[[2]]))
typeof(p)
p +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 6))


covlist[[1]]





pp <- function (n,r=4) {
  x <- seq(-r*pi, r*pi, len=n)
  df <- expand.grid(x=x, y=x)
  df$r <- sqrt(df$x^2 + df$y^2)
  df$z <- cos(df$r^2)*exp(-df$r/6)
  df
}
pp(20)

expand.grid(x = 1:10, y = 1:10)


















