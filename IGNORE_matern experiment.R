library(Rhabit)


library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)

lim <- c(-1, 1, -1, 1)*100
resol <- 1
covlist = list()



#checking effect of sigmasq
sigma = c(0.1, 1, 10 ,100)

for(i in 1:4) {
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = sigma[i], 
                                resol = resol, raster_like = TRUE)
}

p1 = plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1]))
p2 = plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[1]))
p3 = plotRaster(rhabitToRaster(covlist[[3]]), scale.name = expression(c[1]))
p4 = plotRaster(rhabitToRaster(covlist[[4]]), scale.name = expression(c[1]))

grid.arrange(p1,p2,p3,p4)




#checking effect of phi
phi = c(10, 50, 100, 200)

for(i in 1:4) {
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 0.6, rho = phi[i], sigma2 = 0.1, 
                                resol = resol, raster_like = TRUE)
}

p1 = plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1]))
p2 = plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[1]))
p3 = plotRaster(rhabitToRaster(covlist[[3]]), scale.name = expression(c[1]))
p4 = plotRaster(rhabitToRaster(covlist[[4]]), scale.name = expression(c[1]))

grid.arrange(p1,p2,p3,p4)




#checking effect of kappa
kappa = c(0.1, 1, 10 ,100)

for(i in 1:4) {
  covlist[[i]] <- simSpatialCov(lim = lim, nu = kappa[i], rho = 50, sigma2 = 0.1, 
                                resol = resol, raster_like = TRUE)
}

p1 = plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1]))
p2 = plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[1]))
p3 = plotRaster(rhabitToRaster(covlist[[3]]), scale.name = expression(c[1]))
p4 = plotRaster(rhabitToRaster(covlist[[4]]), scale.name = expression(c[1]))

grid.arrange(p1,p2,p3,p4)

##########################################
#checking how parameters affect estimates#
##########################################


simulateLangevin <- function(x_0, cov_list, gammasq, delta, n_obs){
  X = matrix(ncol = 2,nrow = n_obs)
  X[1,] = x_0
  for (i in 2:n_obs) {
    print(X[i-1,])
    X[i,] = X[i-1,] + gammasq*delta*bilinearGrad(as.vector(X[i-1, ]), rhabitToRaster(cov_list))/2 + rnorm(1,mean = c(0,0), gammasq*delta*diag(1,2,2))
  }
  
  return(X)
}

test = simulateLangevin(c(0,0), covlist, 5, 0.01, 1000)
covlist
?vector
?bilinearGrad
#max time for track
Tmax <- 5000
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 50
#speed parameter for Langevin model
speed <- 5




alltimes <- list()
for(i in 1:ntrack)
  alltimes[[i]] <- time



#generate covariates
covlist = list()

allfits2 <- list()
beta = c(4,-0.2)
thin = 5
kappa = c(0.1, 1, 10 ,100)
phi = c(10, 50, 100, 200)
lim <- c(-1, 1, -1, 1)*300


for (i in 1:4) {
  for (j in 1:4) {
    covlist[[1]] <- simSpatialCov(lim = lim, nu = kappa[j], rho = phi[i], sigma2 = 0.1, 
                                  resol = resol, raster_like = TRUE)
    # Include cubed distance to centre of map as covariate
    xgrid <- seq(lim[1], lim[2], by=resol)
    ygrid <- seq(lim[3], lim[4], by=resol)
    xygrid <- expand.grid(xgrid,ygrid)
    dist2 <- ((xygrid[,1])^3+(xygrid[,2])^3)/100
    covlist[[2]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))
    
    
    
    #
    
    
    cl <- makeCluster(getOption("cl.cores",3))
    clusterExport(cl, c("simLangevinMM", "covlist", "time"), envir = environment())
    
    
    
    #generate tracks
    alldat <- parLapply(cl, alltimes, function(time) {
      simLangevinMM(c(4,-0.2), 5, time, c(0, 0), covlist)
    })
    
    stopCluster(cl)
    
    #alldat <- lapply(alltimes, function(time) {
    #  simLangevinMM(c(4,-0.2), 5, time, c(0, 0), covlist)
    #})
    
    
    
    
    
    
    # keep 50000 locations from each track
    subdat <- lapply(alldat, function(dat) as.matrix(dat[1:50000,]))
    # derive gradient for all locations
    subgrad <- lapply(subdat, function(d)
      bilinearGradArray(locs=d[,c("x", "y")], cov_list=covlist))
    
    #thin out tracks and gradients
    thindat <- lapply(subdat, function(dat) dat[seq(1,nrow(dat),by=thin),])
    thingrad <- lapply(subgrad, function(grad) grad[seq(1,nrow(grad),by=thin),,])
    
    
    # Fit model to the 100 thinned tracks
    allfits2[[4*(i-1)+j]] <- t(mapply(function(dat, grad) {
      ID <- dat[,1]
      locs <- as.matrix(dat[,2:3])
      times <- dat[,4]
      fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=grad)
      return(c(fit$betaHat, fit$gamma2Hat))
    }, dat = thindat, grad = thingrad))
    
    print(3*(i-1)+j)
    
    
  }
}



?simLangevinMM
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[1])) + ggtheme
UD <- getUD(covariates=covlist, beta=beta)
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

c1plot
UDplot

detectCores()

#generate tracks








































