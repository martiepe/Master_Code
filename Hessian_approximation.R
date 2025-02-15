# Load packages
library(Rhabit)
library(raster)




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


covlist

###########################
## approximating Hessian ##
###########################

hessian <- function(x, y, covlist){
  n = floor(x)
  m = floor(y)
  
  delta = 1
  #lower left corner of square being interpolated upon
  n = 5
  m = 3
  
  C = covlist[[1]]$z*beta[1] + covlist[[2]]$z*beta[2] + covlist[[3]]$z*beta[3]
  
  f = C[(n+100):(n+103), (m+100):(m+103)]
  
  C = 0
  F11 = matrix(c(f[2,2], f[2,3], f[3,2], f[3,3]), nrow = 2)
  
  
  F21 = matrix(c((f[3, 2] - f[1, 2])/(2*delta),
                 (f[3, 3] - f[1, 3])/(2*delta),
                 (f[4, 2] - f[2, 2])/(2*delta),
                 (f[4, 3] - f[2, 3])/(2*delta)), nrow = 2, byrow = T)
  
  
  F12 = matrix(c((f[2, 3] - f[2, 1])/(2*delta),
                 (f[2, 4] - f[2, 2])/(2*delta),
                 (f[3, 3] - f[3, 1])/(2*delta),
                 (f[3, 4] - f[3, 2])/(2*delta)), nrow = 2, byrow = T)
  
  
  F22 = matrix(c((f[3,3] - f[3, 1] - f[1, 3] + f[1, 1])/(4*delta^2), 
                 (f[3,4] - f[3, 2] - f[1, 4] + f[1, 2])/(4*delta^2),
                 (f[4,3] - f[4, 1] - f[2, 3] + f[2, 1])/(4*delta^2),
                 (f[4,4] - f[4, 2] - f[2, 4] + f[2, 2])/(4*delta^2)), nrow = 2, byrow = T)
  
  
  F = cbind(rbind(F11, F21), rbind(F12, F22))
  
  K = matrix(c(1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1), nrow = 4)
  
  A = t(K) %*% F %*% K

  
  #the zero point of the polynomial is in our case the index of the bottom left vertex
  Hxx = c(0, 0, 2, 6*(x-n)) %*% A %*% c(1, (y-m) , (y-m)^2, (y-m)^3)
  
  Hyy = c(1, (x-n), (x-n)^2, (x-n)^3) %*% A %*% c(0, 0, 2, 6*(y-m))
  
  Hxy = c(0, 1, 2*(x-n), 3*(x-n)^2) %*% A %*% c(0, 1, 2*(y-m), 3*(y-m)^2)
  
  H = cbind(rbind(Hxx, Hxy), rbind(Hxy, Hyy))
  
  return(H)
}

