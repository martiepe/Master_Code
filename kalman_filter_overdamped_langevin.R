# Load packages
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(optimParallel)
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




###################
## Simulate data ##
###################
#max time for track
Tmax <- 10000
#increment between times
dt <- 0.001
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


#####################
## Thinning Tracks ##
#####################

thin = 50
#divided by six because of 5 extra points
delta = dt*thin/(m+1)


X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]


###########################
## approximating Hessian ##
###########################
#distance between grid points


hessian <- function(z, covlist, par){
  n = floor(z[1])
  m = floor(z[2])
  x = z[1]
  y = z[2]
  
  Delta = 1
  #lower left corner of square being interpolated upon
  
  f = covlist[[1]]$z[(n+100):(n+103), (m+100):(m+103)]*par[1] + 
    covlist[[2]]$z[(n+100):(n+103), (m+100):(m+103)]*par[2] + 
    covlist[[3]]$z[(n+100):(n+103), (m+100):(m+103)]*par[3]
  

  F11 = matrix(c(f[2,2], f[2,3], f[3,2], f[3,3]), nrow = 2, byrow = T)
  
  
  F21 = matrix(c((f[3, 2] - f[1, 2])/(2*Delta),
                 (f[3, 3] - f[1, 3])/(2*Delta),
                 (f[4, 2] - f[2, 2])/(2*Delta),
                 (f[4, 3] - f[2, 3])/(2*Delta)), nrow = 2, byrow = T)
  
  
  F12 = matrix(c((f[2, 3] - f[2, 1])/(2*Delta),
                 (f[2, 4] - f[2, 2])/(2*Delta),
                 (f[3, 3] - f[3, 1])/(2*Delta),
                 (f[3, 4] - f[3, 2])/(2*Delta)), nrow = 2, byrow = T)
  
  
  F22 = matrix(c((f[3,3] - f[3, 1] - f[1, 3] + f[1, 1])/(4*Delta^2), 
                 (f[3,4] - f[3, 2] - f[1, 4] + f[1, 2])/(4*Delta^2),
                 (f[4,3] - f[4, 1] - f[2, 3] + f[2, 1])/(4*Delta^2),
                 (f[4,4] - f[4, 2] - f[2, 4] + f[2, 2])/(4*Delta^2)), nrow = 2, byrow = T)
  
  
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





#checking if hessian is computed correctly

func <- function(s, covlist, par){
  x = s[1]
  y = s[2]

  x1 = floor(x)
  x2 = ceiling(x)
  y1 = floor(y)
  y2 = ceiling(y)
  
  
  if(x%%1 == 0){
    x1 = x-1
    x2 = x+1
    
  }
  if(y%%1 == 0){
    y1 = y-1
    y2 = y+1
  }
  
  cov = par[1]*covlist[[1]]$z + par[2]*covlist[[2]]$z + par[3]*covlist[[3]]$z
  

  f11 = cov[x1+101, y1+101]
  f12 = cov[x1+101, y2+101]
  f21 = cov[x2+101, y1+101]
  f22 = cov[x2+101, y2+101]
  
  

  #covariate value at location
  ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))

}




hessian <- function(x, covlist, par){
  func1 <- function(x){
    return(func(x, covlist, par))
  }
  return(numDeriv::hessian(func1, x))
}
  
  




#################
## likelihood ###
#################

#par: B_1, B_2, B_3, gamma^2
#likelihood with missing observations using kalman filter
delta = dt*thin/(m+1)
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
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      #predicted state estimate
      z_p = F_k %*% z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation
      y = X[i, ] - c(H %*% z_p)
      
      #innovation covariance
      S = H %*% P %*% t(H) + R
      
      #optimal Kalman gain
      K = P %*% t(H) %*% solve(S)
      
      #updated state estimate
      z = z_p + K %*% y
      
      #updated estimate covariance
      P = (diag(1, 2, 2) - K %*% H) %*% P
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - H %*% z_p), mean = c(0,0), sigma = S, log = T)
      
      
      for (k in 1:5) {
        #control vector
        u = bilinearGrad(z, covlist) %*% par[1:3]
      
        #predicted state estimate
        z_p = F_k %*% z + B %*% u
      
        #predicted estimate covariance 
        P = F_k %*% P %*% t(F_k) + Q
        
        #innovation covariance
        S = H %*% P %*% t(H) 
      
        #updated state estimate
        z = z_p 
      
        #updated estimate covariance
        P = P
      }
      
    }
  }
  return(l)
  
}


#likelihood without missing observations using kalman filter
delta = dt*thin
lik1 <- function(par){
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
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      #predicted state estimate
      z_p = F_k %*% z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation
      y = X[i, ] - c(H %*% z_p)
      
      #innovation covariance
      S = H %*% P %*% t(H) + R
      
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
#assumptions:
#R = 0
#F = I
#H = I



#simplified likelihood with missing data using kalman filter
lik2 <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #defining transition covariance matrix
    s = delta*par[4]
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) { 
      #control vector
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      #predicted state estimate
      z_p = z + B %*% u 
      
      #updated state estimate
      z = X[i, ]
      
      #adding likelihood contribution of i-th state
      l = l - dnorm(c(X[i, ] - z_p)[1], 0, (m+1)*s, log = T) - dnorm(c(X[i, ] - z_p)[2], 0, (m+1)*s, log = T)
      
      
      for (k in 1:m) {
        #control vector
        u = bilinearGrad(z, covlist) %*% par[1:3]
        
        #predicted state estimate
        z_p = z + B %*% u
        
        #updated state estimate
        z = z_p 
      }
      
    }
  }
  return(l)
  
}


#likelihood using extended kalman filter
#assuming R = 0
delta = dt*thin/(m+1)
lik3 <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #initial covariance guess
    P = 10*Q
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) {
      #control vector
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      
      F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par)

      #predicted state estimate
      z_p = z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation covariance
      S = P 
      
      #updated state estimate
      z =  X[i, ]
      
      #updated estimate covariance
      P = diag(0,2,2)
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - z_p), mean = c(0,0), sigma = S, log = T)
      
      
      for (k in 1:m) {
        #control vector
        u = bilinearGrad(z, covlist) %*% par[1:3]
        
        
        F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par) 
        
        #predicted state estimate
        z_p = z + B %*% u
        
        #predicted estimate covariance 
        P = F_k %*% P %*% t(F_k) + Q
        
        #innovation covariance
        S = P 
        
        #updated state estimate
        z = z_p 
        
        #updated estimate covariance
        P = P
        
      }
      
    }
  }
  return(l)
  
}



#without missing data points
#i want to use this to find out if this likelihood without missing points generate the proper estimates
delta = dt*thin
lik4 <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:1) {
    #defining transition covariance matrix
    Q = diag(delta*par[4],2,2)
    
    #control matrix
    B = diag(delta*par[4]/2,2,2)
    
    #initial covariance guess
    P = 10*Q
    
    #initial state
    z = X[1, ]
    
    for (i in 2:nrow(X)) {
      #control vector
      u = bilinearGrad(z, covlist) %*% par[1:3]
      
      
      F_k = diag(1,2,2) + (delta*par[4]/2)*hessian(z, covlist, par)
      
      #predicted state estimate
      z_p = z + B %*% u 
      
      #predicted estimate covariance 
      P = F_k %*% P %*% t(F_k) + Q
      
      #innovation covariance
      S = P
      
      #updated state estimate
      z =  X[i, ]
      
      #updated estimate covariance
      P = diag(0,2,2)
      
      #adding likelihood contribution of i-th state
      l = l - dmvnorm(c(X[i, ] - z_p), mean = c(0,0), sigma = S, log = T)
    }
  }
  return(l)
  
}




########## USING PARALELL OPTIM ################
m = 20
delta = dt*thin/(m+1)
t1 = Sys.time()
cl <- makeCluster(12)     # set the number of processor cores
clusterExport(cl, varlist = c("delta", "X", "bilinearGrad", "hessian", "covlist", "dmvnorm", "m", "lik3"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o2 = optimParallel(par=c(0,0,0,1), fn=lik3, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o2

beta

t1 = Sys.time()
o = optim(pars, lik4)
Sys.time() - t1
o



#estimate using full data
X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
gradArray = bilinearGradArray(X, covlist)
locs = X
times = alldat[[1]]$t
ID = alldat[[1]]$ID
fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat


















