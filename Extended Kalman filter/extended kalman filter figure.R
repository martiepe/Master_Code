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
Tmax <- 500
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 1
#speed parameter for Langevin model
speed <- 5
thin = 10
#divided by six because of 5 extra points
delta = dt*thin/(m+1)
# Time grids
alltimes <- list()
for(i in 1:ntrack)
  alltimes[[i]] <- time




sim_track <- function(){
  # Generate tracks -- This is very computational and may take a few hours.
  alldat <- lapply(alltimes, function(times) {
    simLangevinMM(beta = beta, gamma2 = speed, times = times,
                loc0 = c(0, 0), cov_list = covlist)
  })

  #changed the start position from runif(2, -20, 20) to c(0,0)

  # Add ID column
  for(zoo in 1:ntrack)
    alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])

  #thinning
  X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
  X = X[(0:(nrow(X)%/%thin -1))*thin +1, ]
  
  return(X)
}



###########################
## approximating Hessian ##
###########################

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



#################
## likelihood ###
#################

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

#beta1
results = matrix(data = NA, nrow = 200, ncol = 5)


for (i in 1:50) {
  X = sim_track
  for (j in 1:4) {
    m = c(0,1,5,10)[j]
    delta = dt*thin/(m+1)
    
    cl <- makeCluster(12)     # set the number of processor cores
    clusterExport(cl, varlist = c("delta", "X", "bilinearGrad", "hessian", "covlist", "dmvnorm", "m", "lik3"))
    setDefaultCluster(cl=cl) # set 'cl' as default cluster
    o = optimParallel(par=c(0,0,0,1), fn=lik3, lower=c(-Inf, -Inf, -Inf, 0))
    stopCluster(cl)
    
    par = o$par
    results[4*(i-1)+j, ] = c(par, m)
    
    print(m)
  }
}



df = data.frame(beta1 = results[,1], beta2 = results[,2], beta3 = results[,3], gammasq = results[,4], m = reults[,5])


p1 = ggplot(df, aes(beta1, m)) +
  geom_boxplot() +
  theme_bw() 

p2 = ggplot(df, aes(beta2, m)) +
  geom_boxplot() +
  theme_bw() 

p3 = ggplot(df, aes(beta3, m)) +
  geom_boxplot() +
  theme_bw() 

p4 = ggplot(df, aes(gammasq, m)) +
  geom_boxplot() +
  theme_bw() 

























