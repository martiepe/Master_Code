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


####################
# Hessian function #
####################

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


beta = c(4, 2, -0.1)
#beta = 50*beta


#simulation
X = simLangevinMM(beta = beta, gamma2 = speed, times = times, loc0 = c(0, 0), cov_list = covlist)



# Eigen decomposition

c = 4.605
axes_lengths = matrix(nrow = 11, ncol = 2)
angle = c()
path = matrix(nrow = 11, ncol = 2)

#mesh nodes
m = thin -1
delta = dt*thin/(m+1)
Q = diag(delta*5,2,2)
B = diag(delta*5/2,2,2)
P = Q
eig <- eigen(P)
axes_lengths[1, ] = sqrt(eig$values * c)
angle = c(angle, atan2(eig$vectors[2,1], eig$vectors[1,1]))
#predicted state estimate
x = c(0,0)
path[1,] = x


for (k in 1:(m+1)) {
  #control vector
  u = bilinearGrad(x, covlist) %*% beta[1:3]
  print(u)
  
  F_k = diag(1,2,2) + (delta*5/2)*hessian(x, covlist, beta) 
  
  #predicted state estimate
  x = x + B %*% u
  path[k+1,] = x
  #predicted estimate covariance 
  P = F_k %*% P %*% t(F_k) + Q
  eig <- eigen(P)
  axes_lengths[k+1, ] = sqrt(eig$values * c)
  angle = c(angle, atan2(eig$vectors[2,1], eig$vectors[1,1]))  
  #innovation covariance
  S = P 
  
  
}

length(angle)

# Build a data frame with ellipse parameters
ellipse_df <- data.frame(
  x = path[,1],
  y = path[,2],
  a = axes_lengths[,1],
  b = axes_lengths[,2],
  angle = angle,
  group = 1:11
)



ggplot() +
  geom_path(aes(X[,1], X[,2]), color = "red") +
  geom_path(aes(path[,1], path[,2]), color = "black") +
  ggforce::geom_ellipse(aes(x0 = x, y0 = y, a = a, b = b, angle = angle),
                        data = ellipse_df, fill = "blue", alpha = 0) +
  labs(x = "x", y = "y")







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

