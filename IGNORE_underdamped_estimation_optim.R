library(Rhabit)
library(mvtnorm)
library(ggplot2)
library(stats)
gamma = 1
sigmasq = 1
delta  = 0.01
x_0 = c(0,0)



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





####################
## Simulate tracks #
####################
ntrack = 10
T_max = 5000
dt = 0.01

gamma = 2




alldat = list()


for (i in 1:ntrack) {
  
  Z = matrix(c(0,0,0,0), ncol = 4)

  x = c(0,0)
  v = c(0,0)

  t = 0
  z = c(x, v)
  
  
  while (t < T_max) {
  grad = bilinearGrad(x, covlist)
  mu_x = x + v*((1-exp(-gamma*dt))/gamma) - (sigmasq/gamma)*(dt - (1-exp(-gamma*dt)))*c(sum(grad[1, ]), sum(grad[2, ]))
  mu_v = v*exp(-gamma*dt) - (sigmasq/gamma)*c(sum(grad[1, ]), sum(grad[2, ]))
  mu = c(mu_x, mu_v)
  
  sigma_11 = sigmasq*(dt + (2/gamma)*exp(-gamma*dt) - exp(-2*gamma*dt)/(2*gamma) - 3/(2*gamma))*diag(1,2,2)
  
  sigma_21 = sigmasq*(1/gamma - (2/gamma)*exp(-gamma*dt) + (1/gamma)*exp(-2*gamma*dt))*diag(1,2,2)
  
  sigma_22 = (2*sigmasq/gamma)*(1 - exp(-2*gamma*dt) )*diag(1,2,2)
  
  sigma = cbind(rbind(sigma_11, sigma_21), rbind(t(sigma_21), sigma_22))
  
  z = rmvnorm(1, mean = mu, sigma = sigma)
  
  Z = rbind(Z,z)
  
  t = t + dt
  
  }
  alldat[[i]] = Z
}

alldat[[1]]



#######################
# defining likelihood #
#######################

###NOTE: R_k and Q_k have been switched?##################

#par = c(beta, gamma, sigmasq)


#defining likelihood function based on Kalman filter
lik <- function(par){
  #log-likelihood
  l = 0
  
  for (j in 1:ntrack) {
    
    thin = 2
    Z_org = alldat[[j]]
    Z = Z_org[(0:(dim(Z_org)[1]%/%thin -1))*thin +1, ]
    
    X = Z[, c(1,2)]
    
    
    
    
    
  #observation matrix
  H = cbind(diag(1,2,2), diag(0,2,2))
  #transition matrix
  F_k = rbind(cbind(diag(1,2,2), ((1-exp(-par[4]*dt))/par[4])*diag(1,2,2)), cbind(diag(0,2,2), exp(-par[4]*dt)*diag(1,2,2)))
  #control matrix
  B = rbind((par[5]/par[4])*(dt - ((1-exp(-par[4]*dt)/par[4])))*diag(1,2,2), - (par[4]/par[4])*(1-exp(-par[4]*dt))*diag(1,2,2))
  #observation covariance
  #R = diag(10^(-6),2,2)
  
  
  #defining transition covariance matrix
  sigma_11 = (2*par[4]/par[4])*(dt + (2/par[4])*exp(-par[4]*dt) - exp(-2*par[4]*dt)/(2*par[4]) - 3/(2*par[4]))*diag(1,2,2)
  sigma_21 = (par[4]/par[4])*(1 - 2*exp(-par[4]*dt) + exp(-2*par[4]*dt))*diag(1,2,2)
  sigma_22 = (2*par[4]/par[4])*(1 - exp(-2*par[4]*dt) )*diag(1,2,2)
  Q = cbind(rbind(sigma_11, sigma_21), rbind(t(sigma_21), sigma_22))
  

  #initial velocity guess
  v_0 = c((X[1,1] - X[2,1])/dt, (X[1,2] - X[2,2])/dt)
  #initial covariance guess
  P = 10*Q
  
  
  #state estimate
  z = c(X[1, ], v_0)
  
  for (i in 2:(dim(X)[1])) {
    
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
    P = (diag(1, 4, 4) - K %*% H) %*% P
    
    #adding likelihood contribution of i-th state
    l = l - dmvnorm(c(X[i, ]), mean = c(H %*% z_p), sigma = S, log = T)
  }
  }
  
  
  
  return(l)
  
}

load("C:\Users\marti\OneDrive\Skrivebord\master\langevin\underdamped_10tracks.RData")


gradArray = bilinearGradArray(Z[, c(1,2)], covlist)


#####using optim#################


dt = dt*thin

#starting values
pars = c(0,0,0,1,1)

optim(pars, lik)







