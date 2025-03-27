library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)






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
Tmax <- 200
dt <- 0.01
time <- seq(0, Tmax, by = dt)
ntrack <- 1
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

##################testing what length to use ##############

X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
gradArray = bilinearGradArray(X, covlist)
locs = X
times = alldat[[1]]$t
ID = alldat[[1]]$ID
fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat

beta

#plotting trails
{
ggplot()+
  geom_path(aes(alldat[[1]]$x[(0:2500)*200+1], alldat[[1]]$y[(0:2500)*200+1]),color = "grey") +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))


A = data.frame(alldat[[1]])
for (i in 2:100) {
  A = rbind(A, alldat[[i]])
}
length(A$ID)
ggplot()+
  geom_path(aes(A$x[(0:100000)*500+1], A$y[(0:100000)*500+1]), colour = A$ID[(0:100000)*500+1]) +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))
}






##############################
## Simulate occurrence data ##
##############################
kappa = 0.02
lambda = kappa*exp(beta[1]*covlist[[1]]$z + beta[2]*covlist[[2]]$z + beta[3]*covlist[[3]]$z)


#intensity of study area
lambda_max = max(lambda)



#simulating number of points
n = rpois(1, 4*lim[2]*lim[2]*lambda_max)

#simulating points
s = matrix(runif(2*n, min = -100, max = 100), ncol = 2)

#proper number of points
N = 0
#proper locations
Y = c()


#trimming number of points to match inhomogeneous poisson process
for (i in 1:n) {
  x = s[i, 1]
  y = s[i, 2]
  
  
  x1 = floor(s[i, 1])
  x2 = ceiling(s[i, 1])
  y1 = floor(s[i, 2])
  y2 = ceiling(s[i, 2])
  
  f11 = lambda[x1+101, y1+101]
  f12 = lambda[x1+101, y2+101]
  f21 = lambda[x2+101, y1+101]
  f22 = lambda[x2+101, y2+101]
  
  #intensity at location
  lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
  
  #acceptance rate
  p = lambda_s/lambda_max
  
  if(p >= runif(1)){
    N = N + 1
    Y = c(Y, s[i, ])
  }
}

Y = matrix(Y, ncol = 2)


N

#plotting simulations
ggplot() +
  geom_point(aes(S[,1],S[,2]), colour = "grey80")+
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))



############### testing what intensity to use ##################

lik <- function(Y, B_1, B_2, B_3, kappa){
  l = 0
  #adding likelihood from occurrence data
  
  #covariate field used to calculate intensity
  cov_field = kappa*exp(B_1*covlist[[1]]$z + B_2*covlist[[2]]$z + B_3*covlist[[3]]$z)
  
  #finding intensity of study area
  lambda_sum = 0
  for (i in -lim[2]:(lim[2]-1)) {
    for (j in -lim[2]:(lim[2]-1)) {
      x1 = i
      x2 = i+1
      y1 = j
      y2 = j+1
      
      f11 = (cov_field[x1+lim[2]+1, y1+lim[2]+1])
      f12 = (cov_field[x1+lim[2]+1, y2+lim[2]+1])
      f21 = (cov_field[x2+lim[2]+1, y1+lim[2]+1])
      f22 = (cov_field[x2+lim[2]+1, y2+lim[2]+1])
      
      
      lambda_sum = lambda_sum + ((y2+y1)*((x2+x1)*f11 + (x2-3*x1)*f21) + (y2-3*y1)*((x2+x1)*f12 + (x2-3*x1)*f22))/4
    }
  }
  l = l - lambda_sum
  
  #intensity of locations
  for (i in 1:nrow(Y)) {
    x = Y[i, 1]
    y = Y[i, 2]
    
    
    x1 = floor(Y[i, 1])
    x2 = ceiling(Y[i, 1])
    y1 = floor(Y[i, 2])
    y2 = ceiling(Y[i, 2])
    
    f11 = cov_field[x1+lim[2]+1, y1+lim[2]+1]
    f12 = cov_field[x1+lim[2]+1, y2+lim[2]+1]
    f21 = cov_field[x2+lim[2]+1, y1+lim[2]+1]
    f22 = cov_field[x2+lim[2]+1, y2+lim[2]+1]
    
    #intensity at location
    lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
    
    l = l + log(lambda_s)
  }
  return(min(-l, 1e8))
}


t1 = Sys.time()
cl <- makeCluster(10)     # set the number of processor cores
clusterExport(cl, varlist = c("covlist", "lik", "lim", "Y", "getUD"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
o2 = optimParallel(par=c(0,0,0,1), fn=lik, lower=c(-Inf, -Inf, -Inf, .0001))
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o2

lik(Y, 1, B_2, B_3, kappa)

#starting values
B_1 = 0
B_2 = 0
B_3 = 0
kappa = 1
#matrix containing sampled parameters
params = matrix(c(B_1, B_2, B_3, kappa),ncol = 4)

#number of samples
n = 0
#number of proposed values
r = 0
#random walk variance
s = 0.0001

Sigma = diag(s, 4,4)

while (n < 100000) {
  
  #log of denominator of acceptance rate
  ld = lik(Y, B_1, B_2, B_3, kappa) + 
    pnorm(B_1, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_2, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_3, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(log(kappa), mean = 0, sd = 10^10,log.p = T) 
  
  #value signifying if proposed value has been accepted
  d = F
  while(d == F){
    prop = rmvnorm(1, mean = c(B_1, B_2, B_3, log(kappa)), sigma = Sigma)
    
    B_1_p = prop[1]
    B_2_p = prop[2]
    B_3_p = prop[3]
    kappa_p = exp(prop[4])
    r = r+1
    
    #compute log of acceptance rate
    a = lik(Y, B_1_p, B_2_p, B_3_p, kappa) +
      pnorm(B_1_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_2_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_3_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(log(kappa_p), mean = 0, sd = 10^10,log.p = T) -
      ld
    #acceptance probability
    A = min(exp(a), 1)
    #print(A)
    #if the proposal is accepted
    if(runif(1) < A){
      params = rbind(params, c(B_1_p, B_2_p, B_3_p, kappa_p))
      B_1 = B_1_p
      B_2 = B_2_p
      B_3 = B_3_p
      kappa = kappa_p
      #proposed value has been accepted
      d = T
      #sample has increased
      n = n + 1
      
      print(n)
    }
  }
}





pars = data.frame(B_1 = params[, 1], B_2 = params[,2], B_3 = params[,3], kappa = params[,4])


p1 = 0
p2 = 0
p3 = 0
p4 = 0


#plotting parameter samples
p1 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,1])) +
  labs(x = "n", y = "Beta_1", title = "Samples of Beta_1") +
  theme_bw()
p2 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,2])) +
  labs(x = "n", y = "Beta_2", title = "Samples of Beta_2") +
  theme_bw()
p3 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,3])) +
  labs(x = "n", y = "Beta_3", title = "Samples of Beta_3") +
  theme_bw()
p4 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,4])) +
  labs(x = "n", y = "kappa", title = "Samples of kappa") +
  theme_bw()

grid.arrange(p1,p2,p3,p4)



#making histograms of samples
p1 <- ggplot() +
  geom_histogram(aes(params[200:10000,1]), color = "black", fill = "grey") +
  labs(x = "Beta_1", y = "count", title = "historgram of Beta_1 samples") +
  geom_vline(xintercept = 4, linetype="dotted", color = "red") +
  theme_bw()
p2 <- ggplot() +
  geom_histogram(aes(params[200:10000,2]), color = "black", fill = "grey") +
  labs(x = "Beta_2", y = "count", title = "historgram of Beta_2 samples") +
  geom_vline(xintercept = 2, linetype="dotted", color = "red") +
  theme_bw()
p3 <- ggplot() +
  geom_histogram(aes(params[200:10000,3]), color = "black", fill = "grey") +
  labs(x = "Beta_3", y = "count", title = "historgram of Beta_3 samples") +
  geom_vline(xintercept = -0.1, linetype="dotted", color = "red") +
  theme_bw()
p4 <- ggplot() +
  geom_histogram(aes(params[200:10000,4]), color = "black", fill = "grey") +
  labs(x = "gammasq", y = "count", title = "historgram of kappa samples") +
  geom_vline(xintercept = 5, linetype="dotted", color = "red") +
  theme_bw()

grid.arrange(p1,p2,p3,p4)










#####################
## MCMC estimation ##
#####################

#number of tracks used in estimation
ntrack = 10
#thinning factor used
thinning = 2
#step length
Delta = dt*thinning

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
lik <- function(X, Y, B_1, B_2, B_3, gammasq, kappa, ntrack, Delta){
  #log-likelihood
  l = 0
  #adding likelihood from track data
  for (i in 1:ntrack) {
    locs = X[[i]]
    len = dim(locs)[1]
    
    
    l = l + sum(dmvnorm(locs[2:len,] - locs[1:(len-1),] - 0.5*gammasq*Delta*(B_1*grad1[1:(len-1),,i] + B_2*grad2[1:(len-1),,i] + B_3*grad3[1:(len-1),,i]),mean = c(0,0), sigma = gammasq*Delta*diag(1, 2, 2), log = T))
  }
  #adding likelihood from occurrence data
  
  #covariate field used to calculate intensity
  cov_field = exp(B_1*covlist[[1]]$z + B_2*covlist[[2]]$z + B_3*covlist[[3]]$z)
  
  #finding intensity of study area
  lambda_sum = 0
  for (i in -lim[2]:(lim[2]-1)) {
    for (j in -lim[2]:(lim[2]-1)) {
      x1 = i
      x2 = i+1
      y1 = j
      y2 = j+1
      
      f11 = kappa*(cov_field[x1+lim[2]+1, y1+lim[2]+1])
      f12 = kappa*(cov_field[x1+lim[2]+1, y2+lim[2]+1])
      f21 = kappa*(cov_field[x2+lim[2]+1, y1+lim[2]+1])
      f22 = kappa*(cov_field[x2+lim[2]+1, y2+lim[2]+1])
      
      
      lambda_sum = lambda_sum + ((y2+y1)*((x2+x1)*f11 + (x2-3*x1)*f21) + (y2-3*y1)*((x2+x1)*f12 + (x2-3*x1)*f22))/4
    }
  }
  l = l - lambda_sum
  
  #intensity of locations
  for (i in 1:nrow(Y)) {
    x = Y[i, 1]
    y = Y[i, 2]
    
    
    x1 = floor(Y[i, 1])
    x2 = ceiling(Y[i, 1])
    y1 = floor(Y[i, 2])
    y2 = ceiling(Y[i, 2])
    
    f11 = kappa*cov_field[x1+lim[2]+1, y1+lim[2]+1]
    f12 = kappa*cov_field[x1+lim[2]+1, y2+lim[2]+1]
    f21 = kappa*cov_field[x2+lim[2]+1, y1+lim[2]+1]
    f22 = kappa*cov_field[x2+lim[2]+1, y2+lim[2]+1]
    
    #intensity at location
    lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
    
    l = l + log(lambda_s)
  }
  
  
  return(l)
}








#starting values
B_1 = 0
B_2 = 0
B_3 = 0
gammasq = 1
kappa = 1
#matrix containing sampled parameters
params = matrix(c(B_1, B_2, B_3, gammasq, kappa),ncol = 5)

#number of samples
n = 0
#number of proposed values
r = 0
#random walk variance
s = 0.0001

Sigma = diag(s, 5,5)

while (n < 100000) {
  
  #log of denominator of acceptance rate
  ld = lik(X, Y, B_1, B_2, B_3, gammasq, kappa, ntrack, Delta) + 
    pnorm(B_1, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_2, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_3, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(log(gammasq), mean = 0, sd = 10^10,log.p = T) +
    pnorm(log(kappa), mean = 0, sd = 10^10,log.p = T) 

  #value signifying if proposed value has been accepted
  d = F
  while(d == F){
    prop = rmvnorm(1, mean = c(B_1, B_2, B_3, log(gammasq), log(kappa)), sigma = Sigma)
    
    B_1_p = prop[1]
    B_2_p = prop[2]
    B_3_p = prop[3]
    gammasq_p = exp(prop[4])
    kappa_p = exp(prop[5])
    r = r+1
    
    #compute log of acceptance rate
    a = lik(X, Y, B_1_p, B_2_p, B_3_p, gammasq_p, kappa, ntrack, Delta) +
      pnorm(B_1_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_2_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_3_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(log(gammasq_p), mean = 0, sd = 10^10,log.p = T) +
      pnorm(log(kappa_p), mean = 0, sd = 10^10,log.p = T) -
      ld
    
    #acceptance probability
    A = min(exp(a), 1)
    print(A)
    #if the proposal is accepted
    if(runif(1) < A){
      params = rbind(params, c(B_1_p, B_2_p, B_3_p, gammasq_p, kappa_p))
      B_1 = B_1_p
      B_2 = B_2_p
      B_3 = B_3_p
      gammasq = gammasq_p
      kappa = kappa_p
      #proposed value has been accepted
      d = T
      #sample has increased
      n = n + 1
      
      print(n)
    }
  }
}






pars = data.frame(B_1 = params[, 1], B_2 = params[,2], B_3 = params[,3], gammasq = params[,4], kappa = params[,5])

save(pars, file = "MCMC_sample3.Rda")

p1 = 0
p2 = 0
p3 = 0
p4 = 0


#plotting parameter samples
p1 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,1])) +
  labs(x = "n", y = "Beta_1", title = "Samples of Beta_1") +
  theme_bw()
p2 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,2])) +
  labs(x = "n", y = "Beta_2", title = "Samples of Beta_2") +
  theme_bw()
p3 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,3])) +
  labs(x = "n", y = "Beta_3", title = "Samples of Beta_3") +
  theme_bw()
p4 <- ggplot() +
  geom_line(aes(1:dim(params)[1], params[,4])) +
  labs(x = "n", y = "gamma_squared", title = "Samples of gamma_squared") +
  theme_bw()

grid.arrange(p1,p2,p3,p4)



#making histograms of samples
p1 <- ggplot() +
  geom_histogram(aes(params[200:10000,1]), color = "black", fill = "grey") +
  labs(x = "Beta_1", y = "count", title = "historgram of Beta_1 samples") +
  geom_vline(xintercept = 4, linetype="dotted", color = "red") +
  theme_bw()
p2 <- ggplot() +
  geom_histogram(aes(params[200:10000,2]), color = "black", fill = "grey") +
  labs(x = "Beta_2", y = "count", title = "historgram of Beta_2 samples") +
  geom_vline(xintercept = 2, linetype="dotted", color = "red") +
  theme_bw()
p3 <- ggplot() +
  geom_histogram(aes(params[200:10000,3]), color = "black", fill = "grey") +
  labs(x = "Beta_3", y = "count", title = "historgram of Beta_3 samples") +
  geom_vline(xintercept = -0.1, linetype="dotted", color = "red") +
  theme_bw()
p4 <- ggplot() +
  geom_histogram(aes(params[200:10000,4]), color = "black", fill = "grey") +
  labs(x = "gammasq", y = "count", title = "historgram of gammasq samples") +
  geom_vline(xintercept = 5, linetype="dotted", color = "red") +
  theme_bw()

grid.arrange(p1,p2,p3,p4)


#calculating acceptence rate
n/r

#finding mean of sampled parameters
B_1 = mean(params[900:length(params[,1]),1])
B_2 = mean(params[900:length(params[,1]),2])
B_3 = mean(params[900:length(params[,1]),3])
gammasq = mean(params[900:length(params[,1]),4])

B_1
B_2
B_3
gammasq

#finding covariance of samples parameters
pars = data.frame(B_1 = params[900:length(params[,1]),1], B_2 = params[900:length(params[,1]),2], B_3 = params[900:length(params[,1]),3], gammasq = params[900:length(params[,1]),4])

Sigma = cov(pars)
Sigma


#making autocorrelation plot
auto1 = acf(params[,1], lag.max = 15, pl = F)
auto2 = acf(params[,2], lag.max = 15, pl = F)
auto3 = acf(params[,3], lag.max = 15, pl = F)
auto4 = acf(params[,4], lag.max = 15, pl = F)

p1 <- ggplot() +
  geom_col(aes(0:15, auto1$acf), color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "lag", y = "correlation", title = "beta_1")
p2 <- ggplot() +
  geom_col(aes(0:15, auto2$acf), color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "lag", y = "correlation", title = "beta_2")
p3 <- ggplot() +
  geom_col(aes(0:15, auto3$acf), color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "lag", y = "correlation", title = "beta_3")
p4 <- ggplot() +
  geom_col(aes(0:15, auto4$acf), color = "black", fill = "grey") +
  theme_bw() +
  labs(x = "lag", y = "correlation", title = "gamma squared")


grid.arrange(p1,p2,p3,p4)

#calculating ESS
library(coda)
effectiveSize(params[,1])








