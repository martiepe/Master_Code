

# Load packages
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
set.seed(NULL)




#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 2
covlist <- list()
for(i in 1:ncov) {
  covlist[[i]] <- simSpatialCov(lim = lim, nu = 3, rho = 50, sigma2 = 0.1, 
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



UD = getUD(covlist, beta)



ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

UDplot


##############################
## Simulate occurrence data ##
##############################

lambda = exp(par[1]*covlist[[1]]$z + par[2]*covlist[[2]]$z + par[3]*covlist[[3]]$z)



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
s = 0
N
n


#plotting simulations
ggplot() +
  geom_point(aes(Y[,1],Y[,2]), colour = "grey80")+
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))









#####################
## MCMC estimation ##
#####################

#likelihood function
lik <- function(par){
  #log-likelihood
  l = 0
  kappa = nrow(Y)
  #covariate field used to calculate intensity
  cov_field = par[4]*exp(par[1]*covlist[[1]]$z + par[2]*covlist[[2]]$z + par[3]*covlist[[3]]$z)
  
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
      
      
      lambda_sum = lambda_sum + (f11 + f12 + f21 + f22) / 4
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
    
    
    if (x == x2) {
      print(x)
    }
    
    if(x2 == x1){
      x1 = x1-1
      x2 = x2+1
    }
    if(y2 == y1){
      y1 = y1-1
      y2 = y2+1
    }
    
    
    
    f11 = cov_field[x1+lim[2]+1, y1+lim[2]+1]
    f12 = cov_field[x1+lim[2]+1, y2+lim[2]+1]
    f21 = cov_field[x2+lim[2]+1, y1+lim[2]+1]
    f22 = cov_field[x2+lim[2]+1, y2+lim[2]+1]
    
    #intensity at location
    lambda_s = ((y2-y)/(y2-y1))*(f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1)) + ((y-y1)/(y2-y1))*(f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1))
    
    l = l + log(lambda_s)
  }

  return(l)
}



library(parallel)

# Define the function to evaluate
func <- function(x, y) {
  lik(c(x,2,-0.1,y))
}
# Create the grid
x_vals <- seq(-8, 8, length.out = 200)
y_vals <- seq(0.25, 3, length.out = 200)
grid <- expand.grid(x = x_vals, y = y_vals)

# Set up parallel cluster
n_cores <- detectCores(logical = FALSE)
cl <- makeCluster(n_cores)

# Export necessary variables and functions to workers
clusterExport(cl, varlist = c("func", "lik", "covlist", "Y", "x_vals", "y_vals", "grid", "lim"))

# Evaluate in parallel
grid$z <- parLapply(cl, seq_len(nrow(grid)), function(i) {
  func(grid$x[i], grid$y[i])
}) %>% unlist()

# Stop the cluster
stopCluster(cl)

# Plot using ggplot2
ggplot(grid, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(title = "Function Evaluation on Grid",
       x = "beta1", y = "kappa", fill = "func(x, y)") +
  theme_minimal()


#starting values
B_1 = 4
B_2 = 2
B_3 = -0.1
kappa = 1
#matrix containing sampled parameters
params = matrix(c(B_1, B_2, B_3, kappa),ncol = 4)

#number of samples
n = 0
#number of proposed values
r = 0
#random walk variance
s = 0.000001

Sigma = diag(s, 4, 4)

while (n < 100000) {
  
  #log of denominator of acceptance rate
  ld = lik(c(B_1, B_2, B_3, kappa)) + 
    pnorm(B_1, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_2, mean = 0, sd = 10^10,log.p = T) + 
    pnorm(B_3, mean = 0, sd = 10^10,log.p = T) 
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
    a = lik(c(B_1_p, B_2_p, B_3_p, kappa_p)) +
      pnorm(B_1_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_2_p, mean = 0, sd = 10^10,log.p = T) + 
      pnorm(B_3_p, mean = 0, sd = 10^10,log.p = T) +
      pnorm(log(kappa_p), mean = 0, sd = 10^10,log.p = T) -
      ld
    #acceptance probability
    A = min(exp(a), 1)
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
  labs(x = "n", y = "kappa", title = "Samples of kappa") +
  theme_bw()

grid.arrange(p1,p2,p3, p4)



