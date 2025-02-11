library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
library(reshape2)
library(gridExtra)
library(mvtnorm)
set.seed(1)






library(INLA)
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
Tmax <- 5000
dt <- 0.01
time <- seq(0, Tmax, by = dt)
ntrack <- 10
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



#plotting trails
length(alldat[[1]]$x)

ggplot()+
  geom_path(aes(alldat[[1]]$x[(0:2500)*200+1], alldat[[1]]$y[(0:2500)*200+1]),color = "grey") +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))


A = data.frame(alldat[[1]])



tail(A)

for (i in 2:100) {
  A = rbind(A, alldat[[i]])
}
length(A$ID)
ggplot()+
  geom_path(aes(A$x[(0:100000)*500+1], A$y[(0:100000)*500+1]), colour = A$ID[(0:100000)*500+1]) +
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))








##############################
## Simulate occurrence data ##
##############################

#getting probabilities from utilization distribution
P = UD$z/sum(UD$z)

#making a grid encoded as positive integers
S = 0:(201*201-1)

#making list of probabilities
p = c()
for (i in P) {
  p = c(p,i)
}

#simulating occurrences
s = sample(S, 50, replace = T, prob = p)

#decoding simulations into Cartesian coordinates
x = s%%201 -100
y = s%/%201 -100

#adding position in each block
u_x = runif(50) - 0.5
u_y = runif(50) - 0.5

x = x + u_x
y = y + u_y

#plotting simulations
ggplot() +
  geom_point(aes(x,y))+
  coord_cartesian(xlim = c(-100, 100), ylim= c(-100, 100))

Y = matrix(c(x,y), ncol = 2)


######################################
# constructing data frame for fitting#
######################################

covlist[[1]]
alldat[[1]]$x[2]


Y = c(alldat[[1]]$x, alldat[[1]]$y)
length(Y)

Z= matrix(ncol = ncov, nrow = length(Y))
ncov = 3
Delta = 0.01


for (i in 1:(length(Y)/2)) {
  
    loc = c(alldat[[j]]$x[i], alldat[[j]]$y[i])
    Z[i,] = bilinearGrad(loc, covlist)[1, ]*Delta/2
    Z[length(Y)/2 + i,] = bilinearGrad(loc, covlist)[2, ]*Delta/2
    
  
}


Z
dim(Z)



###################
#fitting the model#
###################

formula <- d ~ c_1 + c_2 + c_3
imod <- inla(formula, family="gaussian", data=chredlin)

?inla
?control.family
?custom




library(INLA)

# Simulated data
n <- 500001
#covariate

# Data frame for INLA
data <- data.frame(y = Y, c_1 = Z[,1], c_2 = Z[,2], c_3 = Z[,3])

# Formula: model linear predictor
formula <- y ~ c_1 + c_2 + c_3



# Use inla()
result <- inla(
  formula,
  family = "gaussian", # Define Gaussian family
  data = data,
  control.family = list(hyper = list(prec = list(prior = "loggamma")))
)

# Summary of results
summary(result)


l <- lm(formula, data)


((mean(l$residuals^2)*length(Y)^2))/(length(Y)-3)
?control.lp.scale


summary(l)

norm(l$residuals)/(length(Y)-3)


l$residuals

mean((data$y-predict(l))^2)


predict.lm(l)$fit


