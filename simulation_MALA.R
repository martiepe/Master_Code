
# This script reproduces the simulation study of Appendix D of the manuscript
# "The Langevin diffusion as a continuous-time model of animal movement and 
# habitat selection". It uses the package Rhabit, which can be downloaded from
# Github as described below.

# Note that the function "mclapply" from the package parallel does not work on 
# Windows, and would need to be replaced by another parallelized routine.

library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(parallel)
set.seed(1)

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
UDrast <- rasterFromXYZ(rasterToGGplot(UD))

###################
## Simulate data ##
###################
nobs <- 1e3
alldt <- c(5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)
ntrack <- 10
speed <- 5

rates <- rep(NA, length(alldt))
t0 <- Sys.time()
for(iter in 1:length(alldt)) {
    cat("Iteration",iter,"-- dt =",alldt[iter],"\n")
    
    dt <- alldt[iter]
    time <- (1:nobs)*dt
    alltimes <- list()
    for(i in 1:ntrack)
        alltimes[[i]] <- time
    
    # Generate tracks
    alldat <- lapply(alltimes, function(time) {
        Rhabit:::simMALA(beta = beta, gamma2 = speed, times = time, 
                loc0 = runif(2, -50, 50), cov_list = covlist)
    })
    
    rates[iter] <- mean(sapply(alldat, function(dat) dat$acc/(dat$acc+dat$rej)))
    cat("Acceptance rate:",rates[iter],"\n")
    print(Sys.time() - t0)
}

###########################
## Plot acceptance rates ##
###########################
p <- qplot(x=alldt, y=rates, xlab="Time interval", ylab="Acceptance rate") + 
  geom_point() + geom_line() + scale_x_continuous(trans = "log", breaks=alldt) + 
  theme(axis.title = element_text(size=14), axis.text = element_text(size=13))

pdf("MALArates.pdf", width=6, height=4)
plot(p)
dev.off()
p
