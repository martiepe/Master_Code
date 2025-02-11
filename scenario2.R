
# This script reproduces the simulation study of Section 4.2 of the manuscript
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
library(reshape2)
library(gridExtra)
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

(xygrid[,1])

# Compute utilisation distribution
beta <- c(4,2,-0.1)
UD <- getUD(covariates=covlist, beta=beta)

?getUD

# Plot covariates
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

pdf("sim2cov.pdf", width=12, height=3)
grid.arrange(c1plot, c2plot, UDplot, nrow=1)
dev.off()



c1plot
c2plot
UDplot

###################
## Simulate data ##
###################
#max time for track
Tmax <- 5000
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 100
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

##############
## Estimate ## sparser sampling have longer duration for tracks
##############
# Thinning factors
allthin <- c(1,2,5,10,20,50,100)

#find estimates of coefficients for different times
allfits1 <- lapply(as.list(allthin), function(thin) {
    fit <- t(sapply(alldat, function(dat) {
        # Thin and keep 5000 observations per track
        ind <- thin*(1:5000)
        subdat <- as.matrix(dat[ind,])
        
        # Unpack data set
        ID <- subdat[,1]
        locs <- subdat[,2:3]
        times <- subdat[,4]
        
        # Derive covariate gradients at observed locations
        gradarray <- bilinearGradArray(locs = locs, cov_list = covlist)
        
        # Fit Langevin model
        fit <- langevinUD(locs=locs, times=times, ID = ID, grad_array=gradarray)
        return(c(fit$betaHat, fit$gamma2Hat))
    }))
    return(fit)
})

#############################################
## Estimation 2: constant overall duration ## tracks have the same time-lenght, but some are sparser
#############################################
# Thinning factors
allthin <- c(1,2,5,10,20,50,100)

# keep 50000 locations from each track
subdat <- lapply(alldat, function(dat) as.matrix(dat[1:50000,]))
# derive gradient for all locations
subgrad <- lapply(subdat, function(d)
  bilinearGradArray(locs=d[,c("x", "y")], cov_list=covlist))

allfits2 <- list()
for(i in 1:length(allthin)) {
    cat("Iteration",i,"\n")
    
    # Thin data and array of gradients
    thin <- allthin[i]
    thindat <- lapply(subdat, function(dat) dat[seq(1,nrow(dat),by=thin),])
    thingrad <- lapply(subgrad, function(grad) grad[seq(1,nrow(grad),by=thin),,])
    
    # Fit model to the 100 thinned tracks
    allfits2[[i]] <- t(mapply(function(dat, grad) {
        ID <- dat[,1]
        locs <- as.matrix(dat[,2:3])
        times <- dat[,4]
        fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=grad)
        return(c(fit$betaHat, fit$gamma2Hat))
    }, dat = thindat, grad = thingrad))
}

####################
## Plot estimates ##
####################
sameDuration <- T # TRUE to plot results of analysis with constant duration
if(sameDuration) {
    allfits <- allfits2
} else {
    allfits <- allfits1
}

# Prepare data set of estimates
parnames <- c("beta[1]","beta[2]","beta[3]","gamma^2")
truepar <- data.frame(par = parnames, value = c(beta, speed))
zeroline <- data.frame(yint = 0, par=parnames[1:3])

tmp <- list()
for(i in 1:length(allfits)) {
    df <- as.data.frame(allfits[[i]])
    colnames(df) <- parnames
    df2 <- melt(df)
    colnames(df2) <- c("par", "value")
    tmp[[i]] <- cbind(dt = allthin[i]*dt, df2)
}
plotdat <- do.call(rbind.data.frame, tmp)
plotdat$dt <- as.factor(plotdat$dt)

p <- ggplot(plotdat, aes(dt, value)) + geom_hline(data=zeroline, aes(yintercept=yint), lty=2) + 
    geom_hline(data=truepar, mapping=aes(yintercept=value), lty=2, col=2) + 
    facet_wrap("par", scale = "free", labeller = label_parsed) + geom_boxplot(outlier.size = 0.6) +
    theme(strip.text = element_text(size = 12), axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) + xlab("Time interval") + ylab("Estimated value")

ggsave("sim3par1.pdf", p, width = 8, height = 6, units = "in",
       family = "LM Roman 10", device = cairo_pdf)

p
##############################
## Irregular time intervals ##
##############################
# keep 50000 locations from each track
subdat <- lapply(alldat, function(dat) as.matrix(dat[1:50000,]))
# derive gradient for all locations
subgrad <- lapply(subdat, function(d)
    bilinearGradArray(locs=d[,c("x", "y")], cov_list=covlist))

# Thinning factors
allthin <- c(5, 20)

allfits3 <- list()
for(i in 1:length(allthin)) {
    cat("Iteration",i,"\n")
    
    # Thin data and array of gradients
    thin <- allthin[i]
    # Random indices of irregular observations to keep
    whichkeep <- sort(sample(1:nrow(subdat[[1]]), 
                             size = nrow(subdat[[1]])/thin, 
                             replace = FALSE))
    thindat <- lapply(subdat, function(dat) dat[whichkeep,])
    thingrad <- lapply(subgrad, function(grad) grad[whichkeep,,])
    
    # Fit model to the 100 thinned tracks
    allfits3[[i]] <- t(mapply(function(dat, grad) {
        ID <- dat[,1]
        locs <- as.matrix(dat[,2:3])
        times <- dat[,4]
        fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=grad)
        return(c(fit$betaHat, fit$gamma2Hat))
    }, dat = thindat, grad = thingrad))
}

output <- function(x)
    c(mean(x), quantile(x, prob = c(0.025, 0.975)))

# Compare results from regular and irregular data - for thin = 5
reg5 <- t(apply(allfits2[[3]], 2, output))
irreg5 <- t(apply(allfits3[[1]], 2, output)) 

# Compare results from regular and irregular data - for thin = 20
reg20 <- t(apply(allfits2[[5]], 2, output))
irreg20 <- t(apply(allfits3[[2]], 2, output))
