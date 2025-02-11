# Load packages
library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(gridExtra)
library(parallel)
library(reshape2)

set.seed(123)


#defining the two covariates
lim <- c(-1, 1, -1, 1)*150
resol <- 1
ncov <- 5
covlist <- list()


covlist[[1]] <- simSpatialCov(lim = lim, nu = 10, rho = 50, sigma2 = 0.1, 
                                resol = resol, raster_like = TRUE)

covlist[[2]] <- simSpatialCov(lim = lim, nu = 10, rho = 10, sigma2 = 0.1, 
                              resol = resol, raster_like = TRUE)

covlist[[3]] <- simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = 0.1, 
                              resol = resol, raster_like = TRUE)

covlist[[4]] <- simSpatialCov(lim = lim, nu = 0.6, rho = 10, sigma2 = 0.1, 
                              resol = resol, raster_like = TRUE)


# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
covlist[[5]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))









#plotting covariates
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
c3plot <- plotRaster(rhabitToRaster(covlist[[3]]), scale.name = expression(c[3])) + ggtheme
c4plot <- plotRaster(rhabitToRaster(covlist[[4]]), scale.name = expression(c[4])) + ggtheme


grid.arrange(c1plot, c2plot, c3plot, c4plot)



#compute UD
UD <- getUD(covariates=covlist, beta = c(1,1,1,1,-0.1))




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
  simLangevinMM(beta = c(1,1,1,1, -0.1), gamma2 = speed, times = times,
                loc0 = c(0, 0), cov_list = covlist)
})



# Add ID column
for(zoo in 1:ntrack)
  alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])









#############################################
## Estimation 2: constant overall duration ## tracks have the same time-lenght, but some are sparser
#############################################

beta = c(1,1,1,1,-0.1)

# Thinning factors
allthin <- c(1,2,5,10,20,50,100)

# keep 50000 locations from each track
subdat <- lapply(alldat, function(dat) as.matrix(dat[1:50000,]))
# derive gradient for all locations
subgrad <- lapply(subdat, function(d)
  bilinearGradArray(locs=d[,c("x", "y")], cov_list=covlist))

allfits <- list()
for(i in 1:length(allthin)) {
  cat("Iteration",i,"\n")
  
  # Thin data and array of gradients
  thin <- allthin[i]
  thindat <- lapply(subdat, function(dat) dat[seq(1,nrow(dat),by=thin),])
  thingrad <- lapply(subgrad, function(grad) grad[seq(1,nrow(grad),by=thin),,])
  
  # Fit model to the 100 thinned tracks
  allfits[[i]] <- t(mapply(function(dat, grad) {
    ID <- dat[,1]
    locs <- as.matrix(dat[,2:3])
    times <- dat[,4]
    fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=grad)
    return(c(fit$betaHat, fit$gamma2Hat))
  }, dat = thindat, grad = thingrad))
}


allfits[[1]][,1]

cov1 = data.frame(allfits = c(allfits[[1]][1,], allfits[[2]][1,], allfits[[3]][1,], allfits[[4]][1,], allfits[[5]][1,], allfits[[6]][1,], allfits[[6]][1,]),
                  thin = c(rep(1,100), rep(2,100), rep(3,100), rep(4,100), rep(5,100), rep(6,100), rep(7,100)))

ggplot( data = cov1) +
  geom_boxplot(aes(allfits, thin))




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
parnames <- c("beta[1]","beta[2]","beta[3]", "beta[4]", "beta[5]","gamma^2")
truepar <- data.frame(par = parnames, value = c(beta, speed))
zeroline <- data.frame(yint = 0, par=parnames[1:5])

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










