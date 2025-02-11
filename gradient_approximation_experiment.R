library(Rhabit)

# Load packages
library(raster)
library(ggplot2)
library(viridis)
library(gridExtra)
library(reshape2)

set.seed(123)




#######################
## Define covariates ##
#######################
# Generate two random covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 2
covlist <- list()


# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist1 <- sin(xygrid[,1]*xygrid[,2]/100)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100

covlist[[1]] <- list(x=xgrid, y=ygrid, z=matrix(dist1, length(xgrid), length(ygrid)))
covlist[[2]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))



# Compute utilisation distribution
beta <- c(4,-0.1)
UD <- getUD(covariates=covlist, beta=beta)


# Plot covariates
ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                 legend.title = element_text(size=12), legend.text = element_text(size=12))
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + ggtheme
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + ggtheme
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + ggtheme

#pdf("sim2cov.pdf", width=12, height=3)
grid.arrange(c1plot, c2plot, UDplot, nrow=1)
dev.off()



c1plot
c2plot
UDplot








