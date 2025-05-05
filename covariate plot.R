library(ambient)
library(raster)
library(Rhabit)
library(ggplot2)


set.seed(123)
#limits of the covariates
lim <- c(-1, 1, -1, 1)*150
#resolution
resol <- 1
#number of covariates
ncov <- 2
#list of covariates
covlist <- list()
#grids for the covariate values
xgrid <- seq(lim[1], lim[2], by = resol)
ygrid <- seq(lim[3], lim[4], by = resol)
coords <- as.matrix(expand.grid(xgrid, ygrid))

#simulating perlin covariates
for(i in 1:ncov) {
  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}

#including squared distance to origin as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))



# Compute utilisation distribution
beta <- c(4,2,-0.1)
UD <- getUD(covariates=covlist, beta=beta)

# Plot covariates
c1plot <- plotRaster(rhabitToRaster(covlist[[1]]), scale.name = expression(c[1])) + 
  theme_bw()
c2plot <- plotRaster(rhabitToRaster(covlist[[2]]), scale.name = expression(c[2])) + 
  theme_bw()
UDplot <- plotRaster(rhabitToRaster(UD), scale.name = expression(pi)) + 
  theme_bw()


c1plot
c2plot
UDplot
