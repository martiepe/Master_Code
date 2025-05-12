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
beta = c(4,2,-0.1)
speed = 5
dt = 0.01
Tmax = 100

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






bilinearGradVec <- function(loc_mat, cov_list) {
  x_grid <- cov_list[[1]]$x
  y_grid <- cov_list[[1]]$y
  n_cov <- length(cov_list)
  n_obs <- nrow(loc_mat)
  
  ix <- findInterval(loc_mat[,1], x_grid)
  iy <- findInterval(loc_mat[,2], y_grid)
  
  valid <- ix > 0 & ix < length(x_grid) & iy > 0 & iy < length(y_grid)
  
  x1 <- x_grid[ix]
  x2 <- x_grid[ix + 1]
  y1 <- y_grid[iy]
  y2 <- y_grid[iy + 1]
  
  dx <- x2 - x1
  dy <- y2 - y1
  lx <- loc_mat[,1]
  ly <- loc_mat[,2]
  
  grad_array <- array(NA_real_, dim = c(n_cov, n_obs, 2))
  
  for (j in seq_len(n_cov)) {
    f11 <- cov_list[[j]]$z[cbind(ix,     iy)]
    f21 <- cov_list[[j]]$z[cbind(ix + 1, iy)]
    f12 <- cov_list[[j]]$z[cbind(ix,     iy + 1)]
    f22 <- cov_list[[j]]$z[cbind(ix + 1, iy + 1)]
    
    dfdx <- ((y2 - ly) * (f21 - f11) + (ly - y1) * (f22 - f12)) / (dy * dx)
    dfdy <- ((x2 - lx) * (f12 - f11) + (lx - x1) * (f22 - f21)) / (dy * dx)
    
    grad_array[j, valid, 1] <- dfdx[valid]
    grad_array[j, valid, 2] <- dfdy[valid]
  }
  
  grad_array
}


sim_langevin <- function(beta, gammasq, dt, Tmax, loc_0, covlist){
  sigma = sqrt(gammasq*dt)
  n_obs = Tmax/dt
  locs = rbind(matrix(loc_0, nrow=1), mvnfast::rmvn(n_obs-1, c(0,0), sigma = diag(sigma, 2,2), isChol = TRUE)) 
  for (i in 1:(n_obs-1)) {
    locs[i+1, ] = locs[i+1, ] + locs[i, ] + (gammasq*dt/2)*bilinearGrad(locs[i, ], covlist)%*%beta 
  }
  
  return(locs)
}


X = sim_langevin(beta, speed, dt, Tmax, c(0,0), covlist)


UD = getUD(covlist, beta)

z = UD$z

coord_vec <- seq(-150, 150, by = 1)   # length = 301

# ── Long-form data frame for the raster ────────────────────────────
df <- data.frame(
  y     = rep(coord_vec, times = length(coord_vec)),   # row coordinate
  x     = rep(coord_vec, each = length(coord_vec)),    # column coordinate
  value = as.vector(z)                               # cell value
)


path_df <- data.frame(
  x = X[, 1],
  y = X[, 2]
)



ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_path(
    data = path_df,
    aes(x = x, y = y),
    colour = "red",
    linewidth = 1,
    inherit.aes = FALSE   # <‑‑ nothing from the global aes is inherited
  ) +
  scale_y_reverse() +
  coord_equal(xlim  = c(-40, 20), ylim  = c(10, -40), expand = FALSE) +
  theme_minimal()






