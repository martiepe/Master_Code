library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(mvtnorm)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(ambient)
library(mvnfast)
library(optimParallel)





#perlin covariates
lim <- c(-1, 1, -1, 1)*150
resol <- 1
ncov <- 2
covlist <- list()
xgrid <- seq(lim[1], lim[2], by = resol)
ygrid <- seq(lim[3], lim[4], by = resol)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}
# Include squared distance to centre of map as covariate
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



#speed parameter for Langevin model
speed <- 5
set.seed(123)

params = matrix(NA, ncol = 5, nrow = 5*100)

ik = 1
while (ik <= 100) {
  beta <- c(4,2,-0.1)
  dt = 0.01
  n_obs = 5000
  Tmax = n_obs*100*dt
    
  time <- seq(0, Tmax, by = dt)

  # Generate tracks 
  track = sim_langevin(beta, speed, dt, Tmax, c(0,0), covlist)
    
  
  for (jk in 1:5) {
    thin = c(1, 5, 10, 50, 100)[jk]
    delta = dt*thin
    
    
    #thinning tracks
    X = track
    n = nrow(X)
    X = X[(0:(n%/%thin -1))*thin +1, ]
    times = time[(0:(n%/%thin -1))*thin +1]

    
    #shotening tracks to 5000 observations
    X = X[1:5000, ]
    times = times[1:5000]
    ID = rep(1, 5000)
    
    #computing gradient of track observations
    gradArray = bilinearGradVec(X, covlist)
    gradArray <- aperm(gradArray, perm = c(2, 3, 1))
    
    
    fit = langevinUD(locs=X, times=times, ID=ID, grad_array=gradArray)
    
    params[ik*5+jk-5, 1:3] = fit$betaHat
    params[ik*5+jk-5, 4] = fit$gamma2Hat
    params[ik*5+jk-5, 5] = thin
    print(c(fit$betaHat, fit$gamma2Hat))
  }
  
  print(ik)
  
  ik = ik + 1
}




df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], dt = as.factor(dt*params[,5]))
save(df,file="varying_thin_estimates_EM.Rda")

## plotting estimates ##
p1 <- ggplot(data = df, aes(x = dt, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = "Beta_1") +
  theme_bw()

p2 <- ggplot(data = df, aes(x = dt, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = "Beta_2") +
  theme_bw()

p3 <- ggplot(data = df, aes(x = dt, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = "Beta_3") +
  theme_bw()

p4 <- ggplot(data = df, aes(x = dt, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = "gamma^2") +
  theme_bw()


grid.arrange(p1,p2,p3,p4)

