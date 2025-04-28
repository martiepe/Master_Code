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
library(mvtnorm)
library(ambient)
library(mvnfast)
library(optimParallel)





#perlin covariates
lim <- c(-1, 1, -1, 1)*100
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



#number of tracks to be generated
ntrack <- 1
#speed parameter for Langevin model
speed <- 5



set.seed(123)



params = matrix(NA, ncol = 5, nrow = 3*100)

ik = 1
while (ik <= 100) {
  for (jk in 1:3) {
    beta <- c(4,2,-0.1)
    thin = c(10, 50, 100)[jk]
    dt = 0.01
    delta = dt*thin
    N = thin-1
    M = 50
    n_obs = 5000
    Tmax = n_obs*thin*dt
    
    
    
    time <- seq(0, Tmax, by = dt)
    alltimes <- list()
    for(i in 1:ntrack)
      alltimes[[i]] <- time
    
    # Generate tracks 
    alldat <- lapply(alltimes, function(times) {
      simLangevinMM(beta = beta, gamma2 = speed, times = times,
                    loc0 = c(0, 0), cov_list = covlist, silent = TRUE)
    })
    
    # Add ID column
    for(zoo in 1:ntrack)
      alldat[[zoo]] <- cbind(ID = rep(zoo, length(time)), alldat[[zoo]])
    
    
    #thinning tracks
    X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
    n = nrow(X)
    X = X[(0:(n%/%thin -1))*thin +1, ]
    
    
    gradArray = bilinearGradArray(X, covlist)
    
    locs = X
    times = alldat[[1]]$t[(0:(n%/%thin -1))*thin +1]
    ID = alldat[[1]]$ID[(0:(n%/%thin -1))*thin +1]
    
    gradArray = bilinearGradArray(X, covlist)
    fit = langevinUD(locs=X, times=times, ID=ID, grad_array=gradArray)
    
    params[ik*3+jk-3, 1:3] = fit$betaHat
    params[ik*3+jk-3, 3] = fit$gamma2Hat
    params[ik*3+jk-3, 5] = thin
  }
  
  print(ik)
  
  ik = ik + 1
}

df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], thin = as.factor(params[,5]))
save(df,file="varying_thin_estimates.Rda")

## plotting estimates ##

p1 <- ggplot(data = df, aes(x = thin, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = "Beta_1") +
  theme_bw()

p2 <- ggplot(data = df, aes(x = thin, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = "Beta_2") +
  theme_bw()

p3 <- ggplot(data = df, aes(x = thin, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = "Beta_3") +
  theme_bw()

p4 <- ggplot(data = df, aes(x = thin, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = "gamma^2") +
  theme_bw()


grid.arrange(p1,p2,p3,p4)

