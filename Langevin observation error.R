library(Rhabit)
library(raster)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
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



###################
## Simulate data ##
###################
beta <- c(4,2,-0.1)
#max time for track
Tmax <- 2500
#increment between times
dt <- 0.01
#speed parameter for Langevin model
speed <- 5
#covariate coefficients
beta <- c(4,2,-0.1)



sim_langevin_error <- function(beta, gammasq, r, dt, Tmax, loc_0, covlist){
  sigma = sqrt(gammasq*dt)
  n_obs = Tmax/dt
  locs = rbind(matrix(loc_0, nrow=1), mvnfast::rmvn(n_obs-1, c(0,0), sigma = diag(sigma, 2,2), isChol = TRUE)) 
  for (i in 1:(n_obs-1)) {
    locs[i+1, ] = locs[i+1, ] + locs[i, ] + (gammasq*dt/2)*bilinearGrad(locs[i, ], covlist)%*%beta 
  }
  
  locs = locs + mvnfast::rmvn(n_obs, c(0,0), sigma = diag(r, 2, 2), isChol = TRUE)
  
  return(locs)
}

X = sim_langevin_error(beta, speed, 0.05, dt, Tmax, c(0,0), covlist)



#thinning data
thin = 50
n = nrow(X)
X = X[(0:(n%/%thin -1))*thin +1, ]





N = 100
M = 200
delta  = dt*thin


sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
  t(lower.tri(sigma_matrix) * sigma_matrix)
chol_m = (chol(sigma_matrix))



sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
  t(lower.tri(sigma_matrix) * sigma_matrix)
chol_m = (chol(sigma_matrix))


#simulates observation error
E = mvnfast::rmvn(nrow(X), matrix(c(0,0), nrow = 1), sigma = diag(1, 2, 2), isChol = TRUE)

Z_x = mvnfast::rmvn(nrow(X), rep(0,N), sigma = chol_m, isChol = TRUE)
Z_y = mvnfast::rmvn(nrow(X), rep(0,N), sigma = chol_m, isChol = TRUE)

lik <- function(par){
  chol_matrix = sqrt(par[4])*chol_m
  
  result = sapply(seq(M), function(j){
    l_k = 0
    
    #l_k = l_k + mvnfast::dmvn(E, matrix(c(0,0), nrow = 1), sigma = diag(sqrt(par[5]), 2, 2), isChol = TRUE, log = TRUE) 

    Y = X + E*sqrt(par[5])
  
  
    mu_x_all <- rep(Y[1:(nrow(Y) - 1), 1], each = N) + 
      1:N * rep((Y[2:nrow(Y), 1] - Y[1:(nrow(Y) - 1), 1]), each = N) / (N + 1)
    mu_y_all <- rep(Y[1:(nrow(Y) - 1), 2], each = N) + 
      1:N * rep((Y[2:nrow(Y), 2] - Y[1:(nrow(Y) - 1), 2]), each = N) / (N + 1)
  
  
    x_samples <- Z_x*sqrt(par[4])
    y_samples <- Z_y*sqrt(par[4])
  
    l_k = l_k - mvnfast::dmvn(x_samples, rep(0,N), sigma = chol_matrix, isChol = TRUE, log = TRUE) -
      mvnfast::dmvn(y_samples, rep(0,N), sigma = chol_matrix, isChol = TRUE, log = TRUE)
  

    full_track = matrix(data = NA, nrow = nrow(Y)*(N+1) - N, ncol = 2)
    full_track[1, ] = Y[1, ]
  
  
    for (i in 1:(nrow(Y) - 1)) {
    
    
      full_track[((N+1)*(i-1)+2):((N+1)*(i)), 1] = x_samples[i, ] + mu_x_all[(N*(i-1)+1):(N*i)]
      full_track[((N+1)*(i-1)+2):((N+1)*(i)), 2] = y_samples[i, ] + mu_y_all[(N*(i-1)+1):(N*i)]
    
      full_track[(N+1)*i+1, ] = Y[i+1, ]
    }
    

    g = bilinearGradVec(full_track[1:(nrow(full_track)-1), ], covlist)
  

    l_k = l_k + sum(mvnfast::dmvn(full_track[2:(nrow(full_track)), ] - full_track[1:(nrow(full_track)-1), ] - (delta*par[4]/2)*(par[1]*g[1, ,] + par[2]*g[2, , ]+ par[3]*g[3, , ]), c(0,0), sigma = diag(par[4],2,2), log = TRUE))
    
    l_k - log(M)
  })
  
  
  max_p = max(result)
  return(-(max_p + log(sum(exp(result - max_p)))))
  
}


t1 = Sys.time()
lik(c(4,2,-0.1,5,0.05))
Sys.time() - t1


cl <- makeCluster(12)     # set the number of processor cores
clusterExport(cl, varlist = c("rmvn", "dmvn", "bilinearGradVec", "X", "chol_m", "M", "N", "delta", "covlist", "E", "Z_x", "Z_y"))
setDefaultCluster(cl=cl) # set 'cl' as default cluster
t1 = Sys.time()
optimParallel(par=c(0,0,0,1,0.05), fn=lik, lower=c(-Inf, -Inf, - Inf, 0.0001, 0.0001))
Sys.time() - t1
setDefaultCluster(cl=NULL); stopCluster(cl)





