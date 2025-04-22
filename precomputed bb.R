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






###################
## Simulate data ##
###################
beta <- c(4,2,-0.1)
#max time for track
Tmax <- 5000
#increment between times
dt <- 0.01
#time grid
time <- seq(0, Tmax, by = dt)
#number of tracks to be generated
ntrack <- 1
#speed parameter for Langevin model
speed <- 5
#covariate coefficients
beta <- c(4,2,-0.1)
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





#estimate using thinned data to see 
thin = 10
X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
N = nrow(X)
X = X[(0:(N%/%thin -1))*thin +1, ]
#X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)

gradArray = bilinearGradArray(X, covlist)

locs = X
times = alldat[[1]]$t[(0:(N%/%thin -1))*thin +1]
ID = alldat[[1]]$ID[(0:(N%/%thin -1))*thin +1]

fit <- langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)

fit$betaHat
fit$gamma2Hat









###### one brownian bridge importance sampling estimate #####
N = thin-1
M = 50
delta = dt*thin
#find the diffusion constant to be used
gradArray = bilinearGradArray(X, covlist)
locs = X
gammasq = langevinUD(locs=locs, times=times, ID=ID, grad_array=gradArray)$gamma2Hat

#simulate brownian bridges
sigma = matrix(nrow = N, ncol = N)
for (k in 1:N) {
  for (m in 1:k) {
    sigma[k,m] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
    sigma[m,k] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
  }
}




##### Constructing Brownian bridges ############
mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)

B <- array(data = NA, c(2, nrow(X)-1, M, N))
#importance sampling wights
P <- array(data = NA, c(nrow(X)-1, M))
for (i in 1:(nrow(X)-1)) {
  mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
  mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
  
  # Generate all M sample tracks at once
  B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_x, sigma = sigma)
  B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_y, sigma = sigma)
  
  
  P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], mu_x, sigma = sigma) * 
                   mvnfast::dmvn(B[2, i, 1:M, 1:N], mu_y, sigma = sigma))
}



#find the gradient at the bridge nodes
Grad <- array(data = NA, c(ncov +1,nrow(X)-1, M, N, 2))
for (i in 1:(nrow(X)-1)) {
  for(j in 1:M){
    #grad = bilinearGradArray(t(B[1:2, i, j, 1:N]), covlist)
    
    grads <- bilinearGradVec(cbind(B[1, i, j, ], B[2, i, j, ]), covlist)
    
    for (k in 1:(ncov+1)) {
      Grad[k, i, j, 1:N, 1] = grads[k , ,1]
      Grad[k, i, j, 1:N, 2] = grads[k , ,2]
    }
  }
}






#vectorized likelihood
lik <- function(par){
  #log-likelihood
  l = 0
  
  #number of simulations
  
  #for each observation in the track
  for (i in 1:(nrow(X)-1)) {
    L_k = P[i, ]
    
    # Add endpoints to all samples (M x (N+2) matrices)
    # calc initial gradient
    x_samples = B[1, i, , ]
    y_samples = B[2, i, , ]
    
    grad_0 <- array(data = t(gradArray[i, , ]), c(3,1,2))
    
    u_0 <- (delta*par[4]/((N+1)*2)) * 
      (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
    
    full_x <- cbind(X[i,1], x_samples, X[i+1,1])
    full_y <- cbind(X[i,2], y_samples, X[i+1,2])
    # likelihood of all locations
    L_k <- sapply(seq(M), function(j) {
      grads <- Grad[ , i, j, , ]
      us <- (delta*par[4]/((N+1)*2)) * 
        (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
      us <- rbind(u_0, us)
      prod(dmvn((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                   cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us, 
                matrix(c(0,0)),
                diag(delta*par[4]/(N+1), 2, 2)))
    })*L_k
    
    
    l = l + log(sum(L_k/M))
    
  }
  
  return(-l)
  
}








t1 = Sys.time()
lik(c(4,2,-0.1,5))
Sys.time() - t1



t1 = Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("X", "M", "N", "delta", "P", "B", "Grad", "gradArray", "dmvn"))
setDefaultCluster(cl=cl)
o = optimParallel(c(0,0,0,1), lik)
setDefaultCluster(cl=NULL); stopCluster(cl)
Sys.time() - t1
o









#vectorized likelihood and gradient function
lik_grad <- function(par){
  #log-likelihood
  l = 0
  lik_grad = c(0,0,0,0)
  #number of simulations
  
  #for each observation in the track
  for (i in 1:(nrow(X)-1)) {
    L_k = 1
    lik_grad_k = 0
    # Add endpoints to all samples (M x (N+2) matrices)
    # calc initial gradient
    x_samples = B[1, i, , ]
    y_samples = B[2, i, , ]
    
    grad_0 <- array(data = t(gradArray[i, , ]), c(3,1,2))
    
    u_0 <- (delta*par[4]/((N+1)*2)) * 
      (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
    
    full_x <- cbind(X[i,1], x_samples, X[i+1,1])
    full_y <- cbind(X[i,2], y_samples, X[i+1,2])
    # likelihood of all locations
    L_k <- sapply(seq(M), function(j) {
      grads <- Grad[ , i, j, , ]
      us <- (delta*par[4]/((N+1)*2)) * 
        (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
      us <- rbind(u_0, us)
      prod(dmvn((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                   cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us, 
                matrix(c(0,0)),
                diag(delta*par[4]/(N+1), 2, 2)))
      
      
    })*P[i, ]
    
    lik_grad_k <- sapply(seq(M), function(j){
      grads <- Grad[ , i, j, , ]
      us <- (delta*par[4]/((N+1)*2)) * 
        (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
      us <- rbind(u_0, us)
      
      g = (delta*par[4]/((N+1)*2))*cbind(grad_0[,1,],array(aperm(Grad[ , i, j, , ], c(1, 3, 2)), dim = c(3, 2 * N)))
      
      D = matrix(t((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                   cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us), ncol = 1)

      rbind(g%*%D, 2/par[4] + (delta/2)*t(D)%*%D)
      
    })
    
    lik_grad[1] = lik_grad[1] - delta*par[4]*sum(lik_grad_k[1, ]*L_k*P[i, ])/(2*M*sum(L_k))
    lik_grad[2] = lik_grad[2] - delta*par[4]*sum(lik_grad_k[2, ]*L_k*P[i, ])/(2*M*sum(L_k))
    lik_grad[3] = lik_grad[3] - delta*par[4]*sum(lik_grad_k[3, ]*L_k*P[i, ])/(2*M*sum(L_k))
    
    lik_grad[4] = lik_grad[4] - sum(lik_grad_k[1, ]*L_k*P[i, ])/(M*sum(L_k))
    
    
    l = l + log(sum(L_k/M))
    
  }
  
  return(list(l = -l, g = lik_grad))
  
}

lik_grad(c(4,2,-0.1,5))


t1 = Sys.time()
lik_grad(c(4,2,-0.1,5))
Sys.time() - t1

# -1.717325 -14.472553 -27.748698 114.346584
t1 = Sys.time()
grad(lik, x = c( 4.02595981,  2.13715837, -0.07482453,  4.95182969))
Sys.time() - t1




t1 = Sys.time()
lik_grad(c(4,2,-0.1,5))
Sys.time() - t1


t1 = Sys.time()
lik(c(4,2,-0.1,5))
Sys.time() - t1

library(numDeriv)



# Define the function and gradient
func_and_grad <- function(x) {
  f_val <- sum(x^2)  # Example function: sum of squares
  grad_val <- 2 * x   # Gradient: derivative of sum of squares
  list(value = f_val, gradient = grad_val)
}

# Use optim
result <- optim(par = c(1, 2), fn = function(x) func_and_grad(x)$value, gr = function(x) func_and_grad(x)$gradient)


result














