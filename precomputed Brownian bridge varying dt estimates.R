library(Rhabit)
library(raster)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(ambient)
library(mvnfast)


#defining covariates using perlin noise
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


#squared distance to origin
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))


#using bilinear interpolation to find gradient at matrix of locations
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



#speed parameter for Langevin model
speed <- 5


#estimated parameters
params = matrix(NA, ncol = 5, nrow = 3*100)


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
    
    
    
    times <- seq(0, Tmax, by = dt)
    
    # Generate tracks 
    X = sim_langevin(beta, speed, dt, Tmax, c(0,0), covlist)
    
    # Add ID column
    ID = rep(1, length(time))
    
    
    #thinning tracks
    n = nrow(X)
    X = X[(0:(n%/%thin -1))*thin +1, ]

    
    locs = X
    times = alldat[[1]]$t[(0:(n%/%thin -1))*thin +1]
    ID = alldat[[1]]$ID[(0:(n%/%thin -1))*thin +1]
    
    
    #finding estimate for speed parameter used to generate bridges
    gradArray = bilinearGradArray(X, covlist)
    gammasq = langevinUD(locs=X, times=times, ID=ID, grad_array=gradArray)$gamma2Hat
    
    
    
    #covariance matrix of Brownian bridges
    sigma = matrix(nrow = N, ncol = N)
    for (k in 1:N) {
      for (m in 1:k) {
        sigma[k,m] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
        sigma[m,k] = gammasq*delta*(1 - k/(N+1))*(m/(N+1))
      }
    }
    
    
    #interpolating between points to find mean of brownian bridges
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    
    #array of bridges
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #array of importance sampling wights
    P <- array(data = NA, c(nrow(X)-1, M))
    
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_x, sigma = sigma)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, mu_y, sigma = sigma)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], mu_x, sigma = sigma) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], mu_y, sigma = sigma))
    }
    
    
    
    #find the gradient at the bridge nodes
    Grad <- array(data = NA, c(ncov +1,nrow(X)-1, M, N, 2))
    for (i in 1:(nrow(X)-1)) {
      for(j in 1:M){
        grads <- bilinearGradVec(cbind(B[1, i, j, ], B[2, i, j, ]), covlist)
        
        for (k in 1:(ncov+1)) {
          Grad[k, i, j, 1:N, 1] = grads[k , ,1]
          Grad[k, i, j, 1:N, 2] = grads[k , ,2]
        }
      }
    }
    
    
    #paralellized function to compute likelihood and gradient
    lik_grad <- function(par, cl){
      #log-likelihood
      l = 0
      lik_grad = c(0,0,0,0)
      
      #compute log likelihood and gradient between to observations
      compute <- function(i){
        L_k = 1
        lik_grad_k = 0
        # Add endpoints to all samples (M x (N+2) matrices)
        x_samples = B[1, i, , ]
        y_samples = B[2, i, , ]
        
        #gradient at observations
        grad_0 <- array(data = t(gradArray[i, , ]), c(3,1,2))
        u_0 <- (delta*par[4]/((N+1)*2)) * 
          (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
        
        
        #arrays of postion including bridges and observations
        full_x <- cbind(X[i,1], x_samples, X[i+1,1])
        full_y <- cbind(X[i,2], y_samples, X[i+1,2])
        
        #log likelihood of all locations in bridges
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
        
        #log likelihood gradient of all locations in bridges
        lik_grad_k <- sapply(seq(M), function(j){
          grads <- Grad[ , i, j, , ]
          us <- (delta*par[4]/((N+1)*2)) * 
            (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
          us <- rbind(u_0, us)
          
          g = cbind(grad_0[,1,],array(aperm(Grad[ , i, j, , ], c(1, 3, 2)), dim = c(3, 2 * N)))
          
          D = matrix(t((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                          cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us), ncol = 1)
          
          rbind(g%*%D, -(N+1)/par[4] + (N+1)/(2*delta*par[4]^2)*t(D)%*%D + (1/(2*par[4]))* t(t(g)%*%par[1:3]) %*% D)
          
        })
        
        return(c(log(sum(L_k/M)), -sum(lik_grad_k[1, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[2, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[3, ]*L_k)/(2*sum(L_k)), -sum(lik_grad_k[4, ]*L_k)/(sum(L_k))))
      }
      
      results <- parLapply(cl, 1:(nrow(X)-1), compute)
      
      l = sum(unlist(results)[(1:(nrow(X)-1))*5 -4])
      lik_grad[1] = sum(unlist(results)[(1:(nrow(X)-1))*5 -3])
      lik_grad[2] = sum(unlist(results)[(1:(nrow(X)-1))*5 -2])
      lik_grad[3] = sum(unlist(results)[(1:(nrow(X)-1))*5 -1])
      lik_grad[4] = sum(unlist(results)[(1:(nrow(X)-1))*5])
      
      
      
      if(is.nan(l)){
        print("NaN")
        return(list(l = 1e10, g = c(0,0,0,0)))
      }
      
      if(is.infinite(l)){
        print("Inf")
        return(list(l = 1e10, g = c(0,0,0,0)))
      }else{
        #print(l)
        return(list(l = -l, g = lik_grad))
      }
    }
    
    #using paralellized likelihood in optim
    cl <- makeCluster(12)
    clusterExport(cl, c("X", "M", "N", "delta", "P", "B", "Grad", "gradArray", "dmvn",  "lik_grad"))

    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B")

    stopCluster(cl)
    
    print(o$par)
    params[ik*3+jk-3, 1:4] = o$par
    params[ik*3+jk-3, 5] = thin
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], dt = as.factor(dt*params[,5]))
  save(df,file="varying_thin_estimates.Rda")
  

  print(ik)
  
  ik = ik + 1
}



## plotting estimates ##
p1 <- ggplot(data = df, aes(x = dt, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = expression(beta[1]), y = expression(beta[1]), x = expression(Delta[t])) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()

p1

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




