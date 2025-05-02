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



params = matrix(NA, ncol = 5, nrow = 4*100)
for (ik in 1:100) {
  for (jk in 1:4) {
    beta <- c(4,2,-0.1)
    thin = 100
    dt = 0.01
    delta = dt*thin
    N = c(4,9,49,99)[jk]
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
    
    
    #thinning tracks
    X = matrix(c(alldat[[1]]$x, alldat[[1]]$y), ncol = 2)
    n = nrow(X)
    X = X[(0:(n%/%thin -1))*thin +1, ]
    
    
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m = (chol(sigma_matrix))
    
    
    mu_x_all <- rep(X[1:(nrow(X) - 1), 1], each = N) + 
      1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X) - 1), 1]), each = N) / (N + 1)
    mu_y_all <- rep(X[1:(nrow(X) - 1), 2], each = N) + 
      1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X) - 1), 2]), each = N) / (N + 1)
    
    
    #vectorized and paralellized likelihood and gradient function using analytical gradient and no precomputed 
    lik_grad <- function(par, cl){
      #log-likelihood
      l = 0
      lik_grad = c(0,0,0,0)
      
      sigma <- diag(delta * par[4] / (N + 1), 2, 2)
      delta_par <- delta * par[4] / (2 * (N + 1))
      
      
      gamma = sqrt(par[4])
      chol_matrix = gamma*chol_m
      
      
      #number of simulations
      clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                    "chol_matrix", "delta_par",
                                    "par", "covlist", "bilinearGradVec", "delta"), envir = environment())
      clusterEvalQ(cl, library(mvnfast))
      
      
      
      compute <- function(i){
        mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
        mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
        
        
        x_samples <- mvnfast::rmvn(M, mu_x, sigma = chol_matrix, isChol = TRUE)
        y_samples <- mvnfast::rmvn(M, mu_y, sigma = chol_matrix, isChol = TRUE)
        
        
        L_k = 1 / (mvnfast::dmvn(x_samples, mu_x, sigma = chol_matrix, isChol = TRUE) * 
                     mvnfast::dmvn(y_samples, mu_y, sigma = chol_matrix, isChol = TRUE))
        
        
        
        grad_0 <- bilinearGradVec(matrix(X[i, 1:2], ncol = 2), covlist)
        
        u_0 <- (delta*par[4]/((N+1)*2)) * 
          (par[1] * grad_0[1,,] + par[2] * grad_0[2,,] + par[3] * grad_0[3,,])
        
        full_x <- cbind(X[i,1], x_samples, X[i+1,1])
        full_y <- cbind(X[i,2], y_samples, X[i+1,2])
        
        
        #grads is used in both of these
        # likelihood of all locations
        L_k <- sapply(seq(M), function(j) {
          grads <- bilinearGradVec(cbind(x_samples[j, ], y_samples[j, ]), covlist)
          us <- (delta*par[4]/((N+1)*2)) * 
            (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
          us <- rbind(u_0, us)
          prod(dmvn((cbind(full_x[j,0:N+2], full_y[j,0:N+2]) - 
                       cbind(full_x[j,0:N+1], full_y[j,0:N+1])) - us, 
                    matrix(c(0,0)),
                    diag(delta*par[4]/(N+1), 2, 2)))
        }) * L_k
        
        
        lik_grad_k <- sapply(seq(M), function(j){
          grads <- bilinearGradVec(cbind(x_samples[j, ], y_samples[j, ]), covlist)
          us <- (delta*par[4]/((N+1)*2)) * 
            (par[1] * grads[1,,] + par[2] * grads[2,,] + par[3] * grads[3,,]) 
          us <- rbind(u_0, us)
          
          g = cbind(grad_0[,1,],array(aperm(grads, c(1, 3, 2)), dim = c(3, 2 * N)))
          
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
      
      
      return(list(l = -l, g = lik_grad))
      
    }
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(12)
    
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0))
    
    stopCluster(cl)
    
    
    
    print(o$par)
    params[ik*4+jk-4, 1:4] = o$par
    params[ik*4+jk-4, 5] = N
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], N = as.factor(params[,5]))
  save(df,file="varying_N_estimates_stochastic_likelihood.Rda")
  
  
  print(ik)
}



p1 <- ggplot(data = df, aes(x = N, y = beta1)) +
  geom_boxplot() +
  geom_hline(yintercept  = 4, color = "red", linetype = 2) +
  labs(title = "Beta_1") +
  theme_bw()

p2 <- ggplot(data = df, aes(x = N, y = beta2)) +
  geom_boxplot() +
  geom_hline(yintercept  = 2, color = "red", linetype = 2) +
  labs(title = "Beta_2") +
  theme_bw()

p3 <- ggplot(data = df, aes(x = N, y = beta3)) +
  geom_boxplot() +
  geom_hline(yintercept  = -0.1, color = "red", linetype = 2) +
  labs(title = "Beta_3") +
  theme_bw()

p4 <- ggplot(data = df, aes(x = N, y = gammasq)) +
  geom_boxplot() +
  geom_hline(yintercept  = 5, color = "red", linetype = 2) +
  labs(title = "gamma^2") +
  theme_bw()


grid.arrange(p1,p2,p3,p4)

