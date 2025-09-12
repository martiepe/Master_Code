library(mvnfast)
library(parallel)
library(ambient)  # perlin noise
library(Rcpp)
Rcpp::sourceCpp("compute_lik_grad_full.cpp")
library(here)


set.seed(123)
#number of cores used in parallel computations
ncores = 20
#speed parameter for Langevin model
speed = 5
#covariate coefficients
beta = c(4,2,-0.1)



#generating perlin covariate
lim <- c(-1, 1, -1, 1)*250
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



simLMM <- function(delta, gamma2, covlist, beta, loc0, nobs){
  x = mvnfast::rmvn(nobs*delta/0.01, rep(0,2), 0.01*gamma2*diag(1,2,2))
  x[1, ] = matrix(c(0,0), nrow = 1)
  for (i in 2:nrow(x)) {
    grad = bilinearGradVec(matrix(x[i-1, 1:2], nrow=1), covlist)
    x[i, ]  = x[i, ] + x[i-1, ] + (0.01*gamma2/2)*beta %*% grad[,1,]
  }
  thin = delta/0.01
  n = nrow(x)
  x = x[(0:(n%/%thin -1))*thin +1, ]
  return(x)
}

# define cpp_path in global environment
cpp_path <- here("compute_lik_grad_full.cpp")
#vectorized and paralellized likelihood and gradient function using analytical gradient and no precomputed 
lik_grad <- function(par, cl){
  #log-likelihood
  l = 0
  lik_grad = c(0,0,0,0)
  
  sigma <- diag(delta * par[4] / (N + 1), 2, 2)
  delta_par <- delta * par[4] / (2 * (N + 1))
  
  gamma = sqrt(par[4])
  chol_matrix = gamma*chol_m
  
  
  
  
  compute <- function(i){
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    x_samples = sweep(B[1, i, 1:M, 1:N]*gamma, 2, mu_x, "+")
    y_samples = sweep(B[2, i, 1:M, 1:N]*gamma, 2, mu_y, "+")
    
    L_k = P[i, 1:M]*gamma^(2*N)
    
    full_x <- cbind(X[i,1], x_samples, X[i+1,1])
    full_y <- cbind(X[i,2], y_samples, X[i+1,2])
    
    # CALL THE RCPP FUNCTION:
    # Note: i is 1-based here; the C++ expects the i index 1-based and uses X internally.
    res_mat <- compute_lik_grad_full_cpp(full_x, full_y, L_k, X, i, par, delta, N, covlist)
    # res_mat is 5 x M as in your previous lik_grad_k
    
    L_k_new <- res_mat[1, ]
    # now compute the return vector same as before:
    loglike_i <- log(sum(L_k_new / M))
    # gradients: rows 2..4 correspond to g%*%D components; row 5 is the term3
    g2 = -sum(res_mat[2, ] * L_k_new) / (2 * sum(L_k_new))
    g3 = -sum(res_mat[3, ] * L_k_new) / (2 * sum(L_k_new))
    g4 = -sum(res_mat[4, ] * L_k_new) / (2 * sum(L_k_new))
    g5 = -sum(res_mat[5, ] * L_k_new) / (sum(L_k_new))
    
    return(c(loglike_i, g2, g3, g4, g5))
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



print("varying delta_t, fixed number of observations")
#varying delta_t, fixed number of observations
params = matrix(NA, ncol = 6, nrow = 5*100)
for (ik in 1:100) {
  for (jk in 1:5) {
    beta <- c(4,2,-0.1)
    thin = c(5, 10, 20, 50, 100)[jk]
    dt = 0.01
    delta = dt*thin
    N = thin-1
    M = 50
    n_obs = 5000
    Tmax = n_obs*thin*dt
    
    #simulating track
    X = simLMM(delta,speed,covlist,beta,c(0,0),n_obs)
    
    
    #brownian bridge covariance matrix
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m = (chol(sigma_matrix))
    
    
    #brownian bridge endpoints
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    #brownian bridge array
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      # Generate all M sample tracks at once
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE))
    }
    
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    
    
    clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                  "chol_m",
                                  "covlist", "bilinearGradVec", "delta", "B", "P"), envir = environment())
    clusterEvalQ(cl, library(mvnfast))
    cpp_path <- here("compute_lik_grad_full.cpp")
    
    # export the variable cpp_path (the name as a string)
    clusterExport(cl, varlist = "cpp_path")
    
    # load Rcpp and source the file on all workers
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
    
    
    
    t = Sys.time()
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    print(delta)
    print(t)
    print(o$convergence)
    print(o$counts)
    print(o$par)
    params[ik*5+jk-5, 1:4] = o$par
    params[ik*5+jk-5, 5] = delta
    params[ik*5+jk-5,6] = t
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], dt = as.factor(delta), time = params[,6])
  save(df,file="varying_thin_estimates.Rda")
  
  
  print(ik)
  
  ik = ik + 1
}



print("varying delta_t, fixed maximum time")
#varying delta_t, fixed maximum time
params = matrix(NA, ncol = 6, nrow = 5*100)
for (ik in 1:100) {
  for (jk in 1:5) {
    beta <- c(4,2,-0.1)
    thin = c(5, 10, 20, 50, 100)[jk]
    dt = 0.01
    Tmax = 500
    
    delta = dt*thin
    N = thin-1
    M = 50
    n_obs = Tmax/(dt*thin)
    
    
    #simulating track
    X = simLMM(delta,speed,covlist,beta,c(0,0),n_obs)
    
    
    #brownian bridge covariance matrix
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m = (chol(sigma_matrix))
    
    
    #brownian bridge endpoints
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    #brownian bridge array
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      # Generate all M sample tracks at once
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE))
    }
    
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                  "chol_m",
                                  "covlist", "bilinearGradVec", "delta", "B", "P"), envir = environment())
    clusterEvalQ(cl, library(mvnfast))
    cpp_path <- here("compute_lik_grad_full.cpp")
    
    # export the variable cpp_path (the name as a string)
    clusterExport(cl, varlist = "cpp_path")
    
    # load Rcpp and source the file on all workers
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
    
    
    t = Sys.time()
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    print(delta)
    print(t)
    print(o$convergence)
    print(o$counts)
    print(o$par)
    params[ik*5+jk-5, 1:4] = o$par
    params[ik*5+jk-5, 5] = delta
    params[ik*5+jk-5,6] = t
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], dt = as.factor(delta), time = params[,6])
  save(df,file="varying_thin_estimates_fixed_Tmax.Rda")
  
  
  print(ik)
  
  ik = ik + 1
}




print("varying M")
#varying M
params = matrix(NA, ncol = 6, nrow = 5*100)
for (ik in 1:100) {
  beta <- c(4,2,-0.1)
  thin = 100
  dt = 0.01
  delta = dt*thin
  N = 49
  n_obs = 5000
  Tmax = n_obs*thin*dt
  
  #simulating track
  X = simLMM(delta,speed,covlist,beta,c(0,0),n_obs)
  
  
  for (jk in 1:5) {
    M = c(5, 10, 50, 100, 200)[jk]
    
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m = (chol(sigma_matrix))
    
    
    #brownian bridge endpoints
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    #brownian bridge array
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      # Generate all M sample tracks at once
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE))
    }
    
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                  "chol_m",
                                  "covlist", "bilinearGradVec", "delta", "B", "P"), envir = environment())
    clusterEvalQ(cl, library(mvnfast))
    cpp_path <- here("compute_lik_grad_full.cpp")
    
    # export the variable cpp_path (the name as a string)
    clusterExport(cl, varlist = "cpp_path")
    
    # load Rcpp and source the file on all workers
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
    
    
    t = Sys.time()
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    print(M)
    print(t)
    print(o$convergence)
    print(o$counts)
    print(o$par)
    params[ik*5+jk-5, 1:4] = o$par
    params[ik*5+jk-5, 5] = M
    params[ik*5+jk-5,6] = t
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], M = as.factor(params[, 5]), time = params[, 6])
  save(df,file="varying_M_estimates_stochastic likelihood.Rda")
  
  
  print(ik)
}


print("varying N")
#varying N
params = matrix(NA, ncol = 6, nrow = 4*100)
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
    
    
    
    #simulating track
    X = simLMM(delta,speed,covlist,beta,c(0,0),n_obs)
    
    
    sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m = (chol(sigma_matrix))
    
    
    #brownian bridge endpoints
    mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
    mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
    #brownian bridge array
    B <- array(data = NA, c(2, nrow(X)-1, M, N))
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
      mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
      
      # Generate all M sample tracks at once
      B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
      
      
      P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE) * 
                       mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE))
    }
    
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                  "chol_m",
                                  "covlist", "bilinearGradVec", "delta", "B", "P"), envir = environment())
    clusterEvalQ(cl, library(mvnfast))
    cpp_path <- here("compute_lik_grad_full.cpp")
    
    # export the variable cpp_path (the name as a string)
    clusterExport(cl, varlist = "cpp_path")
    
    # load Rcpp and source the file on all workers
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
    
    
    t = Sys.time()
    o = optim(par = c(0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    print(N)
    print(t)
    print(o$convergence)
    print(o$counts)
    print(o$par)
    params[ik*5+jk-5, 1:4] = o$par
    params[ik*5+jk-5, 5] = N
    params[ik*5+jk-5,6] = t
    
    
  }
  
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], gammasq = params[,4], N = as.factor(params[,5]), time = params[,6])
  save(df,file="varying_N_estimates.Rda")
  
  
  print(ik)
}

