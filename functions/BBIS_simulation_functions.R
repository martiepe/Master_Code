# simulate langevin track (High-res)
sumLTrack <- function(delta, gamma2, covlist, beta, loc0, nobs) {
  x = mvnfast::rmvn(nobs*delta/0.01, rep(0,2), 0.01*gamma2*diag(1,2,2))
  x[1, ] <- loc0 
  for (i in 2:nrow(x)) {
    grad = bilinearGradVec(matrix(x[i-1, 1:2], nrow=1), covlist)
    x[i, ]  = x[i, ] + x[i-1, ] + (0.01*gamma2/2)*beta %*% grad[,1,]
  }
  return(x)
}

# thin high-resolution track
thinTrack <- function(x, delta) {
  thin = delta/0.01
  n = nrow(x)
  x = x[(0:(n%/%thin -1))*thin +1, ]
  return(x)
}

# simulate thinned Langevin movement model
simLMM <- function(delta, gamma2, covlist, beta, loc0, nobs){ 
  x <- sumLTrack(delta, gamma2, covlist, beta, loc0, nobs)  # Simulate high-resolution track
  x <- thinTrack(x, delta) # Thin the track
  return(x)
}

# define cpp_path in global environment
library(Rcpp)
Rcpp::sourceCpp("compute_lik_grad_full.cpp")

# bilinear interpolation gradient function
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

# vectorized and paralellized likelihood and gradient function using analytical gradient and no precomputed 
# Generalized likelihood function
lik_grad <- function(par, cl, n_cov){
  # par has length n_cov + 1: first n_cov are beta coefficients, last is variance
  
  l = 0
  lik_grad = rep(0, n_cov + 1)
  
  sigma <- diag(delta * par[n_cov + 1] / (N + 1), 2, 2)
  delta_par <- delta * par[n_cov + 1] / (2 * (N + 1))
  
  gamma = sqrt(par[n_cov + 1])
  chol_matrix = gamma * chol_m
  
  compute <- function(i){
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    x_samples = sweep(B[1, i, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples = sweep(B[2, i, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    L_k = P[i, 1:M] * gamma^(2 * N)
    
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])
    
    # Call C++ function
    # res_mat is (n_cov + 2) x M: 1 row for likelihood, n_cov rows for beta gradients, 1 row for variance gradient
    
    
    
    
    
    
    
    
    res_mat <- compute_lik_grad_full_cpp(full_x, full_y, L_k, X, i, par, delta, N, covlist)
    
    L_k_new <- res_mat[1, ]
    loglike_i <- log(sum(L_k_new / M))
    
    # Compute gradients for all beta coefficients
    grad_beta <- numeric(n_cov)
    for(j in 1:n_cov) {
      grad_beta[j] = -sum(res_mat[j + 1, ] * L_k_new) / (2 * sum(L_k_new))
    }
    
    # Gradient for variance parameter
    grad_sigma = -sum(res_mat[n_cov + 2, ] * L_k_new) / (sum(L_k_new))
    
    return(c(loglike_i, grad_beta, grad_sigma))
  }
  
  results <- parLapply(cl, 1:(nrow(X) - 1), compute)
  
  # Extract results
  results_mat <- do.call(rbind, results)
  l = sum(results_mat[, 1])
  for(j in 1:(n_cov + 1)) {
    lik_grad[j] = sum(results_mat[, j + 1])
  }
  
  if(is.nan(l)){
    print("NaN")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  }
  
  if(is.infinite(l)){
    print("Inf")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  } else {
    return(list(l = -l, g = lik_grad))
  }
}

# Wrapper function to fit the model
fit_langevin_bbis <- function(X, covlist, delta, N = 4, M = 50, ncores = 10, 
                              init_beta = NULL, init_sigma = 1,
                              lower_sigma = 0.0001, cpp_path = NULL, verpose = TRUE) {
  # prep covariates
  if (inherits(covlist, "SpatRaster")) {
    # convert from spatraster
    covlist <- spatRast_to_covlist(covlist)
  } 
  
  # Determine number of covariates
  n_cov <- length(covlist)
  
  # Set initial parameters if not provided
  if(is.null(init_beta)) {
    init_beta <- rep(0, n_cov)
  }
  
  if(length(init_beta) != n_cov) {
    stop("init_beta length must match number of covariates")
  }
  
  par_init <- c(init_beta, init_sigma)
  
  # Brownian bridge covariance matrix
  sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
  sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
    t(lower.tri(sigma_matrix) * sigma_matrix)
  chol_m <- chol(sigma_matrix)
  
  # Brownian bridge endpoints
  mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 
    1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
  mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 
    1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
  
  # Brownian bridge array
  B <- array(data = NA, c(2, nrow(X)-1, M, N))
  # Importance sampling weights
  P <- array(data = NA, c(nrow(X)-1, M))
  
  # Generate bridges
  if(verpose) cat("Generating Brownian bridges...\n")
  for (i in 1:(nrow(X)-1)) {
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
    B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
    
    P[i, 1:M] <- 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE) * 
                      mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE))
  }
  
  # Set up parallel cluster
  if(verpose) cat("Setting up parallel cluster with", ncores, "cores...\n")
  cl <- makeCluster(ncores)
  
  # Export variables to cluster
  clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                "chol_m", "covlist", "delta", "B", "P", "n_cov"), 
                envir = environment())
  
  # Load required libraries on workers
  clusterEvalQ(cl, library(mvnfast))
  
  # Source C++ file if provided
  if(!is.null(cpp_path)) {
    clusterExport(cl, varlist = "cpp_path")
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
  }
  
  # Run optimization
  if(verpose) cat("Running optimization...\n")
  t_start <- Sys.time()
  
  # Set bounds: no bounds on betas, lower bound on sigma
  lower_bounds <- c(rep(-Inf, n_cov), lower_sigma)
  
  o <- optim(par = par_init, 
             fn = function(x) lik_grad(x, cl, n_cov)$l, 
             gr = function(x) lik_grad(x, cl, n_cov)$g, 
             method = "L-BFGS-B", 
             lower = lower_bounds)
  
  t_elapsed <- Sys.time() - t_start
  
  # Stop cluster
  stopCluster(cl)
  
  # Prepare results
  results <- list(
    par = o$par,
    beta = o$par[1:n_cov],
    sigma = o$par[n_cov + 1],
    convergence = o$convergence,
    value = o$value,
    counts = o$counts,
    time = t_elapsed,
    delta = delta,
    N = N,
    M = M,
    n_cov = n_cov
  )
  
  # Print summary
  if(verpose) {
    cat("\n--- Optimization Results ---\n")
    cat("Time elapsed:", format(t_elapsed), "\n")
    cat("Convergence:", o$convergence, "\n")
    cat("Function evaluations:", o$counts[1], "\n")
    cat("Gradient evaluations:", o$counts[2], "\n")
    cat("\nEstimated parameters:\n")
    for(i in 1:n_cov) {
      cat(sprintf("  beta[%d] = %.4f\n", i, results$beta[i]))
    }
    cat(sprintf("  sigma   = %.4f\n", results$sigma))
  }
  
  return(results)
}














###### For varying time, multiple tracks ##############


# vectorized and paralellized likelihood and gradient function using analytical gradient and no precomputed 
# Generalized likelihood function
lik_grad <- function(par, cl, n_cov){
  # par has length n_cov + 1: first n_cov are beta coefficients, last is variance
  
  l = 0
  lik_grad = rep(0, n_cov + 1)
  
  sigma <- diag(delta * par[n_cov + 1] / (N + 1), 2, 2)
  delta_par <- delta * par[n_cov + 1] / (2 * (N + 1))
  
  gamma = sqrt(par[n_cov + 1])
  chol_matrix = gamma * chol_m
  
  compute <- function(i){
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    x_samples = sweep(B[1, i, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples = sweep(B[2, i, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    L_k = P[i, 1:M] * gamma^(2 * N)
    
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])
    
    # Call C++ function
    # res_mat is (n_cov + 2) x M: 1 row for likelihood, n_cov rows for beta gradients, 1 row for variance gradient
    
    
    
    
    
    
    
    
    res_mat <- compute_lik_grad_full_cpp(full_x, full_y, L_k, X, i, par, delta, N, covlist)
    
    L_k_new <- res_mat[1, ]
    loglike_i <- log(sum(L_k_new / M))
    
    # Compute gradients for all beta coefficients
    grad_beta <- numeric(n_cov)
    for(j in 1:n_cov) {
      grad_beta[j] = -sum(res_mat[j + 1, ] * L_k_new) / (2 * sum(L_k_new))
    }
    
    # Gradient for variance parameter
    grad_sigma = -sum(res_mat[n_cov + 2, ] * L_k_new) / (sum(L_k_new))
    
    return(c(loglike_i, grad_beta, grad_sigma))
  }
  

  
  results <- parLapply(cl, 1:(nrow(X) - 1), compute)
  
  # Extract results
  results_mat <- do.call(rbind, results)
  l = sum(results_mat[, 1])
  for(j in 1:(n_cov + 1)) {
    lik_grad[j] = sum(results_mat[, j + 1])
  }
  
  if(is.nan(l)){
    print("NaN")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  }
  
  if(is.infinite(l)){
    print("Inf")
    return(list(l = 1e10, g = rep(0, n_cov + 1)))
  } else {
    return(list(l = -l, g = lik_grad))
  }
}
















# Wrapper function to fit the model
fit_langevin_bbis <- function(X, covlist, delta_max, N = 4, M = 50, ncores = 10, 
                              init_beta = NULL, init_sigma = 1,
                              lower_sigma = 0.0001, cpp_path = NULL, verpose = TRUE) {
  # prep covariates
  if (inherits(covlist, "SpatRaster")) {
    # convert from spatraster
    covlist <- spatRast_to_covlist(covlist)
  } 
  
  # Determine number of covariates
  n_cov <- length(covlist)
  
  # Set initial parameters if not provided
  if(is.null(init_beta)) {
    init_beta <- rep(0, n_cov)
  }
  
  if(length(init_beta) != n_cov) {
    stop("init_beta length must match number of covariates")
  }
  
  par_init <- c(init_beta, init_sigma)
  
  #number of observations in each track
  n = c()
  for (k in unique(X$ID)) {
    n = c(n, sum(X$ID == k))
  }
  
  
  #transition indexes
  transitions = 1:(n[1]-1)
  for (k in 2:length(n)) {
    transitions = c(transitions, (n[k-1]-1):(n[k]-1))
  }
  
  #finding number of nodes needed for each transition
  N = ceiling((X$time[transitions+1] - X$time[transitions])/delta_max)
  
  #time step between the nodes between each transition
  delta = (X$time[transitions+1] - X$time[transitions])/N
  
  #generating brownian bridges for each transition
  
  
  
  
  
  
  
  
  
  
  
  # Brownian bridge covariance matrix
  sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
  sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
    t(lower.tri(sigma_matrix) * sigma_matrix)
  chol_m <- chol(sigma_matrix)
  
  # Brownian bridge endpoints
  mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 
    1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
  mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 
    1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
  
  # Brownian bridge array
  B <- array(data = NA, c(2, nrow(X)-1, M, N))
  # Importance sampling weights
  P <- array(data = NA, c(nrow(X)-1, M))
  
  # Generate bridges
  if(verpose) cat("Generating Brownian bridges...\n")
  for (i in 1:(nrow(X)-1)) {
    mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
    mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
    
    B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
    B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0, N), sigma = chol_m, isChol = TRUE)
    
    P[i, 1:M] <- 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE) * 
                      mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0, N), sigma = chol_m, isChol = TRUE))
  }
  
  # Set up parallel cluster
  if(verpose) cat("Setting up parallel cluster with", ncores, "cores...\n")
  cl <- makeCluster(ncores)
  
  # Export variables to cluster
  clusterExport(cl, varlist = c("X", "N", "M", "mu_x_all", "mu_y_all",
                                "chol_m", "covlist", "delta", "B", "P", "n_cov"), 
                envir = environment())
  
  # Load required libraries on workers
  clusterEvalQ(cl, library(mvnfast))
  
  # Source C++ file if provided
  if(!is.null(cpp_path)) {
    clusterExport(cl, varlist = "cpp_path")
    clusterEvalQ(cl, {
      library(Rcpp)
      sourceCpp(cpp_path)
    })
  }
  
  # Run optimization
  if(verpose) cat("Running optimization...\n")
  t_start <- Sys.time()
  
  # Set bounds: no bounds on betas, lower bound on sigma
  lower_bounds <- c(rep(-Inf, n_cov), lower_sigma)
  
  o <- optim(par = par_init, 
             fn = function(x) lik_grad(x, cl, n_cov)$l, 
             gr = function(x) lik_grad(x, cl, n_cov)$g, 
             method = "L-BFGS-B", 
             lower = lower_bounds)
  
  t_elapsed <- Sys.time() - t_start
  
  # Stop cluster
  stopCluster(cl)
  
  # Prepare results
  results <- list(
    par = o$par,
    beta = o$par[1:n_cov],
    sigma = o$par[n_cov + 1],
    convergence = o$convergence,
    value = o$value,
    counts = o$counts,
    time = t_elapsed,
    delta = delta,
    N = N,
    M = M,
    n_cov = n_cov
  )
  
  # Print summary
  if(verpose) {
    cat("\n--- Optimization Results ---\n")
    cat("Time elapsed:", format(t_elapsed), "\n")
    cat("Convergence:", o$convergence, "\n")
    cat("Function evaluations:", o$counts[1], "\n")
    cat("Gradient evaluations:", o$counts[2], "\n")
    cat("\nEstimated parameters:\n")
    for(i in 1:n_cov) {
      cat(sprintf("  beta[%d] = %.4f\n", i, results$beta[i]))
    }
    cat(sprintf("  sigma   = %.4f\n", results$sigma))
  }
  
  return(results)
}










