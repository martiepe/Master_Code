# Case study
# prep workspace ####
library(here)
library(mvnfast)
library(parallel)
library(terra)
library(dplyr)
library(Rhabit)
library(raster)
# custom functions
source(here("functions/utility_functions.R"))
sourceDir("functions")


# explore
library(tidyterra)
library(ggplot2)
ggplot() + 
  geom_spatraster(data= habitat) +
  facet_wrap(~lyr) +
  geom_point(data = tracks, aes(x,y))



####################
## Prepare tracks ##
####################
# Load regularised track (from SSLpreprocessing.R)
tracks <- read.csv(here("post-submission/case_study/SSLpreddat.csv"))
ID <- as.integer(tracks$ID)
times <- as.POSIXct(tracks$time)
times <- as.numeric(times)
times <- (times-min(times))/3600
X <- matrix(c(tracks$x, tracks$y)/1000, ncol=2) # convert to km




#####################
## Load covariates ##
#####################
hbfull <- brick(here("post-submission/case_study/aleut_habitat.grd"), values=TRUE)
covlist0 <- list(bathy = hbfull$bathy,
                 slope = hbfull$slope,
                 d2site = hbfull$d2site,
                 d2shelf = hbfull$d2shelf)

# Convert to km
for(i in 1:length(covlist0)) {
  extent(covlist0[[i]]) <- extent(c(xmin(covlist0[[i]]), xmax(covlist0[[i]]), 
                                    ymin(covlist0[[i]]), ymax(covlist0[[i]]))/1000)
  projection(covlist0[[i]]) <- gsub("units=m", "units=km", projection(covlist0[[i]]))
}

ncov <- length(covlist0)
# Resample covariates to the same grid
for(i in 2:ncov)
  covlist0[[i]] <- resample(covlist0[[i]],covlist0[[1]])

rasterToRhabit <- function(cov) {
  # Extract limits and resolution from raster
  lim <- as.vector(extent(cov))
  res <- res(cov)
  # Define x and y grids
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by=res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by=res[2])
  
  # Put covariate values in the right matrix format
  z <- t(apply(as.matrix(cov),2,rev))
  
  return(list(x = xgrid, y = ygrid, z = z))
}

covlist <- lapply(covlist0, rasterToRhabit)



# Evaluate covariate gradients at observed locations
gradarray <- bilinearGradArray(locs = X, cov_list = covlist)


ncores = detectCores()-1



#calculates likelihood and gradient of using the brownian bridge Langevin likelihood approximation
lik_grad <- function(par, cl){
  #number of beta parameters
  p <- length(par) - 1L
  
  gamma <- sqrt(par[p + 1L])

  #compute likelihood and gradient contribution between two observations 
  compute <- function(i){
    if (ID[i] != ID[i + 1]) return(rep(0, p + 2L))
    
    N     <- ceiling((times[i + 1] - times[i]) / delta_max)
    #N = max(2, N)
    delta <- (times[i + 1] - times[i])
    
    if (N == 1L) {
      ## --- Maruyama one-step shortcut; no IS, no C++ ---
      ## increment
      y <- c(X[i + 1, 1] - X[i, 1],
             X[i + 1, 2] - X[i, 2])
      
      ## gradients at the start location; bilinearGradVec returns (p, n_obs, 2)
      G_arr <- bilinearGradVec(matrix(X[i, ], nrow = 1), covlist)
      ## make a 2 x p with rows (df/dx, df/dy)
      G_i <- t(drop(G_arr[, 1, ]))  # 2 x p
      
      ## if any NA/Inf (e.g., boundary), skip this segment consistently
      if (!all(is.finite(G_i))){
        print("inf")
        return(rep(0, p + 2L))
      }
      
      beta <- par[1:p]
      s    <- gamma^2
      
      mu <- as.vector((delta * s / 2) * (G_i %*% beta))  # length 2
      e  <- y - mu
      
      ## log-likelihood for d=2, Σ = δ s I_2
      loglike_i <- -log(2 * pi) - log(delta * s) - sum(e * e) / (2 * delta * s)
      
      ## scores: wrt beta and s (= gamma^2)
      g_beta <- 0.5 * as.numeric(crossprod(G_i, e))  # p-vector
      g_s    <- -1 / s + (sum(e * e)) / (2 * delta * s^2) + as.numeric(t(e) %*% (G_i %*% beta)) / (2 * s)
      
      if (!is.finite(loglike_i) || !all(is.finite(g_beta)) || !is.finite(g_s)) {
        return(rep(0, p + 2L))
      }
      return(c(loglike_i, -g_beta, -g_s))
    }
    
    
    # brownian bridge endpoints (same construction you had)
    mu_x <- rep(X[i, 1], each = N) + 1:N * rep((X[i + 1, 1] - X[i, 1]), each = N) / (N + 1)
    mu_y <- rep(X[i, 2], each = N) + 1:N * rep((X[i + 1, 2] - X[i, 2]), each = N) / (N + 1)
    
    #scaled bridges
    x_samples <- sweep(B[[i]][1, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples <- sweep(B[[i]][2, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    #bridge contribtution to likelihood
    l_k <- P[i, 1:M]  + log(gamma)*(2 * N)
    
    #bridges with endpoints
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])
    
    #likelihood and gradient calculationsusing a c++ function
    res_mat <- compute_log_lik_grad_full_cpp(full_x, full_y, l_k, X, i, par, delta, N, covlist)
    
    logw <- res_mat[1, ]
    a    <- max(logw)
    w    <- exp(logw - a)
    Wsum <- sum(w)
    if (!is.finite(Wsum) || Wsum <= 0) return(rep(0, p + 2L))
    
    loglike_i <- a + log(Wsum) - log(M)
    
    norm_w <- w / Wsum
    
    g_beta <- numeric(p)
    for (j in 1:p) {
      g_beta[j] <- -sum(res_mat[1 + j, ] * norm_w) / 2
    }
    
    g_s <- -sum(res_mat[p + 2L, ] * norm_w)
    
    c(loglike_i, g_beta, g_s)
    
    
  }
  
  #likelihood and gradient between all observations
  results_list <- parLapply(cl, 1:(nrow(X) - 1L), compute)

  res_mat <- do.call(rbind, results_list)  
  
  l_sum      <- sum(res_mat[, 1], na.rm = TRUE)
  g_betasum  <- colSums(res_mat[, 2:(p + 1L), drop = FALSE], na.rm = TRUE)
  g_ssum     <- sum(res_mat[, p + 2L], na.rm = TRUE)
  
  if (!is.finite(l_sum)) {
    return(list(l = 1e10, g = rep(0, p + 1L)))
  }
  list(l = -l_sum, g = c(g_betasum, g_ssum))
}


#cached function that stores the likelihood and gradient for the last parameters used
make_lik_grad_cached <- function(cl) {
  last_par <- NULL
  last_res <- NULL
  
  eval_lik_grad <- function(par) {
    # recompute only if par changed
    if (is.null(last_par) || length(par) != length(last_par) || any(par != last_par)) {
      last_par <<- par
      last_res <<- lik_grad(par, cl)  # one expensive call
    }
    last_res
  }
  
  fn <- function(par) {
    eval_lik_grad(par)$l  # return likelihood
  }
  
  gr <- function(par) {
    eval_lik_grad(par)$g  # return gradient
  }
  
  list(fn = fn, gr = gr)
}

wrapped <- make_lik_grad_cached(cl)



#sea lion estimates on grid using M=25 bridges
M = 25
deltas = exp(seq(log(0.01), log(25), length.out = 30))
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))
for (k in 1:length(deltas)) {
  for (d in 1:10) {
    delta_max = deltas[k]
    
    #brownian bridge array
    B <- list()
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
    if(ID[i] == ID[i+1]){
      N = ceiling((times[i+1] - times[i])/delta_max)
      #N = max(2, N)
      if(N != 1){
        delta = (times[i+1] - times[i])
        
        #brownian bridge covariance matrix
        sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
        sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
          t(lower.tri(sigma_matrix) * sigma_matrix)
        chol_m = (chol(sigma_matrix))
        
        
        chol_m <- chol(sigma_matrix)
        
        b <- array(data = NA, c(2, M, N))
        
        
        # Generate all M sample tracks at once
        b[1, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
        b[2, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
        
        B[[i]] = b
        
        P[i, 1:M] = -mvnfast::dmvn(b[1, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE) -
          mvnfast::dmvn(b[2, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE)
      }
      
      
    } 
  }
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X",  "M",
                                  "covlist", "bilinearGradVec", "B", "P", 
                                  "times", "ID", "delta_max"), envir = environment())
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
    o = optim(par = c(0,0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf,-Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    o
    print(t)
    
    params[(k-1)*10 + d, ] = c(o$par, delta_max)
    print(d)
  }
  df25 = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], beta4 = params[,4], gammasq = params[,5],delta_max = params[,6])
  
  save(df, file = "sea_lion_deltamax_studyM=25.Rda")
  
}


#sea lion estimates on grid using M=50 bridges
M = 50
deltas = exp(seq(log(0.01), log(25), length.out = 30))
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))
for (k in 1:length(deltas)) {
  for (d in 1:10) {
    delta_max = deltas[k]
    
    #brownian bridge array
    B <- list()
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      if(ID[i] == ID[i+1]){
        N = ceiling((times[i+1] - times[i])/delta_max)
        #N = max(2, N)
        if(N != 1){
          delta = (times[i+1] - times[i])
          
          #brownian bridge covariance matrix
          sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
          sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
            t(lower.tri(sigma_matrix) * sigma_matrix)
          chol_m = (chol(sigma_matrix))
          
          
          chol_m <- chol(sigma_matrix)
          
          b <- array(data = NA, c(2, M, N))
          
          
          # Generate all M sample tracks at once
          b[1, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          b[2, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          
          B[[i]] = b
          
          P[i, 1:M] = -mvnfast::dmvn(b[1, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE) -
            mvnfast::dmvn(b[2, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE)
        }
        
        
      } 
    }
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X",  "M",
                                  "covlist", "bilinearGradVec", "B", "P", 
                                  "times", "ID", "delta_max"), envir = environment())
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
    o = optim(par = c(0,0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf,-Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    o
    print(t)
    
    params[(k-1)*10 + d, ] = c(o$par, delta_max)
    print(d)
  }
  df50 = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], beta4 = params[,4], gammasq = params[,5],delta_max = params[,6])
  
  save(df, file = "sea_lion_deltamax_studyM=50.Rda")
  
}


#sea lion estimates on grid using M=100 bridges
M = 100
deltas = exp(seq(log(0.01), log(25), length.out = 30))
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))
for (k in 1:length(deltas)) {
  for (d in 1:10) {
    delta_max = deltas[k]
    
    #brownian bridge array
    B <- list()
    #importance sampling wights
    P <- array(data = NA, c(nrow(X)-1,M))
    #generating bridges
    for (i in 1:(nrow(X)-1)) {
      if(ID[i] == ID[i+1]){
        N = ceiling((times[i+1] - times[i])/delta_max)
        #N = max(2, N)
        if(N != 1){
          delta = (times[i+1] - times[i])
          
          #brownian bridge covariance matrix
          sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
          sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
            t(lower.tri(sigma_matrix) * sigma_matrix)
          chol_m = (chol(sigma_matrix))
          
          
          chol_m <- chol(sigma_matrix)
          
          b <- array(data = NA, c(2, M, N))
          
          
          # Generate all M sample tracks at once
          b[1, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          b[2, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
          
          B[[i]] = b
          
          P[i, 1:M] = -mvnfast::dmvn(b[1, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE) -
            mvnfast::dmvn(b[2, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE, log = TRUE)
        }
        
        
      } 
    }
    
    #using paralellized and vectorized likelihood in optim
    cl <- makeCluster(ncores)
    
    clusterExport(cl, varlist = c("X",  "M",
                                  "covlist", "bilinearGradVec", "B", "P", 
                                  "times", "ID", "delta_max"), envir = environment())
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
    o = optim(par = c(0,0,0,0,1), fn = function(x) lik_grad(x, cl)$l, gr = function(x) lik_grad(x, cl)$g, method = "L-BFGS-B", lower = c(-Inf, -Inf,-Inf, -Inf, 0.0001))
    t = Sys.time()-t
    stopCluster(cl)
    o
    print(t)
    
    params[(k-1)*10 + d, ] = c(o$par, delta_max)
    print(d)
  }
  df100 = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], beta4 = params[,4], gammasq = params[,5],delta_max = params[,6])
  
  save(df, file = "sea_lion_deltamax_studyM=100.Rda")
  
}




# Plots

load("sea_lion_deltamax_studyM=25.Rda")
load("sea_lion_deltamax_studyM=50.Rda")
load("sea_lion_deltamax_studyM=100.Rda")


library(ggplot2)
library(dplyr)
#beta1
z <- qnorm(0.95)  # 90% interval
summarise_gauss <- function(df, label) {
  df %>%
    group_by(delta_max) %>%
    summarise(
      mu = median(beta1),
      sd = sd(beta1),
      .groups = "drop"
    ) %>%
    mutate(
      lo = mu - z * sd,
      hi = mu + z * sd,
      M = label
    )
}

sum25  <- summarise_gauss(df25,  "M=25")
sum50  <- summarise_gauss(df50,  "M=50")
sum100 <- summarise_gauss(df100, "M=100")

sum_all <- bind_rows(sum25, sum50, sum100)

ggplot() +
  geom_point(data = df25,  aes(delta_max, beta1, color = "M=25"),  alpha = 0.1) +
  geom_point(data = df50,  aes(delta_max, beta1, color = "M=50"),  alpha = 0.1) +
  geom_point(data = df100, aes(delta_max, beta1, color = "M=100"), alpha = 0.1) +
  #geom_ribbon(data = sum_all, aes(x = delta_max, ymin = lo, ymax = hi, fill = M, group = M), alpha = 0.15, colour = NA) +
  geom_line(data = sum_all, aes(delta_max, mu, color = M)) +
  #geom_line(data = sum_all, aes(delta_max, lo, color = M), linetype = "dashed") +
  #geom_line(data = sum_all, aes(delta_max, hi, color = M), linetype = "dashed") +
  
  scale_x_log10() +
  labs(x = expression(Delta[max]), y = expression(beta[1]), color = NULL) +
  theme_bw()




#beta2
z <- qnorm(0.95)  # 90% interval
summarise_gauss <- function(df, label) {
  df %>%
    group_by(delta_max) %>%
    summarise(
      mu = median(beta2),
      sd = sd(beta2),
      .groups = "drop"
    ) %>%
    mutate(
      lo = mu - z * sd,
      hi = mu + z * sd,
      M = label
    )
}

sum25  <- summarise_gauss(df25,  "M=25")
sum50  <- summarise_gauss(df50,  "M=50")
sum100 <- summarise_gauss(df100, "M=100")

sum_all <- bind_rows(sum25, sum50, sum100)

ggplot() +
  geom_point(data = df25,  aes(delta_max, beta2, color = "M=25"),  alpha = 0.1) +
  geom_point(data = df50,  aes(delta_max, beta2, color = "M=50"),  alpha = 0.1) +
  geom_point(data = df100, aes(delta_max, beta2, color = "M=100"), alpha = 0.1) +
  #geom_ribbon(data = sum_all, aes(x = delta_max, ymin = lo, ymax = hi, fill = M, group = M), alpha = 0.15, colour = NA) +
  geom_line(data = sum_all, aes(delta_max, mu, color = M)) +
  #geom_line(data = sum_all, aes(delta_max, lo, color = M), linetype = "dashed") +
  #geom_line(data = sum_all, aes(delta_max, hi, color = M), linetype = "dashed") +
  
  scale_x_log10() +
  labs(x = expression(Delta[max]), y = expression(beta[2]), color = NULL) +
  theme_bw()




#beta3
z <- qnorm(0.95)  # 90% interval
summarise_gauss <- function(df, label) {
  df %>%
    group_by(delta_max) %>%
    summarise(
      mu = median(beta3),
      sd = sd(beta3),
      .groups = "drop"
    ) %>%
    mutate(
      lo = mu - z * sd,
      hi = mu + z * sd,
      M = label
    )
}

sum25  <- summarise_gauss(df25,  "M=25")
sum50  <- summarise_gauss(df50,  "M=50")
sum100 <- summarise_gauss(df100, "M=100")

sum_all <- bind_rows(sum25, sum50, sum100)

ggplot() +
  geom_point(data = df25,  aes(delta_max, beta3, color = "M=25"),  alpha = 0.1) +
  geom_point(data = df50,  aes(delta_max, beta3, color = "M=50"),  alpha = 0.1) +
  geom_point(data = df100, aes(delta_max, beta3, color = "M=100"), alpha = 0.1) +
  #geom_ribbon(data = sum_all, aes(x = delta_max, ymin = lo, ymax = hi, fill = M, group = M), alpha = 0.15, colour = NA) +
  geom_line(data = sum_all, aes(delta_max, mu, color = M)) +
  #geom_line(data = sum_all, aes(delta_max, lo, color = M), linetype = "dashed") +
  #geom_line(data = sum_all, aes(delta_max, hi, color = M), linetype = "dashed") +
  
  scale_x_log10() +
  labs(x = expression(Delta[max]), y = expression(beta[3]), color = NULL) +
  theme_bw()




#beta4
z <- qnorm(0.95)  # 90% interval
summarise_gauss <- function(df, label) {
  df %>%
    group_by(delta_max) %>%
    summarise(
      mu = median(beta4),
      sd = sd(beta4),
      .groups = "drop"
    ) %>%
    mutate(
      lo = mu - z * sd,
      hi = mu + z * sd,
      M = label
    )
}

sum25  <- summarise_gauss(df25,  "M=25")
sum50  <- summarise_gauss(df50,  "M=50")
sum100 <- summarise_gauss(df100, "M=100")

sum_all <- bind_rows(sum25, sum50, sum100)

ggplot() +
  geom_point(data = df25,  aes(delta_max, beta4, color = "M=25"),  alpha = 0.1) +
  geom_point(data = df50,  aes(delta_max, beta4, color = "M=50"),  alpha = 0.1) +
  geom_point(data = df100, aes(delta_max, beta4, color = "M=100"), alpha = 0.1) +
  #geom_ribbon(data = sum_all, aes(x = delta_max, ymin = lo, ymax = hi, fill = M, group = M), alpha = 0.15, colour = NA) +
  geom_line(data = sum_all, aes(delta_max, mu, color = M)) +
  #geom_line(data = sum_all, aes(delta_max, lo, color = M), linetype = "dashed") +
  #geom_line(data = sum_all, aes(delta_max, hi, color = M), linetype = "dashed") +
  
  scale_x_log10() +
  labs(x = expression(Delta[max]), y = expression(beta[4]), color = NULL) +
  theme_bw()



#gammasq
z <- qnorm(0.95)  # 90% interval
summarise_gauss <- function(df, label) {
  df %>%
    group_by(delta_max) %>%
    summarise(
      mu = median(gammasq),
      sd = sd(gammasq),
      .groups = "drop"
    ) %>%
    mutate(
      lo = mu - z * sd,
      hi = mu + z * sd,
      M = label
    )
}

sum25  <- summarise_gauss(df25,  "M=25")
sum50  <- summarise_gauss(df50,  "M=50")
sum100 <- summarise_gauss(df100, "M=100")

sum_all <- bind_rows(sum25, sum50, sum100)

ggplot() +
  geom_point(data = df25,  aes(delta_max, gammasq, color = "M=25"),  alpha = 0.1) +
  geom_point(data = df50,  aes(delta_max, gammasq, color = "M=50"),  alpha = 0.1) +
  geom_point(data = df100, aes(delta_max, gammasq, color = "M=100"), alpha = 0.1) +
  #geom_ribbon(data = sum_all, aes(x = delta_max, ymin = lo, ymax = hi, fill = M, group = M), alpha = 0.15, colour = NA) +
  geom_line(data = sum_all, aes(delta_max, mu, color = M)) +
  #geom_line(data = sum_all, aes(delta_max, lo, color = M), linetype = "dashed") +
  #geom_line(data = sum_all, aes(delta_max, hi, color = M), linetype = "dashed") +
  
  scale_x_log10() +
  labs(x = expression(Delta[max]), y = expression(gamma^2), color = NULL) +
  theme_bw()


