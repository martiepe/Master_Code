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

# import data ####
tracks <- read.csv(here("post-submission/case_study/SSLpreddat.csv"))
habitat <- rast(here("post-submission/case_study/aleut_habitat.grd"))

# fit model ####
x <- tracks[tracks$ID == tracks$ID,]

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
tracks <- read.csv("C:/Users/marti/Desktop/Master_Code/post-submission/case_study/SSLpreddat.csv")
ID <- as.integer(tracks$ID)
time <- as.POSIXct(tracks$time)
time <- as.numeric(time)
time <- (time-min(time))/3600
xy <- matrix(c(tracks$x, tracks$y)/1000, ncol=2) # convert to km




#####################
## Load covariates ##
#####################
hbfull <- brick("C:/Users/marti/Desktop/Master_Code/post-submission/case_study/aleut_habitat.grd", values=TRUE)
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

# Crop covariates to area of interest
#border <- 30
#lim <- c(min(xy[,1])-border,max(xy[,1])+border,min(xy[,2])-border,max(xy[,2])+border)
#covlist0 <- lapply(covlist0, crop, y=extent(lim))

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
gradarray <- bilinearGradArray(locs = xy, cov_list = covlist)


############## THIRD ATTEMPT ######################
#simulating track
#X = simLMM(delta,speed,covlist,beta,c(0,0),n_obs)
#ID = c(rep(1, 2500), rep(2, 2500))
#times = delta*(1:n_obs)

##generating perlin covariate
#library(ambient)
#lim <- c(-1, 1, -1, 1)*250
#resol <- 1
#ncov <- 3
#covlist <- list()
#xgrid <- seq(lim[1], lim[2], by = resol)
#ygrid <- seq(lim[3], lim[4], by = resol)
#coords <- as.matrix(expand.grid(xgrid, ygrid))
#for(i in 1:ncov) {
#  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
#  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
#}
## Include squared distance to centre of map as covariate
#xgrid <- seq(lim[1], lim[2], by=resol)
#ygrid <- seq(lim[3], lim[4], by=resol)
#xygrid <- expand.grid(xgrid,ygrid)
#dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
#covlist[[ncov+1]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))
#
#
#
#thin = 100
#dt = 0.01
#delta = dt*thin
#N = thin-1
#M = 50
#n_obs = 5000
#Tmax = n_obs*thin*dt
#
#
##simulating track
#X = simLMM(delta,speed,covlist,c(4,2,1,-0.1),c(0,0),n_obs)
#ID = c(rep(1, 2500), rep(2, 2500))
#times = delta*(1:n_obs)

ncores = 10
X = xy 
times = time 



lik_grad <- function(par, cl){
  # par = [beta_1, ..., beta_p, s]
  p <- length(par) - 1L
  if (p < 1L) stop("par must be [betas..., s] with at least one covariate.")
  
  gamma <- sqrt(par[p + 1L])
  
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
    
    x_samples <- sweep(B[[i]][1, 1:M, 1:N] * gamma, 2, mu_x, "+")
    y_samples <- sweep(B[[i]][2, 1:M, 1:N] * gamma, 2, mu_y, "+")
    
    l_k <- P[i, 1:M]  + log(gamma)*(2 * N)
    
    full_x <- cbind(X[i, 1], x_samples, X[i + 1, 1])
    full_y <- cbind(X[i, 2], y_samples, X[i + 1, 2])
    
    # C++ must return a (1 + p + 1) x M matrix now
    # res_mat from compute_log_lik_grad_full_cpp:
    # row 1      : log w_j
    # rows 2..(p+1): score components for beta_j
    # row p+2      : score for diffusion parameter s
    
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
  
  results_list <- parLapply(cl, 1:(nrow(X) - 1L), compute)
  #results_list <- lapply(1:(nrow(X) - 1L), compute)
  res_mat <- do.call(rbind, results_list)  # (n_seg) x (p+2)
  
  # fold
  l_sum      <- sum(res_mat[, 1], na.rm = TRUE)
  g_betasum  <- colSums(res_mat[, 2:(p + 1L), drop = FALSE], na.rm = TRUE)
  g_ssum     <- sum(res_mat[, p + 2L], na.rm = TRUE)
  
  if (!is.finite(l_sum)) {
    return(list(l = 1e10, g = rep(0, p + 1L)))
  }
  #print(c(par, -l_sum, g_betasum, g_ssum))
  list(l = -l_sum, g = c(g_betasum, g_ssum))
}



M = 500
delta_max = 0.01
#brownian bridge array
B <- list()
#importance sampling wights
P <- array(data = NA, c(nrow(X)-1,M))
#generating bridges
for (i in 1:(nrow(X)-1)) {
  if(ID[i] == ID[i+1]){
    N = ceiling((times[i+1] - times[i])/delta_max)
    N = max(2, N)
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


M = 100
deltas = exp(seq(log(0.01), log(25), length.out = 30))
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))
for (k in 30:length(deltas)) {
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
  
}

df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], beta4 = params[,4], gammasq = params[,5],delta_max = params[,6])
save(df, file = "post-submission/case_study/sea_lion_deltamax_studyM=100.Rda")





M = 500
delta_max = 0.05
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))

for (d in 1:10) {
  
  
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
  
  params[d, ] = c(o$par, delta_max)
  print(d)
}
  
(params[1:10, 1])


ggplot() +
  geom_point(data = df,
             mapping = aes(x = delta_max, y = beta4),
             alpha = 0.5,
             color = "#F8766D") +
  geom_point(aes(x = rep(0.05, 10), y = params[1:10, 4])) +
  labs(x = "delta_max") +
  scale_x_log10() +
  theme_bw()

















# Varying number of bridges, recreating bridges within likelihood computation
library(Rcpp)
cpp_path <- here("compute_lik_grad_full.cpp")
sourceCpp(cpp_path)

n_trans <- nrow(X) - 1L

Var_max <- 1e-3   # per-transition Var(log p̂_i) target
M_max   <- 2000L  # safety cap

# initial number of bridges per transition (start at 1)
M_vec <- rep(2L, n_trans)

lik_grad <- function(par, cl) {
  # par = [beta_1, ..., beta_p, s]
  p <- length(par) - 1L
  if (p < 1L) stop("par must be [betas..., s] with at least one covariate.")
  
  gamma <- sqrt(par[p + 1L])
  
  compute <- function(i) {
    # default output if we skip this transition
    zero_contrib <- rep(0, p + 2L)
    
    # cross-ID transitions: no contribution, keep M as is
    if (ID[i] != ID[i + 1]) {
      return(list(contrib = zero_contrib, M_used = 0L))
    }
    
    N     <- ceiling((times[i + 1] - times[i]) / delta_max)
    delta <- (times[i + 1] - times[i])
    
    ## ----------------- N = 1: Maruyama shortcut (unchanged) -----------------
    if (N == 1L) {
      y <- c(X[i + 1, 1] - X[i, 1],
             X[i + 1, 2] - X[i, 2])
      
      G_arr <- bilinearGradVec(matrix(X[i, ], nrow = 1), covlist)
      G_i   <- t(drop(G_arr[, 1, ]))  # 2 x p
      
      if (!all(is.finite(G_i))) {
        return(list(contrib = zero_contrib, M_used = 0L))
      }
      
      beta <- par[1:p]
      s    <- gamma^2
      
      mu <- as.vector((delta * s / 2) * (G_i %*% beta))  # length 2
      e  <- y - mu
      
      loglike_i <- -log(2 * pi) - log(delta * s) - sum(e * e) / (2 * delta * s)
      
      g_beta <- 0.5 * as.numeric(crossprod(G_i, e))  # p-vector
      g_s    <- -1 / s +
        (sum(e * e)) / (2 * delta * s^2) +
        as.numeric(t(e) %*% (G_i %*% beta)) / (2 * s)
      
      if (!is.finite(loglike_i) || !all(is.finite(g_beta)) || !is.finite(g_s)) {
        return(list(contrib = zero_contrib, M_used = 0L))
      }
      
      contrib <- c(loglike_i, -g_beta, -g_s)
      return(list(contrib = contrib, M_used = 0L))  # M_used irrelevant for N=1
    }
    
    ## ----------------- N > 1: bridges + adaptive M_i -----------------
    
    # Brownian bridge means for interior points 1..N
    mu_x <- X[i, 1] + (1:N) * (X[i + 1, 1] - X[i, 1]) / (N + 1)
    mu_y <- X[i, 2] + (1:N) * (X[i + 1, 2] - X[i, 2]) / (N + 1)
    
    # Brownian bridge covariance
    t_idx        <- 1:N
    sigma_matrix <- delta * outer(1 - t_idx / (N + 1), t_idx / (N + 1))
    sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
      t(lower.tri(sigma_matrix) * sigma_matrix)
    chol_m <- chol(sigma_matrix)
    
    # starting M for this segment (from global vector)
    M0 <- M_vec[i]
    if (is.na(M0) || M0 < 1L) M0 <- 1L
    if (M0 > M_max) M0 <- M_max
    
    # running aggregators in log-space
    A      <- -Inf          # max log w
    T1     <- 0.0           # Σ exp(logw_k - A)
    T2     <- 0.0           # Σ exp(2(logw_k - A))
    S_beta <- numeric(p)    # Σ score_beta_k * exp(logw_k - A)
    S_s    <- 0.0           # Σ score_s_k * exp(logw_k - A)
    
    loglike_i_final <- NA_real_
    g_beta_final    <- rep(NA_real_, p)
    g_s_final       <- NA_real_
    
    # set seed once; we will draw M0 bridges in batch, then continue the RNG
    set.seed(i+1)
    
    ## ---------- Batch: first M0 bridges ----------
    # b1, b2: M0 x N
    b1 <- mvnfast::rmvn(M0, rep(0, N), sigma = chol_m, isChol = TRUE)
    b2 <- mvnfast::rmvn(M0, rep(0, N), sigma = chol_m, isChol = TRUE)
    
    x_samples <- sweep(b1 * gamma, 2L, mu_x, "+")  # M0 x N
    y_samples <- sweep(b2 * gamma, 2L, mu_y, "+")  # M0 x N
    
    log_phi1 <- mvnfast::dmvn(b1, rep(0, N), sigma = chol_m,
                              isChol = TRUE, log = TRUE)
    log_phi2 <- mvnfast::dmvn(b2, rep(0, N), sigma = chol_m,
                              isChol = TRUE, log = TRUE)
    log_L_batch <- -log_phi1 - log_phi2 + (2 * N) * log(gamma)  # length M0
    
    full_x_batch <- cbind(rep(X[i, 1], M0), x_samples, rep(X[i + 1, 1], M0))
    full_y_batch <- cbind(rep(X[i, 2], M0), y_samples, rep(X[i + 1, 2], M0))
    
    res_batch <- compute_log_lik_grad_full_cpp(
      full_x_batch, full_y_batch,
      log_L_k = log_L_batch,
      X       = X,
      i_index = i,
      par     = par,
      delta   = delta,
      N       = N,
      covlist = covlist
    )
    
    logw_batch <- res_batch[1, ]                     # length M0
    score_b_mat <- res_batch[2:(p + 1L), , drop = FALSE]  # p x M0
    score_s_vec <- res_batch[p + 2L, ]              # length M0
    
    if (any(!is.finite(logw_batch)) ||
        any(!is.finite(score_b_mat)) ||
        any(!is.finite(score_s_vec))) {
      return(list(contrib = zero_contrib, M_used = M0))
    }
    
    # initialise aggregators from batch
    A  <- max(logw_batch)
    z  <- exp(logw_batch - A)          # length M0
    T1 <- sum(z)
    T2 <- sum(z * z)
    S_beta <- score_b_mat %*% z        # p-vector
    S_s    <- sum(score_s_vec * z)
    
    M_i <- M0
    
    # compute loglike and grad from current aggregates
    loglike_i <- A + log(T1) - log(M_i)
    g_beta    <- -(as.numeric(S_beta) / T1) / 2
    g_s       <- -(S_s / T1)
    
    loglike_i_final <- loglike_i
    g_beta_final    <- g_beta
    g_s_final       <- g_s
    
    # delta-method variance after batch
    if (M_i >= 2L) {
      cv2        <- M_i * (T2 / (T1 * T1)) - 1.0
      var_loghat <- cv2 / M_i
      if (is.finite(var_loghat) && var_loghat <= Var_max) {
        contrib <- c(loglike_i_final, g_beta_final, g_s_final)
        return(list(contrib = contrib, M_used = M_i))
      }
    }
    
    ## ---------- Incremental: add step_M bridges at a time ----------
    step_M = 50
    repeat {
      if (M_i >= M_max) {
        warning(sprintf("Segment %d reached M_max = %d", i, M_max))
        break
      }
      
      # how many new bridges this round (cap at M_max)
      add_M <- min(step_M, M_max - M_i)
      if (add_M <= 0L) break
      
      # simulate add_M new bridges; RNG continues from where batch left off
      b1_new <- mvnfast::rmvn(add_M, rep(0, N), sigma = chol_m, isChol = TRUE)  # add_M x N
      b2_new <- mvnfast::rmvn(add_M, rep(0, N), sigma = chol_m, isChol = TRUE)
      
      x_samples_new <- sweep(b1_new * gamma, 2L, mu_x, "+")    # add_M x N
      y_samples_new <- sweep(b2_new * gamma, 2L, mu_y, "+")
      
      log_phi1_new <- mvnfast::dmvn(b1_new, rep(0, N), sigma = chol_m,
                                    isChol = TRUE, log = TRUE)
      log_phi2_new <- mvnfast::dmvn(b2_new, rep(0, N), sigma = chol_m,
                                    isChol = TRUE, log = TRUE)
      log_L_new <- -log_phi1_new - log_phi2_new + (2 * N) * log(gamma)  # length add_M
      
      full_x_new <- cbind(rep(X[i, 1], add_M), x_samples_new, rep(X[i + 1, 1], add_M))
      full_y_new <- cbind(rep(X[i, 2], add_M), y_samples_new, rep(X[i + 1, 2], add_M))
      
      res_new <- compute_log_lik_grad_full_cpp(
        full_x_new, full_y_new,
        log_L_k = log_L_new,
        X       = X,
        i_index = i,
        par     = par,
        delta   = delta,
        N       = N,
        covlist = covlist
      )
      
      logw_new     <- res_new[1, ]                         # length add_M
      score_b_new  <- res_new[2:(p + 1L), , drop = FALSE]  # p x add_M
      score_s_new  <- res_new[p + 2L, ]                    # length add_M
      
      if (any(!is.finite(logw_new)) ||
          any(!is.finite(score_b_new)) ||
          any(!is.finite(score_s_new))) {
        break
      }
      
      # update aggregators one new bridge at a time (no extra C++ calls)
      for (k in 1:add_M) {
        logw_k    <- logw_new[k]
        score_b_k <- score_b_new[, k]
        score_s_k <- score_s_new[k]
        
        if (logw_k <= A) {
          z  <- exp(logw_k - A)
          T1 <- T1 + z
          T2 <- T2 + z * z
          S_beta <- S_beta + score_b_k * z
          S_s    <- S_s + score_s_k * z
        } else {
          z_scale <- exp(A - logw_k)
          T1      <- T1 * z_scale + 1.0
          T2      <- T2 * (z_scale * z_scale) + 1.0
          S_beta  <- S_beta * z_scale + score_b_k
          S_s     <- S_s * z_scale + score_s_k
          A       <- logw_k
        }
      }
      
      M_i <- M_i + add_M
      
      # recompute loglike and gradient
      loglike_i <- A + log(T1) - log(M_i)
      g_beta    <- -(as.numeric(S_beta) / T1) / 2
      g_s       <- -(S_s / T1)
      
      loglike_i_final <- loglike_i
      g_beta_final    <- g_beta
      g_s_final       <- g_s
      
      # MC variance with M_i bridges
      if (M_i >= 2L) {
        cv2        <- M_i * (T2 / (T1 * T1)) - 1.0
        var_loghat <- cv2 / M_i
        if (is.finite(var_loghat) && var_loghat <= Var_max) {
          break
        }
      }
    }
    
    
    contrib <- c(loglike_i_final, g_beta_final, g_s_final)
    list(contrib = contrib, M_used = M_i)
  }
  
  # parallel over segments
  results_list <- parLapply(cl, 1:(nrow(X) - 1L), compute)
  
  # extract contributions and updated M values
  contrib_mat <- do.call(rbind, lapply(results_list, `[[`, "contrib"))  # (n_seg) x (p+2)
  M_used_new  <- vapply(results_list, function(x) as.integer(x$M_used), integer(1))
  
  
  # update global M_vec (non-decreasing)
  M_vec <<- pmax(M_vec, M_used_new)
  
  # fold
  l_sum     <- sum(contrib_mat[, 1], na.rm = TRUE)
  g_betasum <- colSums(contrib_mat[, 2:(p + 1L), drop = FALSE], na.rm = TRUE)
  g_ssum    <- sum(contrib_mat[, p + 2L], na.rm = TRUE)
  
  if (!is.finite(l_sum)) {
    return(list(l = 1e10, g = rep(0, p + 1L)))
  }
  print(par)
  list(l = -l_sum, g = c(g_betasum, g_ssum))
}


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

M_vec <- rep(2L, n_trans)


Var_max <- 1/5342
M_max     <- Inf      # already used above
delta_max <- 0.05
ncores    <- 10         # or whatever

# parallel setup
cl <- makeCluster(ncores)

clusterExport(cl, varlist = c(
  "X", "ID", "times",
  "covlist", "bilinearGradVec",
  "delta_max", "Var_max", "M_max",
  "M_vec"  # M_vec is captured in compute, but exporting doesn't hurt
), envir = environment())

clusterEvalQ(cl, library(mvnfast))

cpp_path <- here("compute_lik_grad_full.cpp")

clusterExport(cl, varlist = "cpp_path")

clusterEvalQ(cl, {
  library(Rcpp)
  sourceCpp(cpp_path)
})

# build cached wrappers
wrapped <- make_lik_grad_cached(cl)

t <- Sys.time()
o <- optim(
  par   = c(7.5e-05, -0.000175, 0.001, 0.0009, 12.5),
  fn    = wrapped$fn,
  gr    = wrapped$gr,
  method = "L-BFGS-B",
  lower  = c(-Inf, -Inf, -Inf, -Inf, 0.0001),
  control = list(factr = 1e2, pgtol = 1e-8)
)
t <- Sys.time() - t

stopCluster(cl)

o
print(t)
hist(M_vec)

















load("post-submission/case_study/sea_lion_deltamax_studyM=100.Rda")


M = 500
delta_max = 0.05
params = matrix(NA, ncol = 6, nrow = 10*length(deltas))

for (d in 1:10) {
  
  
  
  params[d, ] = c(o$par, delta_max)
  print(d)
  
  
  ggplot() +
  geom_point(data = df,
             mapping = aes(x = delta_max, y = beta4),
             alpha = 0.5,
             color = "#F8766D") +
  geom_point(aes(x = rep(0.05, 10), y = params[1:10, 4])) +
  labs(x = "delta_max") +
  scale_x_log10() +
  theme_bw()
  
}








# [1]  9.279162e-04 -4.472065e-04 -1.811526e-04  6.279265e-05  1.257023e+01

# [1]  9.812248e-04 -6.172214e-04 -1.928540e-04  7.664456e-05  1.253654e+01





















