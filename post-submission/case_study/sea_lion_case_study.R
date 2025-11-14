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
    N = max(2, N)
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



M = 50
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




deltas = c(5,1,0.5, 0.2, 0.1, 0.05, 0.01)
params = matrix(NA, ncol = 5, nrow = length(deltas))
for (k in 1:length(deltas)) {
  delta_max = deltas[k]
  
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
  
  params[k, ] = o$par
}

ggplot() +
  geom_line(aes(x = deltas, y = params[,5])) +
  labs(x = "delta_max") +
  theme_bw()

library(ggplot2)

# Evaluate covariate gradients at observed locations
gradarray <- bilinearGradArray(locs = xy, cov_list = covlist)


# Fit model
fit <- langevinUD(locs = xy, times = time, ID = ID, 
                  grad_array = gradarray)


fit$betaHat
fit$gamma2Hat



#M = 50, deltamax = 1
c(6.965460e-04,  4.829506e-04, -1.071735e-04,  2.696888e-05,  1.242916e+01)

#M = 100, delta_max = 0.1
c(1.058950e-03, -8.897097e-04, -1.848609e-04 , 6.769654e-05 , 1.200593e+01)

#M = 200, delta_max = 0.05
c(-3.578747e-03, -9.469007e-06, -9.559275e-05, -1.009843e-04 , 1.109573e+00)

#M = 50, delta_max = 1000
c(1.170703e-04 , 2.130799e-03, -2.425340e-05,  3.498373e-06,  1.234881e+01)

#michelot 2019 function
c(1.386441e-04,  9.205626e-02, -2.469842e-05,  3.566714e-06,  12.35642)








