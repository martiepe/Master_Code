# Simulation study
# prep workspace ---------------------------------------------------------- ####
library(here)
# custom functions
source(here("functions/utility_functions.R"))  # custom general perpose functions
sourceDir("functions")  # custom function to load all functions in folder
load_lib(mvnfast, parallel, terra, dplyr,
         ambient, Rcpp)  # custom function to install & load packages
Rcpp::sourceCpp("compute_lik_grad_full.cpp")
cpp_path <- here("compute_lik_grad_full.cpp")

# Define parameters up front ---------------------------------------------- ####
set.seed(123)

## track pars 
speed <- 5              # speed parameter for Langevin model
dt    <- 0.01           # temporal resolution of simulated tracks
beta  <- c(4, 2, -0.1)  # covariate coefficients

## default estimation pars
ncores <- 20     # number of cores used in parallel computations
thin   <- 100    # thinning
N      <- thin-1 # default nodes
M      <- 50     # default number of bridges
n_obs  <- 5000   # default number of observations
n_sim  <- 100    # number of simulations per simulation study

## covariate pars
res  <- 1  # resolution of covariates 
ncov <- 2  # number of covariates
ext  <- c(-1, 1, -1, 1)*250  # extent of study area
perlin_f <- 0.05  # Perlin noise frequency

# simulate covariates with Perlin noise ----------------------------------- ####
covlist <- list()
xgrid <- seq(ext[1], ext[2], by = res)
ygrid <- seq(ext[3], ext[4], by = res)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals <- 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = perlin_f)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}

# Include squared distance to centre of map as covariate
xgrid <- seq(ext[1], ext[2], by = res)
ygrid <- seq(ext[3], ext[4], by = res)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x = xgrid, y = ygrid,
                     z = matrix(dist2, length(xgrid), length(ygrid)))


# Sim 1: varying delta_t, fixed number of observations -------------------- ####
print("varying delta_t, fixed number of observations")
params <- matrix(NA, ncol = 6, nrow = 5*n_sim)
sim_var <- c(5, 10, 20, 50, 100)
for (ik in 1:n_sim) {
  for (jk in seq_along(sim_var)) {
    # set up simulation parameters
    beta_sim <- beta
    thin_sim <- sim_var[jk]
    dt_sim <- dt
    delta <- dt_sim*thin_sim
    N_sim <- thin_sim-1
    M_sim <- M
    n_obs_sim <- n_obs
    Tmax <- n_obs_sim*thin_sim*dt_sim
    
    # simulating track
    X <- simLMM(delta, speed, covlist, beta_sim, c(0,0), n_obs_sim)
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                      ncores = ncores, cpp_path = cpp_path)  
    # store results
    params[ik*5+jk-5, 1:4] <- out$par
    params[ik*5+jk-5, 5] <- delta
    params[ik*5+jk-5, 6] <- as.numeric(out$time, units = "secs")
  }
  
  df <- data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], 
                   gammasq = params[,4], dt = as.factor(params[,5]), 
                   time = params[,6])
  save(df, file = "post-submission/simulation_studies/varying_thin_estimates.Rda")
}

# Sim 2: varying delta_t, fixed maximum time ------------------------------ ####
print("varying delta_t, fixed maximum time")
params <- matrix(NA, ncol = 6, nrow = 5*n_sim)
sim_var <- c(5, 10, 20, 50, 100)
for (ik in 1:n_sim) {
  for (jk in seq_along(sim_var)) {
    # set up simulation parameters
    beta_sim <- beta
    thin_sim <- sim_var[jk]
    dt_sim <- dt
    delta <- dt_sim*thin_sim
    N_sim <- thin_sim-1
    M_sim <- M
    Tmax <- 500
    n_obs_sim <- Tmax/(dt_sim*thin_sim)
    
    # simulating track
    X <- simLMM(delta, speed, covlist, beta_sim, c(0,0), n_obs_sim)
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    # store results
    params[ik*5+jk-5, 1:4] <- out$par
    params[ik*5+jk-5, 5] <- delta
    params[ik*5+jk-5, 6] <- as.numeric(out$time, units = "secs")
  }
  # save output
  df <- data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3],
                   gammasq = params[,4], dt = as.factor(params[,5]), 
                   time = params[,6])
  save(df,file = "post-submission/simulation_studies/varying_thin_estimates_fixed_Tmax.Rda")
}

# Sim 3: varying number of bridges (M) ------------------------------------ ####
print("varying M")
params <- matrix(NA, ncol = 6, nrow = 5*n_sim)
sim_var <- c(5, 10, 50, 100, 200)
for (ik in 1:n_sim) {
  beta_sim <- beta
  thin_sim <- thin
  dt_sim <- dt
  delta <- dt*thin
  N <- 49
  n_obs_sim <- n_obs
  Tmax <- n_obs_sim*thin_sim*dt_sim
  
  # simulating track
  X <- simLMM(delta, speed, covlist, beta, c(0,0), n_obs_sim)
  
  for (jk in seq_along(sim_var)) {
    M_sim <- sim_var[jk]
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    
    params[ik*5+jk-5, 1:4] = out$par
    params[ik*5+jk-5, 5] = M_sim
    params[ik*5+jk-5,6] = as.numeric(out$time, units = "secs")
  }
  # save output
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3],
                  gammasq = params[,4], M = as.factor(params[, 5]),
                  time = params[, 6])
  save(df,file = "post-submission/simulation_studies/varying_M_estimates_stochastic likelihood.Rda")
}

# Sim 4: varying number of nodes (N) -------------------------------------- ####
print("varying N")
params <- matrix(NA, ncol = 6, nrow = 5*n_sim)
sim_var <- c(4,9,49,99)
for (ik in 1:n_sim) {
  for (jk in seq_along(sim_var)) {
    beta_sim <- beta
    thin_sim <- thin
    dt_sim <- dt
    delta <- dt_sim*thin_sim
    N <- sim_var[jk]
    M_sim <- M
    n_obs_sim <- n_obs
    Tmax <- n_obs_sim*thin_sim*dt_sim
    
    # simulating track
    X <- simLMM(delta, speed, covlist, beta_sim, c(0,0), n_obs_sim)
    
    # fit model
    out <- fit_langevin_bbis(X, covlist, delta, N = N_sim, M = M_sim,
                             ncores = ncores, cpp_path = cpp_path)  
    # store results
    params[ik*4+jk-4, 1:4] = out$par
    params[ik*4+jk-4, 5] = N_sim
    params[ik*4+jk-4,6] = as.numeric(out$time, units = "secs")
  }
  # save output
  df = data.frame(beta1 = params[,1], beta2 = params[,2], beta3 = params[,3], 
                  gammasq = params[,4], N = as.factor(params[,5]), 
                  time = params[,6])
  save(df, file = "post-submission/simulation_studies/varying_N_estimates.Rda")
}
    