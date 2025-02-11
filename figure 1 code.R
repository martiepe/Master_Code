rm(list = ls())
library(Rhabit)
library(abind)
load("scenario1_utils.RData")
# Loading functions and parameters ----------------------------------------

#change by to 0.01
simulation_times <- seq(0, by = 1,
                        length.out = 15000 + 1)# WARNING Long trajectories
# The data will be subsampled every 50 point for estimation

# Simulation --------------------------------------------------------------

# MAY TAKE SOME TIME

# Simulation of  independant data sets

# IMPORTANT: a directory simulated_data/scenario1 should be created



data_path <- "simulated_data/scenario1/"

dir.create("simulated_data/scenario1", recursive = T)

file.create("simAnalyticData_seed.txt")

n_trajectories <- 100 # Like in the article, could be changed

# Simulating and saving independant trajectories

lapply(1:n_trajectories, function(seed){
  set.seed(seed)# reproducibility
  complete_sample <- simLangevinMM(beta = beta_true, gamma2 = 1,
                                   times = simulation_times, loc0 = c(0, 0),
                                   grad_fun = gradient_list, keep_grad = T,
                                   silent = T)
  write.table(x = complete_sample, file = paste0(data_path, "simAnalyticData_seed", seed, ".txt"), sep = ";", row.names = F, col.names = T)
  })



# Imaginary covariates grid -----------------------------------------------

x_grid <- y_grid <- seq(-16, 16, length = 8) # Grid on which covariate are sampled
#should cover all data
covariates_list <- lapply(fun_cov_list, 
                          function(foo){
                            covar_sampling_loc <- expand.grid(x_grid, y_grid)
                            z = matrix(apply(covar_sampling_loc, 1, foo), 
                                       nrow = length(x_grid))
                            return(list(x = x_grid, y = y_grid, z = z))
                          }) 

# Estimation function ------------------------------------------------------

estim_on_sample <- function(seed){
  complete_sample <- read.table(paste0(data_path, "simAnalyticData_seed",
                                       seed, ".txt"),
                                header = T, sep = ";")
  firstpoint <- 1
  sel <- seq(firstpoint, by = 50, length.out = 301) #subsampling
  xy <- as.matrix(complete_sample[sel, c("x", "y")])
  times <- complete_sample[sel, "t"]
  
  # computations of true gradient
  true_grad_array <- do.call(what = function(...) abind(..., along = 3),
                             args = lapply(gradient_list, 
                                           function(foo){
                                             t(apply(xy, 1, foo))
                                           }))
  #Interpolated covArray
  interp_grad_array <- bilinearGradArray(locs = xy, cov_list = covariates_list)
  interp_grad_array[,,3] <- true_grad_array[,,3] # The distance gradient is always explicit
  fit_true_grad <- Rhabit::langevinUD(locs = xy, times = times, 
                                      grad_array = true_grad_array)
  fit_interp_grad <- Rhabit::langevinUD(locs = xy, times = times, 
                                        grad_array = interp_grad_array)
  # Defining output
  params <- rep(c(rownames(fit_true_grad$betaHatVariance), "gamma2Hat"), 2)
  vals <- c(fit_true_grad$betaHat, fit_true_grad$gamma2Hat,
            fit_interp_grad$betaHat, fit_interp_grad$gamma2Hat)
  meth <- rep(c("true", "interp"), c(4, 4))
  seeds <- rep(seed, length(meth))
  result <- data.frame(param = params, estimation = vals, gradient = meth,
                       seed = seeds)
  return(result)
}

# DO NOT RUN ---------------------------------------------------------------


# The mclapply works on Linux and MacOS. 
# In windows, should be replaced  by a simple lapply or a "for" loop

library(parallel)
scenario1_result <- do.call(rbind.data.frame,
                            lapply(1:n_trajectories, estim_on_sample))
write.table(scenario1_result, file = "scenario1_result.txt", row.names = F)


scenario1_result

scenario1_result$estimation[(0:99)*8 + 1]

###########################
## Plot acceptance rates ##
###########################
p <- qplot(x=alldt, y=rates, xlab="Time interval", ylab="Acceptance rate") + 
  geom_point() + geom_line() + scale_x_continuous(trans = "log", breaks=alldt) + 
  theme(axis.title = element_text(size=14), axis.text = element_text(size=13))

pdf("MALArates.pdf", width=6, height=4)
plot(p)
dev.off()