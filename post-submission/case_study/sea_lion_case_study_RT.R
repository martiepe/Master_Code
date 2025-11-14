# Case study
# prep workspace ####
library(here)
library(mvnfast)
library(parallel)
library(terra)
library(dplyr)
# custom functions
source(here("functions/utility_functions.R"))
sourceDir("functions")

# import data ####
hbfull <- rast(here("post-submission/case_study/aleut_habitat.grd"))
tracks <- read.csv(here("post-submission/case_study/SSLpreddat.csv")) %>% 
  mutate(time = as.POSIXct(time)) %>% 
  vect(geom = c("x","y"), crs = crs(hbfull))

# Change the resolution and extent from m to km
crs <- gsub("units=m","units=km",crs(hbfull, T))
r <- rast(nrows = nrow(hbfull), ncols = ncol(hbfull),
          ext = as.vector(ext(hbfull))/1000,
          crs = crs) # define template raster

# standardise projections
hbfull <- project(hbfull, r)     # transform raster
tracks <- project(tracks, r) %>% # transform tracks
  as.data.frame(geom = "XY")

#######################
## fit Langevin BBIS ##
#######################
ncores <- 10
M <- 50
dt_max <- 1
out <- fit_langevin_bbis(tracks, hbfull, 
                         M = M,
                         dt_max = dt_max, 
                         dt_units = "hours",
                         ncores = 10, 
                         fixed_sampling = FALSE) 

# explore effect of dmax
deltas = c(5,1,0.5, 0.2, 0.1, 0.05, 0.01)
params = matrix(NA, ncol = 5, nrow = length(deltas))
for (k in 1:length(deltas)) {
  # fit
  out <- fit_langevin_bbis(tracks, hbfull, 
                           M = M,
                           dt_max = deltas[k], 
                           dt_units = "hours",
                           ncores = 10, 
                           fixed_sampling = FALSE) 
  # store outputs
  params[k, ] = o$par
}

ggplot() +
  geom_line(aes(x = deltas, y = params[,5])) +
  labs(x = "delta_max") +
  theme_bw()

library(ggplot2)
