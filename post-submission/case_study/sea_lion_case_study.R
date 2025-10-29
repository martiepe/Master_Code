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

