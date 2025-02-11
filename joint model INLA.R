install.packages("inlabru")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 
install.packages("graph")

library(devtools)
githubinstall("HKprocess")
install_github("HristosTyr/version 0.1-1")

library(HKprocess)

install.packages("HKprocess.zip")

BiocManager::install("Rgraphviz")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HKprocess")

library(HKprocess)
library(Rgraphviz)

a#Setting things u"inlabru"#Setting things up
library(inlabru)
library(INLA)
library(ggplot2)
bru_safe_sp(force = TRUE)

data(gorillas_sf, package = "inlabru")


matern <- inla.spde2.pcmatern(gorillas_sf$mesh,
                              prior.range = c(0.1, 0.01),
                              prior.sigma = c(1, 0.01)
)

cmp <- ~
  Common(geometry, model = matern) +
  Difference(geometry, model = matern) +
  Intercept(1)

