install.packages("inlabru")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 
install.packages("graph")


#Setting things u"inlabru"#Setting things up
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

?Common


fml.major <- geometry ~ Intercept + Common + Difference / 2
fml.minor <- geometry ~ Intercept + Common - Difference / 2




lik_minor <- like("cp",
                  formula = fml.major,
                  data = gorillas_sf$nests[gorillas_sf$nests$group == "major", ],
                  samplers = gorillas_sf$boundary,
                  domain = list(geometry = gorillas_sf$mesh)
)
lik_major <- like("cp",
                  formula = fml.minor,
                  data = gorillas_sf$nests[gorillas_sf$nests$group == "minor", ],
                  samplers = gorillas_sf$boundary,
                  domain = list(geometry = gorillas_sf$mesh)
)





jfit <- bru(cmp, lik_major, lik_minor,
            options = list(
              control.inla = list(int.strategy = "eb",
                                  parallel.linesearch = TRUE),
              bru_max_iter = 1
            )
)


install.packages("patchwork")

library(patchwork)
pl.major <- ggplot() +
  gg(gorillas_sf$mesh,
     mask = gorillas_sf$boundary,
     col = exp(jfit$summary.random$Common$mean)
  )
pl.minor <- ggplot() +
  gg(gorillas_sf$mesh,
     mask = gorillas_sf$boundary,
     col = exp(jfit$summary.random$Difference$mean)
  )
(pl.major + scale_fill_continuous(trans = "log")) +
  (pl.minor + scale_fill_continuous(trans = "log")) &
  theme(legend.position = "right")

















