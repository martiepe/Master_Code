# generate conceptual plot of BBIS #
# load packages 
library(here)
library(mvnfast)
library(parallel)
library(ambient)  # perlin noise
library(Rcpp)
Rcpp::sourceCpp("compute_lik_grad_full.cpp")
library(here)
library(terra)
library(dplyr)
library(ggplot2)
library(ggbrace)
library(patchwork)

# define custom functions
{
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
  
  # thin highresolution track
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
}


# simulate covariate raster layers
seed <- runif(1,1,1e6); print(seed)
set.seed(265240.4)
# simulate data


#generating perlin covariates
lim <- c(-1, 1, -1, 1)*100
resol <- 1
ncov <- 2
covlist <- list()
xgrid <- seq(lim[1], lim[2], by = resol)
ygrid <- seq(lim[3], lim[4], by = resol)
coords <- as.matrix(expand.grid(xgrid, ygrid))
for(i in 1:ncov) {
  vals = 3*noise_perlin(c(length(xgrid), length(ygrid)), frequency = 0.05)
  covlist[[i]] = list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid)))
}
# Include squared distance to centre of map as covariate
xgrid <- seq(lim[1], lim[2], by=resol)
ygrid <- seq(lim[3], lim[4], by=resol)
xygrid <- expand.grid(xgrid,ygrid)
dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/(100)
covlist[[3]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))


# define parameters
##  simulation parameters 
speed <- 5  # speed parameter for Langevin model
beta <- c(-10,0,0)  # covariate selection coefficients
# beta <- c(4,2,-0.1)  # covariate selection coefficients
n_tracks <- 1  # number of tracks
# bridge parameters
thin <- 20 # c(5, 10, 20, 50, 100)
dt <- 0.01
M <- 50
n_obs <- 2

# generate data
delta = dt*thin
N = thin-1
Tmax = n_obs*thin*dt

# simulating track
x <- sumLTrack(delta,speed,covlist,beta,c(0,0),n_obs)  # Simulate high-resolution track
X <- thinTrack(x, delta)

# brownian bridge covariance matrix
sigma_matrix <- delta * outer(1 - 1:N/(N+1), 1:N/(N+1))
sigma_matrix <- lower.tri(sigma_matrix, TRUE) * sigma_matrix +
  t(lower.tri(sigma_matrix) * sigma_matrix)
chol_m = (chol(sigma_matrix))

# brownian bridge endpoints
mu_x_all <- rep(X[1:(nrow(X)-1), 1], each = N) + 1:N * rep((X[2:nrow(X), 1] - X[1:(nrow(X)-1), 1]), each = N) / (N+1)
mu_y_all <- rep(X[1:(nrow(X)-1), 2], each = N) + 1:N * rep((X[2:nrow(X), 2] - X[1:(nrow(X)-1), 2]), each = N) / (N+1)
# brownian bridge array
B <- array(data = NA, c(2, nrow(X)-1, M, N))
# importance sampling wights
P <- array(data = NA, c(nrow(X)-1,M))
# generating bridges
for (i in 1:(nrow(X)-1)) {
  mu_x <- mu_x_all[((i - 1) * N + 1):(i * N)]
  mu_y <- mu_y_all[((i - 1) * N + 1):(i * N)]
  
  # Generate all M sample tracks at once
  B[1, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
  B[2, i, 1:M, 1:N] <- mvnfast::rmvn(M, rep(0,N), sigma = chol_m, isChol = TRUE)
  
  P[i, 1:M] = 1/(mvnfast::dmvn(B[1, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE) * 
                   mvnfast::dmvn(B[2, i, 1:M, 1:N], rep(0,N), sigma = chol_m, isChol = TRUE))
}

# generate plot ####
# prep plot data

# define step to visualize
time_step <- 1
obs_start <- X[time_step, ]
obs_end <- X[time_step + 1, ]

# Create the latent track between these two points
# The latent track is sampled at 0.01 intervals
latent_start_idx <- (time_step - 1) * thin + 1
latent_end_idx <- time_step * thin + 1

# Reconstruct the latent path (you'll need to save this from simLMM)
# For now, we'll extract it from the bridge midpoints plus endpoints
latent_x <- x[latent_start_idx:latent_end_idx,1]*-1
latent_y <- x[latent_start_idx:latent_end_idx,2]*-1.15
latent_track <- data.frame(x = latent_x, y = latent_y, 
                           t = seq(0, delta, length.out = N + 2),
                           type = "Latent Track",
                           linewidth = 0.75, alpha = 1)

# Select a subset of bridges to plot (e.g., 8 out of M bridges)
n_bridges_plot <- 4
bridge_indices <- seq(1, M, length.out = n_bridges_plot)

# Prepare bridge data
bridge_data <- data.frame()
for (m in bridge_indices) {
  bridge_x <- c(obs_start[1], 
                B[1, time_step, m, ] + mu_x_all[((time_step - 1) * N + 1):(time_step * N)], 
                obs_end[1])
  bridge_y <- c(obs_start[2], 
                B[2, time_step, m, ] + mu_y_all[((time_step - 1) * N + 1):(time_step * N)], 
                obs_end[2])
  
  bridge_df <- data.frame(
    x = bridge_x*-1,
    y = bridge_y*-1.15,
    t = seq(0, delta, length.out = N + 2),
    bridge_id = as.factor(m),
    type = "Brownian Bridge",
    linewidth = 0.5, alpha = 0.4
  )
  bridge_data <- rbind(bridge_data, bridge_df)
}

# Observed points
obs_points <- data.frame(
  x = c(obs_start[1], obs_end[1])*-1,
  y = c(obs_start[2], obs_end[2])*-1.15,
  t = c(0, delta),
  label = c("italic(t)[0]", "italic(t)[1]"),
  type = "Observed",
  linewidth = 1,
  alpha = 1
)
obs_text <- obs_points
obs_text$y <- obs_text$y +c(-0.2, 0.25)

all_tracks <- bind_rows(obs_points, bridge_data, latent_track)
# Create the plot
## define colors
col_obs <- "#5A78C1"
col_latent <- "#CD7338"
col_bridge <- "#9B9B9B"
## define common theme
theme <- theme_minimal() + 
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.margin = margin(t = 0, unit = "pt"),
        legend.box.margin = margin(t = -10, unit = "pt"),
        panel.border = element_rect(fill = NA, color = "gray70"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# define y-axis limits
y_diff <- abs(diff(range(all_tracks$y)))
y_top <- max(all_tracks$y)+y_diff*0.3
y_bot <- min(all_tracks$y)+y_diff*0.1
y_mid <- y_bot + (y_top-y_bot)/2
y_ran <- y_top-y_bot
lab_pos <- data.frame(par = c("h","dt","M","N"),
                      pos = c(y_bot+y_ran*1,
                              y_bot+y_ran*0.9,
                              y_bot+y_ran*0.75,
                              y_bot+y_ran*0.6))
cur_dif <- 0.1
#
# space-space plot ####
p1 <-  ggplot(obs_points, aes(x = x, y = y, col = type)) +
  # paths
  geom_path(linewidth = 1) +
  geom_path(data = bridge_data, 
            aes(group = bridge_id),
            alpha = 0.4, linewidth = 0.5) +
  geom_path(data = latent_track, linewidth = 0.75) +
  # points
  geom_point(data = bridge_data, size = 0.5) +
  geom_point(data = latent_track,  size = 1) +
  geom_point(size = 2, shape = 21, fill = col_obs) +
  labs(x = "X", y = "Y") +
  coord_fixed()  + 
  theme  + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank())

# space-time plot ####
p2 <- ggplot(obs_points, aes(x = t, y = y, col = type)) +
  # paths
  geom_path(data = obs_points, linewidth = 1) +
  geom_path(data = bridge_data, aes(group = bridge_id),
            alpha = 0.4, linewidth = 0.5) +
  geom_path(data = latent_track, linewidth = 0.75) +
  # points
  geom_point(data = bridge_data, size = 0.5) +
  geom_point(data = latent_track, size = 1) +
  geom_point(data = obs_points, size = 2, shape = 21, fill = col_obs) +
  # Curly bracket for h
  annotate("text", x = dt/2-0.0025, 
           y = lab_pos[1,2],#y_top-y_diff*0.00+0.2, 
           hjust = 0, vjust = 0,
           label = paste0("italic(h) == ", dt), parse = TRUE,
           color = col_latent, size = 4) +
  stat_brace(aes(x = c(0, dt), 
                 y = rep(lab_pos[1,2]-cur_dif),2), #rep(y_top-y_diff*0.00, 2)), 
             width = y_ran/50, col = col_latent, linewidth = 0.5) +
  # Curly bracket for delta_t
  annotate("text", x = delta/2-0.0025, 
           y = lab_pos[2,2], #y_top-y_diff*0.2 + 0.2, 
           hjust = 0, vjust = 0,
           label = paste0("Delta[italic(t)] == ", delta), parse = TRUE,
           color = col_obs, size = 4) +
  stat_brace(aes(x = c(0, 0.2), 
                 y = rep(lab_pos[2,2]-cur_dif),2),#rep(y_top-y_diff*0.2, 2)), 
             width = y_ran/20, col = col_obs, linewidth = 0.5) +
  # M annotation with line
  geom_path(data = data.frame(x = c(0.5,1.25,2,2.75,3.5)/4*delta/10 + max(obs_points$t) * 0.65,
                              y = c(0.3,0,1,-0.75,0)/4*y_ran/20 + lab_pos[3,2]),
            aes(x = x, y = y),
            col = col_bridge, linewidth = 0.5, alpha = 0.5) + 
  # annotate("segment", 
  #          x = max(obs_points$t) * 0.65, 
  #          xend = max(obs_points$t) * 0.75,
  #          y = rep(lab_pos[3,2],2),#rep(y_top-y_diff*0.25, 2),
  #          yend =  rep(lab_pos[3,2],2),#rep(y_top-y_diff*0.25, 2),
  #          color = col_bridge, linewidth = 0.5) +
  annotate("text", 
           x = max(obs_points$t) * 0.8, 
           y = lab_pos[3,2], # rep(y_top-y_diff*0.25, 2),
           label = paste0("italic(M) == ", n_bridges_plot), parse = TRUE,
           hjust = 0, size = 4) +
  # N annotation with point
  annotate("point", 
           x = max(obs_points$t) * 0.7, 
           y = rep(lab_pos[4,2],,2), #rep(y_top-y_diff*0.35, 2),
           color = col_bridge, size = 0.5) +
  annotate("text", 
           x = max(obs_points$t) * 0.8, 
           y =  lab_pos[4,2], #rep(y_top-y_diff*0.35, 2),
           label = paste0("italic(N) == ", N), parse = TRUE,
           hjust = 0, size = 4) +
  # text
  geom_text(data = obs_text, aes(label = label), col = col_obs,
            size = 5, parse = TRUE) + 
  theme +
  scale_x_continuous(
    breaks = seq(min(bridge_data$t),
                 max(bridge_data$t),
                 by = 0.05),
    minor_breaks = seq(min(bridge_data$t),
                       max(bridge_data$t),
                       by = 0.01)
  ) +
  coord_cartesian(clip = "off") + 
  labs(x = "Time", y = "Y") + 
  theme(axis.text.y = element_blank())

# Combine plots with shared legend at bottom center  ####
combined_plot <- (p1 | p2) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  scale_y_continuous(limits = c(min(all_tracks$y), y_top),
                     expand = expansion(mult = c(0.1, 0.1))) &
  scale_color_manual(
    name = "",
    values = c("Observed" = col_obs, "Latent Track" = col_latent, "Brownian Bridge" = col_bridge),
    breaks = c("Observed", "Latent Track", "Brownian Bridge")
  ) &
  theme(legend.position = "bottom",
        legend.margin = margin(t = 0, b = 0, unit = "pt"),
        legend.box.margin = margin(t = -10, b = -5, unit = "pt"),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(hjust = 0, vjust = 0),
        plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))
# save plot
ggsave(
  here("figures/fig1_concept.pdf"),
  combined_plot, 
  device = cairo_pdf, 
  width = 7, 
  height = 3.5
); browseURL( here("figures/fig1_concept.pdf"))

