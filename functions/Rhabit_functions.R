## Rhabit functions
# View(Rhabit:::simLangevinMM)
simLangevinMM <- function (beta, gamma2 = 1, times, loc0, cov_list = NULL, grad_fun = NULL, 
          silent = F, keep_grad = F) 
{
  checkCovGrad(cov_list, grad_fun)
  nb_obs <- length(times)
  xy <- matrix(NA, nb_obs, 2)
  xy[1, ] <- loc0
  dt <- diff(times)
  J <- length(beta)
  grad_array <- NULL
  if (keep_grad) {
    grad_array <- matrix(ncol = 2 * J, nrow = nb_obs)
    colnames(grad_array) <- paste(rep(paste0("grad_c", 1:J), 
                                      rep(2, J)), rep(c("x", "y"), J), sep = "_")
  }
  computeGradient <- function(loc, cov_list) {
    bilinearGrad(loc, cov_list)
  }
  if (!is.null(grad_fun)) {
    computeGradient <- function(loc, cov_list) {
      sapply(grad_fun, function(foo) {
        foo(loc)
      })
    }
  }
  for (t in 2:nb_obs) {
    if (!silent) 
      cat("\rSimulating Langevin process...", round(100 * 
                                                      t/nb_obs), "%")
    cov_list_tmp <- lapply(cov_list, getGridZoom, x0 = xy[t - 
                                                            1, ])
    grad_val <- computeGradient(loc = xy[t - 1, ], cov_list = cov_list_tmp)
    if (keep_grad) 
      grad_array[t - 1, ] <- as.numeric(grad_val)
    grad <- grad_val %*% beta
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 
                                                       1]))
    xy[t, ] <- xy[t - 1, ] + 0.5 * grad * gamma2 * dt[t - 
                                                        1] + rand_part
  }
  if (!silent) 
    cat("\n")
  main_df <- data.frame(x = xy[, 1], y = xy[, 2], t = times)
  if (keep_grad) {
    cov_list_tmp <- lapply(cov_list, getGridZoom, x0 = xy[nb_obs - 
                                                            1, ])
    grad_val <- computeGradient(loc = xy[t - 1, ], cov_list = cov_list_tmp)
    grad_array[nb_obs, ] <- as.numeric(grad_val)
    main_df <- cbind.data.frame(main_df, grad_array)
  }
  return(main_df)
}

# View(Rhabit:::checkCovGrad)
checkCovGrad <- function (cov_list, grad_fun) 
{
  if (is.null(cov_list) & is.null(grad_fun)) {
    stop("Either cov_list of grad_fun must be non NULL")
  }
  else if (is.null(cov_list)) {
    checkList(grad_fun, "grad_fun")
  }
  else if (is.null(grad_fun)) {
    checkList(cov_list, "cov_list")
  }
  else {
    check2Lists(cov_list, grad_fun)
  }
}

# View(Rhabit:::checkList)
checkList <- function (my_list, name = "my_list") 
{
  if (!inherits(my_list, "list")) {
    stop(paste(name, "should be a list"))
  }
  else {
    test_null <- sapply(my_list, is.null)
    if (any(test_null)) 
      stop(paste("Element(s):", paste(which(test_null), 
                                      collapse = " and "), "has/have null element in both cov_list and grad_fun"))
  }
  return(NULL)
}

# View(Rhabit:::check2Lists)
check2Lists <- function (list1, list2) 
{
  if (!(inherits(list1, "list") & inherits(list2, "list"))) {
    stop("cov_list and grad_fun should be either NULL or lists")
  }
  else {
    test_null <- mapply(function(x, y) {
      is.null(x) & is.null(y)
    }, list1, list2)
    if (any(test_null)) 
      stop(paste("Element(s):", paste(which(test_null), 
                                      collapse = " and "), "has/have null element in both cov_list and grad_fun"))
  }
  return(NULL)
}

# View(Rhabit:::bilinearGrad)
bilinearGrad <- function (loc, cov_list) 
{
  J <- length(cov_list)
  grad_val <- sapply(1:J, function(j) {
    x_grid <- cov_list[[j]]$x
    y_grid <- cov_list[[j]]$y
    cov_mat <- cov_list[[j]]$z
    cell <- gridCell(loc = loc, xgrid = x_grid, ygrid = y_grid, 
                     covmat = cov_mat)
    x <- cell$coords[1:2]
    y <- cell$coords[3:4]
    f <- cell$values
    dfdx <- ((y[2] - loc[2]) * (f[2, 1] - f[1, 1]) + (loc[2] - 
                                                        y[1]) * (f[2, 2] - f[1, 2]))/((y[2] - y[1]) * (x[2] - 
                                                                                                         x[1]))
    dfdy <- ((x[2] - loc[1]) * (f[1, 2] - f[1, 1]) + (loc[1] - 
                                                        x[1]) * (f[2, 2] - f[2, 1]))/((y[2] - y[1]) * (x[2] - 
                                                                                                         x[1]))
    return(c(dfdx, dfdy))
  })
  return(grad_val)
}

# View(Rhabit:::gridCell)
gridCell <- function (loc, xgrid, ygrid, covmat) 
{
  ix <- findInterval(loc[1], xgrid)
  iy <- findInterval(loc[2], ygrid)
  coords <- c(xgrid[ix], xgrid[ix + 1], ygrid[iy], ygrid[iy + 
                                                           1])
  values <- covmat[ix:(ix + 1), iy:(iy + 1)]
  return(list(coords = coords, values = values))
}
# View(Rhabit:::getGridZoom)
getGridZoom <- function (covar, x0, lag_inter = 2) 
{
  if (is.null(covar)) 
    return(NULL)
  else {
    x_pos <- findInterval(x0[1], covar$x)
    y_pos <- findInterval(x0[2], covar$y)
    x_inter <- max(1, x_pos - lag_inter):min(length(covar$x), 
                                             x_pos + lag_inter)
    y_inter <- max(1, y_pos - lag_inter):min(length(covar$y), 
                                             y_pos + lag_inter)
    return(list(x = covar$x[x_inter], y = covar$y[y_inter], 
                z = covar$z[x_inter, y_inter]))
  }
}


# bilinearGradArray from Rhabit
bilinearGradArray <- function (locs, cov_list) {
  if (!inherits(locs, "matrix")) 
    stop("'locs' must be a matrix")
  grad <- unlist(lapply(1:nrow(locs), function(i) bilinearGrad(loc = locs[i, 
    ], cov_list = cov_list)))
  gradarray <- array(grad, c(2, length(cov_list), nrow(locs)))
  gradarray <- aperm(gradarray, c(3, 1, 2))
  return(gradarray)
}

#langevinUD
langevinUD <- function (locs, times, ID = NULL, grad_array, with_speed = TRUE, 
  alpha = 0.95, leverage = FALSE) 
{
  if (!(inherits(locs, "matrix") & typeof(locs) %in% c("double", 
    "integer"))) 
    stop("locs must be a numeric matrix")
  if (inherits(grad_array, "matrix")) {
    grad_array <- array(grad_array, dim = c(n, 2, 1))
    warning("gradsArray was a matrix, and has been transformed to an array")
  }
  n <- nrow(locs)
  if (is.null(ID)) 
    ID <- rep(1, n)
  if (length(times) != n) 
    stop("Length of times must be the nrow of locs")
  if (length(ID) != n) 
    stop("Length of ID must be the nrow of locs")
  if (dim(grad_array)[1] != n) 
    stop("The first dimension of gradientArray must be of size nrow(locs)")
  J <- dim(grad_array)[3]
  break_ind <- which(ID[-1] != ID[-n])
  start_ind <- c(1, break_ind + 1)
  end_ind <- c(break_ind, n)
  grad_mat <- 0.5 * rbind(grad_array[-end_ind, 1, ], grad_array[-end_ind, 
    2, ])
  time_lag <- rep(times[-start_ind] - times[-end_ind], 2)
  sq_time_lag <- sqrt(time_lag)
  loc_increment <- as.numeric(locs[-start_ind, ] - locs[-end_ind, 
    ])
  Y <- loc_increment/sq_time_lag
  nu_hat_var <- solve(t(grad_mat * time_lag) %*% grad_mat)
  DF <- length(Y) - J
  if (DF <= 4) {
    stop("nrow(locs) must be strictly larger than 2 + dim(grad_array)[3] / 2")
  }
  nu_hat <- nu_hat_var %*% t(grad_mat) %*% loc_increment
  predictor <- sq_time_lag * grad_mat %*% nu_hat
  if (with_speed) {
    gamma2_hat <- colSums((Y - predictor)^2)/DF
    beta_hat <- nu_hat/gamma2_hat * (DF - 2)/DF
    beta_hat_var <- (2 * beta_hat %*% t(beta_hat)/(DF - 
      4) + nu_hat_var/gamma2_hat * (1 + 2/(DF - 4)))
  }
  else {
    beta_hat <- nu_hat
    beta_hat_var <- nu_hat_var
  }
  quantlow <- (1 - alpha)/2
  quantup <- 1 - (1 - alpha)/2
  conf_interval_beta <- t(sapply(1:length(beta_hat), function(j) {
    beta_hat[j] + c(1, -1) * stats::qnorm(quantlow) * sqrt(beta_hat_var[j, 
      j])
  }))
  conf_interval_gamma2 <- gamma2_hat * DF/stats::qchisq(c(quantup, 
    quantlow), DF)
  conf_interval <- rbind(conf_interval_beta, conf_interval_gamma2)
  rownames(beta_hat_var) <- colnames(beta_hat_var) <- paste0("beta", 
    1:J)
  rownames(conf_interval) <- c(rownames(beta_hat_var), "gamma2")
  colnames(conf_interval) <- c(quantlow, quantup)
  r_square <- 1 - colSums((Y - predictor)^2)/sum(Y^2)
  res <- matrix(Y - predictor, ncol = 2)
  if (leverage) {
    Z <- sq_time_lag * grad_mat
    H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    lever <- diag(H)
    res <- res/(sqrt(gamma2_hat * (1 - diag(H))))
  }
  else {
    lever <- NULL
  }
  AIC <- AICEuler(beta = as.numeric(beta_hat), gamma2 = gamma2_hat, 
    locs = locs, times = times, ID = ID, grad_array = grad_array)
  return(list(betaHat = as.numeric(beta_hat), gamma2Hat = gamma2_hat, 
    betaHatVariance = beta_hat_var, CI = conf_interval, 
    predicted = matrix(predictor, ncol = 2), R2 = r_square, 
    residuals = res, lever = lever, AIC = AIC))
}

# AICEuler
AICEuler <- function (beta, gamma2 = 1, locs, times, ID = NULL, grad_array) 
{
  npar <- length(beta)
  if (gamma2 != 1) 
    npar <- npar + 1
  nllk <- nllkEuler(beta = beta, gamma2 = gamma2, locs = locs, 
    times = times, ID = ID, grad_array = grad_array)
  AIC <- 2 * (npar + nllk)
  return(AIC)
}


# nllkEuler
nllkEuler <- function (beta, gamma2 = 1, locs, times, ID = NULL, grad_array) 
{
  n <- nrow(locs)
  gradmat <- 0.5 * gamma2 * apply(grad_array, 2, function(mat) mat %*% 
    beta)
  if (is.null(ID)) 
    ID <- rep(1, nrow(locs))
  break_ind <- which(ID[-1] != ID[-n])
  start_ind <- c(1, break_ind + 1)
  end_ind <- c(break_ind, n)
  dt <- times[-start_ind] - times[-end_ind]
  llk <- sum(dnorm(locs[-start_ind, ], locs[-end_ind, ] + 
    dt * gradmat[-end_ind, ], sqrt(gamma2 * dt), log = TRUE))
  return(-llk)
}
