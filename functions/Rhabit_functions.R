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


