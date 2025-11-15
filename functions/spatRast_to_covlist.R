#' spatRast_to_covlist
#'
#' @param rast terra spatRaster Description of argument
#'
#' @return returns covariates as list formatted for BBIS fitting
#'
#' @export
spatRast_to_covlist <- function(rast) {
  xgrid  <- terra::xFromCol(rast, col = 1:terra::ncol(rast))
  ygrid  <- terra::yFromRow(rast, row = terra::nrow(rast):1)
  
  covlist <- lapply(seq_len(terra::nlyr(rast)), function(i) {
    z <- terra::as.matrix(rast[[i]], wide = TRUE)
    z <- t(apply(z, 2, rev))  # flip vertically then transpose (matches raster::as.matrix path)
    list(x = xgrid, y = ygrid, z = z)
  })
  names(covlist) <- names(rast)
  covlist
}

