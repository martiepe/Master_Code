#' spatRast_to_covlist
#'
#' @param rast terra spatRaster Description of argument
#'
#' @return returns covariates as list formatted for BBIS fitting
#'
#' @export
spatRast_to_covlist <- function (rast) 
 {
     covlist <- list()
     xgrid <- terra::xFromCol(rast, col = 1:terra::ncol(rast))
     ygrid <- terra::yFromRow(rast, row = 1:terra::nrow(rast))
     for (i in 1:nlyr(rast)) {
         layer_matrix <- terra::as.matrix(rast[[i]], wide = TRUE)
         covlist[[i]] <- list(x = xgrid, y = ygrid, z = t(layer_matrix))
     }
     return(covlist)
 }

