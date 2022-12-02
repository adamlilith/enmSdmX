#' @name madEnv
#'
#' @title Climate rasters for Madagascar
#'
#' @description Climate rasters for Madagascar from WorldClim version 1.4. Values of these rasters have been rounded, so please do not use them for "official" work.
#'
#' @docType data
#'
#' @usage rastFile <- system.file('extdata', 'madEnv.tif', package='enmSdmX')
#' madEnv <- rast(rastFile)
#'
#' @format An object of class \code{'SpatRaster'}.
#'
#' @keywords datasets, Madagascar
#'
#' @source \href{WorldClim}{www.worldclim.org}
#' 
#' @examples
#' rastFile <- system.file('extdata', 'madEnv.tif', package='enmSdmX')
#' madEnv <- rast(rastFile)
#' plot(madEnv)
NULL
