#' @name madClim2030
#'
#' @title Future climate rasters for Madagascar
#'
#' @description Rasters representing average climate across 2021-2040 modeled with CanESM5 for SSP 585 for Madagascar from WorldClim version 2.1. Values of these rasters have been rounded to one digit, so \emph{please do not use these for "official" work}. Please also note that CanESM5 in CMIP6 is known to run "too hot", but is useful here to aid illustration.
#'
#' @docType data
#'
#' @format An object of class \code{'SpatRaster'}.
#'
#' @keywords climate Madagascar
#'
#' @source \href{https://worldclim.org}{WorldClim}
#' 
#' @examples
#'
#' library(terra)
#' rastFile <- system.file('extdata', 'madClim2030.tif', package='enmSdmX')
#' madClimFut <- rast(rastFile)
#' plot(madClimFut)
#'
NULL
