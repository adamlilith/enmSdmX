#' @name canada
#'
#' @title Vector outline of Canada
#'
#' @description This \code{SpatVector} represents the outline of Canada in WGS84 (unprojected) coordinates. This is the "low resolution" (less accurate) version from GADM.
#'
#' @docType data
#'
#' @format An object of class \code{'SpatVector'}.
#'
#' @keywords Canada
#'
#' @source \href{https://gadm.org/index.html}{Database of Global Administrative Areas (GADM)}
#' 
#' @examples
#'
#' library(terra)
#' canFile <- system.file('extdata', 'canada_level0_gadm41.gpkg', package='enmSdmX')
#' canada <- vect(canFile)
#' plot(canada)
#'
NULL
