#' Convert SpatVector to Spatial*
#' 
#' This function converts a \code{SpatVector} object from the \pkg{terra} package to a \code{Spatial} object of the appropriate class (\code{SpatialPoints}, \code{SpatialPointsDataFrame}, \code{SpatialPolygons}, or \code{SpatialPolygonsDataFrame}) from the \pkg{sp} package. Note that \pkg{sp} is to be retired in 2023, so this function is to be come useful only for legacy applications.
#'
#' @param x		\code{SpatVector} object.
#' @return Object of class \code{Spatial}.
#' @examples
#'
#' f <- system.file('ex/lux.shp', package='terra')
#' v <- terra::vect(f)
#' spat <- svToSpatial(v)
#' class(spat)
#'
#' @export

svToSpatial <- function(x) {

	x <- sf::st_as_sf(x)
	x <- as(x, 'Spatial')
	x

}
