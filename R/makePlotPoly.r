#' Create spatial polygon same size as a plot
#'
#' This function creates a "rectangular" \code{SpatVector} object with the same dimensions as a plot window. It is especially useful for cropping subsequent rasters or vector objects to the plot window. A plot must be made before calling this function.
#'
#' @param x Either \code{NULL} (default), an object of class \code{crs}, a coordinate reference string (PROJ6 WKT string), or an object with a coordinate reference system. If any of these is provided, the \code{SpatVector} object will have this CRS.
#'
#' @return \code{SpatVector}
#' @seealso \link{extentToVect}
#' @examples
#'
#' if (FALSE) {
#'
#' library(sf)
#'
#' data(mad0)
#' plot(st_geometry(mad0))
#' outline <- makePlotPoly(mad0)
#' plot(outline, col='cornflowerblue', lty='dotted')
#' plot(st_geometry(mad0), add=TRUE)
#' 
#' }
#'
#' @export

makePlotPoly <- function(x = NULL) {

	corners <- graphics::par('usr')
	xCoords <- c(corners[1], corners[2], corners[2], corners[1])
	yCoords <- c(corners[3], corners[3], corners[4], corners[4])
	xy <- matrix(c(xCoords, yCoords), ncol=2)
	out <- terra::vect(xy, type='polygon')
	if (!is.null(x)) {
		
		if (inherits(x, 'sf')) {
			proj <- sf::st_crs(x)
			proj <- as.character(proj$wkt)

			out <- sf::st_as_sf(out)
			out <- sf::st_set_crs(out, proj)
		} else if (inherits(x, c('SpatVector', 'SpatRaster'))) {
			proj <- terra::crs(x)
			out <- terra::crs(out) <- proj
		} else if (inherits(x, c('crs', 'character'))) {
			proj <- x
		}
		
		
	}
	
	out
	
}
