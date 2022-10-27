#' Create SpatialPolygon same size as a plot
#'
#' This function creates a "rectangular" \code{SpatVector} object with the same dimensions as a plot window. It is especially useful for cropping subsequent rasters or \code{SpatVector} objects to the plot window. A plot must be made before calling this function.
#'
#' @param x Either \code{NULL} (default), an object of class \code{crs}, a coordinate reference string (proj6 WKT string), or an object with a coordinate reference system. If any of these is provided, the \code{SpatVector} object will have this CRS.
#'
#' @return \code{SpatVector}
#' @seealso \link{extentToVect}
#' @examples
#'
#' data(mad0)
#' poly <- makePlotPoly(mad0)
#' plot(poly, border='blue', lty='dotted')
#' plot(mad0, add=TRUE)
#'
#' @export

makePlotPoly <- function(x = NULL) {

	usr <- par('usr')
	extent <- terra::ext(usr)
	corners <- extent@ptr$vector
	xCoords <- c(corners[1], corners[2], corners[2], corners[1])
	yCoords <- c(corners[3], corners[3], corners[4], corners[4])
	xy <- matrix(c(xCoords, yCoords), ncol=2)
	poly <- terra::vect(xy, type='polygon')
	if (!is.null(x)) {
		if (inherits(x, c('crs', 'character'))) {
			terra::project(poly) <- x
		} else {
			terra::crs(poly) <- terra::crs(x)
		}
	}
	
	poly
	
}
