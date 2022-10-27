#' Convert extent to polygon
#'
#' This function returns a \code{SpatVector} polygon representing an extent. The input can be a \code{SpatExtent} object or an object from which a \code{SpatExtent} (extent) can be obtained.
#'
#' @param x A \code{sf}, \code{SpatVector}, or \code{SpatRaster} object, \emph{or} a vector of four numeric values representing (in this order): x-coordinate of western side of the extent, x-coordinate of eastern side, y-coordinate of the southern side, and y-coordinate of the northern side. If numeric coordinates are supplied, the output will not have a CRS assigned to it unless supplied in \code{...}
#' @param ... Arguments to supply to \code{vect}.
#'
#' @return A \code{SpatVector}.
#' @seealso \link{createPlotPoly}
#' @examples
#'
#' data(mad0)
#' madExtent <- extentToVect(mad0)
#' plot(madExtent, border='blue', lty='dotted')
#' plot(mad0, add=TRUE)
#'
#' @export

extentToVect <- function(x, ...) {

	if (inherits(x, c('numeric'))) {

		xCorners <- c(x[1L], x[2L], x[2L], x[1L])
		yCorners <- c(x[3L], x[3L], x[4L], x[4L])
		corners <- cbind(xCorners, yCorners)
		out <- terra::vect(corners, type='polygon', ...)
		
	} else if (inherits(x, c('SpatVector', 'SpatRaster'))) {

		proj <- terra::crs(x)

		x <- terra::ext(x)@ptr$vector
		xCorners <- c(x[1L], x[2L], x[2L], x[1L])
		yCorners <- c(x[3L], x[3L], x[4L], x[4L])
		corners <- cbind(xCorners, yCorners)
		out <- terra::vect(corners, type='polygon', crs=proj, ...)
	
	}
	
	out
	
}
