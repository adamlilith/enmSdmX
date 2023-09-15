#' "Friendly" wrapper for terra::global() for calculating raster statistics
#'
#' Calculate "global" statistics across all the values in a raster. This function is a wrapper for \code{\link[terra]{global}}.  That function, by default, sets \code{na.rm = FALSE}, so any cell that is \code{NA} can cause the summary statistic to also be \code{NA} (usually undesirable). The function also returns a \code{data.frame}, so often needs a further line of code to get the actual value(s). This function sets \code{na.rm = TRUE} by default, and returns a numeric vector (not a \code{data.frame}).
#'
#' @param x A \code{SpatRaster}.
#' @param fun A function or the name of a function (in quotes). See \code{\link[terra]{global}} for more details.
#' @param na.rm If \code{TRUE} (default), then the function in \code{fun} will ignore \code{NA} cells.
#' @param ... Additional arguments to pass to \code{fun}.
#' @param weights Either \code{NULL} or a \code{SpatRaster}.
#' 
#' @return A numeric vector, one value per layer in \code{x}.
#'
#' @examples
#' 
#' library(terra)
#'
#' r <- rast(ncols=10, nrows=10)
#' values(r) <- 1:ncell(r)
#'
#' 
#' global(r, 'sum') # terra
#' globalx(r, 'sum') # enmSdmX
#'
#' global(r, "mean", na.rm=TRUE)[1, 1] # terra... same as enmSdmX::globalx
#' 
#' @export 
globalx <- function(x, fun, na.rm = TRUE, ..., weights = NULL) {

	out <- terra::global(x = x, fun = fun, na.rm = na.rm, ..., weights = weights)
	out <- unlist(out)
	out

}
