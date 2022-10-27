#' Sample random points from a raster with/out replacement
#'
#' This function returns coordinates randomly located on a raster where cells can be sampled with replacement (if desired) and where the probability of selection is proportionate to the cell value, cell area, or the product of cell value times cell area.
#'
#' @param x \code{SpatRaster} object.
#'
#' @param n Positive integer. Number of points to draw.
#' @param adjArea If \code{TRUE} (default) then adjust probabilities so sampling accounts for cell area.
#' @param replace If \code{TRUE} (default) then sample with replacement.
#' @param prob If \code{TRUE} (default) then sample cells with probabilities proportional to cell values. If \code{adjArea} is also \code{TRUE} then probabilities are drawn proportional to the product of cell area * the value of the cell.
#'
#' @return 2-column matrix with longitude and latitude of random points.
#'
#' @seealso \code{\link[terra]{spatSample}}
#'
#' @examples
#' library(terra)
#' r <- rast()
#' nc <- ncell(r)
#' r[] <- 1:nc
#' rands1 <- sampleRast(r, 10000)
#' rands2 <- sampleRast(r, 10000, adjArea=FALSE)
#' rands3 <- sampleRast(r, 10000, prob=FALSE)
#' rands4 <- sampleRast(r, 10000, adjArea=FALSE, prob=FALSE)
#' par(mfrow=c(2, 2))
#' plot(r, main='adjArea = TRUE & prob = TRUE')
#' points(rands1, pch='.')
#' plot(r, main='adjArea = FALSE & prob = TRUE')
#' points(rands2, pch='.')
#' plot(r, main='adjArea = TRUE & prob = FALSE')
#' points(rands3, pch='.')
#' plot(r, main='adjArea = FALSE & prob = FALSE')
#' points(rands4, pch='.')
#'
#' @export

sampleRast <- function(x, n, adjArea = TRUE, replace = TRUE, prob = TRUE) {

	if (!inherits(x, 'SpatRaster')) x <- terra::rast(x)
	val <- as.vector(x)

	# adjust probabilities for cell area and/or cell values
	if (adjArea) {

		areas <- terra::cellSize(x, mask=TRUE)
		areas <- areas * (x * 0 + 1) # because "mask" argument does not work as documented 2021-12-26

		areas <- as.vector(areas)
		prob <- if (prob) {
			val * areas
		} else {
			areas
		}

	} else if (!adjArea & prob) {

		if (any(val < 0, na.rm=TRUE)) stop('Some probabilities are < 0.')
		prob <- val

	} else if (!adjArea & !prob) {

		prob <- rep(1, length(val))

	}

	cellNum <- 1:terra::ncell(x)
	cellNum <- cellNum[!is.na(val)]
	prob <- prob[!is.na(val)]

	# centers of cell
	sites <- sample(cellNum, size=n, replace=replace, prob=prob)
	xy <- terra::xyFromCell(x, sites)
	
	# move within cells
	from <- -0.5 * terra::res(x)[1]
	to <- 0.5 * terra::res(x)[1]
	rands1 <- runif(n, from, to)
	from <- -0.5 * terra::res(x)[2]
	to <- 0.5 * terra::res(x)[2]
	rands2 <- runif(n, from, to)

	xy[ , 'x'] <- xy[ , 'x'] + rands1
	xy[ , 'y'] <- xy[ , 'y'] + rands2

	xy

}
