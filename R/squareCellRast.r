#' Create a raster with square cells
#' 
#' This function creates a raster from an object with an extent (i.e., another raster or similar spatial object) with square cells. The user can specify cell resolution (linear dimension) \emph{or} the approximate number of cells desired.
#'
#' @param x An object with a spatial extent property (e.g., a \code{SpatRaster} or a \code{SpatVector}).
#' @param numCells Positive integer, approximate number of cells desired. If this is specified, then \code{res} is ignored. If this number of cells cannot be fit into the desired extent exactly, then the actual number of cells will be larger.
#' @param res Positive numeric. Size of a cell in the units of the projection of \code{x} (typically meters). Ignored if \code{numCells} is not \code{NULL}.
#' @param vals Numeric, value to assign to cells. Note that if this is shorter than the number of cells in the output, then values will be recycled. If longer, then values will be truncated. The default is to assign all 0s.
#'
#' @return \code{SpatRaster} object. The raster will have an extent of the same size or larger than the extent of \code{x}.
#'
#' @examples
#'
#' library(sf)
#' library(terra)
#'
#' # project outline of Madagascar to equal-area:
#' data(mad0)
#' mad0Ea <- st_transform(mad0, getCRS('madAlbers'))
#'
#' n <- 101
#' cellSize_meters <- 10E4
#' byNumCells <- squareCellRast(mad0Ea, numCells=n)
#' byCellSize <- squareCellRast(mad0Ea, res=cellSize_meters)
#' 
#' oldPar <- par(mfrow=c(1, 2))
#'
#' main1 <- paste0('Cells: ', n, ' desired, ', ncell(byNumCells), ' actual')
#' plot(byNumCells, main = main1)
#' plot(mad0Ea, add = TRUE)
#'
#' main2 <- paste0('Cells ', cellSize_meters, ' m on a side')
#' plot(byCellSize, main = main2)
#' plot(mad0Ea, add = TRUE)
#'
#' par(oldPar)
#'
#' # Note that in this example they look the same, but the one on the left
#' # has one less row than the one on the right.
#'
#' @export

squareCellRast <- function(
	x,
	numCells = NULL,
	res = NULL,
	vals = NULL
) {

	if (is.null(numCells) & is.null(res)) stop('Either "numCells" or "res" must be specified.')
	if (!is.null(numCells) & !is.null(res)) warning('Both "numCells" and "res" are specified. Ignoring argument "res".')
	
	if (inherits(x, 'sf')) x <- terra::vect(x)
	
	crs <- terra::crs(x)
	ext <- terra::ext(x)@ptr$vector
	names(ext) <- c('xmin', 'xmax', 'ymin', 'ymax')
		
	longDist <- ext['xmax'] - ext['xmin']
	latDist <- ext['ymax'] - ext['ymin']

	# create raster with approximate number of desired cells
	if (!is.null(numCells)) res <- sqrt((longDist * latDist) / numCells)

	# number of rows/columns
	numRows <- ceiling(latDist / res)
	numCols <- ceiling(longDist / res)

	# expand extent to accommodate rows/columns
	padLat <- res * (numRows - (latDist / res)) / 2
	padLong <- res * (numCols - (longDist / res)) / 2
	
	ext <- c(
		ext['xmin'] - padLong,
		ext['xmax'] + padLong,
		ext['ymin'] - padLat,
		ext['ymax'] + padLat
	)
	
	ext <- terra::ext(ext)
	out <- terra::rast(ext, nrows=numRows, ncols=numCols, crs=crs)
	numCells <- terra::ncell(out)
	vals <- if (is.null(vals)) {
		rep(0, numCells)
	} else {
		rep(vals, length.out=numCells)
	}
	terra::values(out) <- vals
	out
	
}
