#' Get or assign values to cells in a raster
#'
#' These functions get values from a raster at specific cells, or values to specific cells.
#'
#' @param x A \code{SpatRaster}.
#' @param val One or more values. If more the number of cells specified is greater than the number of values in \code{val}, then values in \code{val} will be recycled.
#' @param cell Cell indices. There must be one per value in \code{val}.
#' @param format The type of cell indexing used. This can be either "raster" for row indexing (default) or "matrix" for column indexing. Row indexing (the default for rasters), starts with cell "1" in the upper left, cell "2" is to its right, and so on. Numbering then wraps around to the next row. Column indexing (the default for matrices) has the cell "1" in the upper left corner of the matrix. The cell "2" is below it, and so on. The numbering then wraps around to the top of the next column.
#'
#' @return A data frame (\code{getValueByCell}) with cell numbers (in row format), or a \code{SpatRaster} (\code{setValueByCell}).
#'
#' @seealso \code{\link[terra]{setValues}}, \code{\link[terra]{values}}
#'
#' @examples
#'
#' x <- rast(nrow=10, ncol=10)
#' x[] <- round(10 * runif(100))
#' 
#' cell=c(1, 20, 40, 80)
#' getValueByCell(x, cell = cell)
#' getValueByCell(x, cell = cell, format = 'matrix')
#' 
#' y <- setValueByCell(x, val = 20, cell = cell)
#' plot(y)
#' z <- setValueByCell(x, val = 30, cell = cell, format = 'matrix')
#' 
#' plot(c(x, y, z))
#' 
#' @export

getValueByCell <- function(x, cell, format = 'raster') {

	# convert to row cells
	if (format == 'matrix') {
		dims <- dim(x)[1:2]
		cell <- omnibus::rowColIndexing(dims, cell=cell, dir='row')
	}
	
	# extract
	x <- as.data.frame(x, cell=TRUE)
	x[x$cell %in% cell, , drop=FALSE]

}

setValueByCell <- function(x, val, cell, format = 'raster') {

	# errors
	if (terra::nlyr(x) > 1L) stop('"x" can contain only one raster.')

	lv <- length(val)
	lc <- length(cell)
	if (lv > lc) stop('The length of "val" must be as long as or shorter than the length of "cell".')
		
	# values
	if (lv < lc) val <- rep(val, length.out=lc)

	# order
	ord <- order(cell)
	cell <- cell[ord]
	val <- val[ord]

	# convert to row cells
	dims <- dim(x)[1:2]
	if (format == 'matrix') {
		cell <- omnibus::rowColIndexing(dims, cell=cell, dir='row')
	}

	m <- as.data.frame(x, cell=TRUE)
	m[m$cell %in% cell, 2] <- val
	out <- as.matrix(x, wide=TRUE)
	m$newCell <- omnibus::rowColIndexing(dims, cell=m$cell, dir='col')
	m <- m[order(m$newCell), ]
	
	out[m$newCell] <- m[ , 2L, drop=TRUE]
	out <- terra::rast(out, extent=terra::ext(x), crs=terra::crs(x))
	out

}
