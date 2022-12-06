#' Generate rasters with cell values equal to cell longitude or latitude
#'
#' This function generates a raster stack with two rasters, one with cell values equal to the cell's longitude and the other with cell values equal to the cell's latitude.
#'
#' @param x \code{SpatRaster} object. The output will have the same resolution, extent, and coordinate projection system as \code{x}.
#' @param m Any of:
#'	\itemize{
#'		\item \code{TRUE} (default): Calculate longitude and latitude only for cells that are not \code{NA}.
#'		\item \code{FALSE}: Calculate longitude and latitude for all cells.
#'		\item A \code{SpatRaster} object: Force any cells that are \code{NA} in this raster to also be \code{NA} in the output.
#'	}
#' @param filePath String or \code{NULL}. If a string, then this is the path (not including file name) to which to write the raster stack with longitude/latitude rasters. If \code{NULL} then no file is written.
#' @param ... Arguments to pass to \code{writeRaster} (if \code{filePath} is not \code{NULL}).
#' @return Object of class \code{SpatRaster}.
#' @examples
#' 
#' library(terra)
#'
#' # generate long/lat rasters for the world
#' x <- rast() # raster with 1 deg resolution and extent equal to entire world
#' x[] <- 1:ncell(x)
#' longLat <- rastLongLag(x)
#' plot(longLat)
#'
#' # demonstrate masking
#' # randomly force some cells to NA
#' v <- 1:ncell(x)
#' n <- 10000
#' v[sample(v, n)] <- NA
#' x[] <- v
#' longLatTRUE <- rastLongLag(x, m = TRUE)
#' longLatFALSE <- rastLongLag(x, m = FALSE)
#' rasts <- c(x, longLatTRUE, x, longLatFALSE)
#' names(rasts) <- c('x', 'long_m_TRUE', 'lat_m_TRUE',
#' 	'x', 'long_m_FALSE', 'lat_m_FALSE')
#' plot(rasts)
#'
#' @export

rastLongLag <- function(
	x,
	m = TRUE,
	filePath = NULL,
	...
) {

	if (terra::nlyr(x) > 1L) x <- x[[1]]

	# mask raster
	if (inherits(m, 'logical')) {
		if (m) {
			m <- x * 0L + 1L
		} else if (!m) {
			m <- x
			m[] <- 1L
		}
	} else {
		m <- x
		m <- m * 0L + 1L
	}

	# initiate lat/long rasters
	lat <- terra::rast(terra::ext(x), nrow=nrow(x), ncol=ncol(x))
	long <- terra::rast(terra::ext(x), nrow=nrow(x), ncol=ncol(x))

	if (!is.null(filePath)) {

		long <- terra::writeStart(x=long, filename=paste0(filePath, '/longitude'), ...)
		lat <- terra::writeStart(x=lat, filename=paste0(filePath, '/latitude'), ...)

	}

	# write rasters
	if (!is.null(filePath)) {

		# for each block, calculate latitude and longitude
		for (countRow in 1L:nrow(x)) {

			# initiate list to store vectors of lat/long
			theseLong <- rep(NA, ncol(x))
			theseLat <- rep(NA, ncol(x))

			# assign latitudes and longitudes
			theseLong <- terra::xFromCol(object=x, col=1:ncol(x))
			theseLat <- rep(terra::yFromRow(object=x, row=countRow), ncol(x))

			# mask
			theseLong <- theseLong * m[countRow, ]
			theseLat <- theseLat * m[countRow, ]

			# remember output
			terra::writeValues(x=long, v=theseLong, start=countRow)
			terra::writeValues(x=lat, v=theseLat, start=countRow)

		} # for each block

		# stop writing
		terra::crs(lat) <- terra::crs(long) <- terra::crs(x)

		names(long) <- 'longitude'
		names(lat) <- 'latitude'

		long <- terra::writeStop(long)
		lat <- terra::writeStop(lat)

	# do in memory--no writing!
	} else {

		long[] <- rep(terra::xFromCol(object=x, col=1:ncol(x)), nrow(x))
		lat[] <- rep(terra::yFromRow(object=x, row=1:nrow(x)), each=ncol(x))
	
		terra::crs(long) <- terra::crs(lat) <- terra::crs(x)
		
		long <- long * m
		lat <- lat * m

	}

	names(long) <- 'longitude'
	names(lat) <- 'latitude'

	ll <- c(long, lat)
	ll

}
