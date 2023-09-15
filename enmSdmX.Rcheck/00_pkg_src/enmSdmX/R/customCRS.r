#' Custom coordinate reference system WKT2 string
#'
#' These functions take as input either a spatial object or coordinate pair and a custom WKT2 (well-known text) coordinate reference system string centered on the object or coordinate. Projections include:
#' \itemize{
#'	\item Albers conic equal-area
#'	\item Lambert azimuthal equal-area
#'	\item Vertical near-side (i.e., as the world appears from geosynchronous orbit)
#' }
#' Please note that these are \emph{NOT} standard projections, so do not have an EPSG or like code.
#'
#' @param x 	Either an object of class \code{SpatRaster}, \code{SpatVector}, or \code{sf}, \emph{or} a numeric vector with two values (longitude and latitude of the center of the projection), \emph{or} a two-column matrix/data frame with the centroid of the projection.
#'
#' @return A WKT2 (well-known text) string.
#'
#' @seealso \code{\link{getCRS}}, \code{\link{customAlbers}}, \code{\link{customLambert}}, \code{\link{customVNS}}
#'
#' @example man/examples/customCRS_examples.r
#'
#' @export

customAlbers <- function(x) {

	cent <- .getCentroid(x)
	long <- cent$long
	lat <- cent$lat

	out <- paste0(
		'PROJCRS["Albers_Equal_Area_Conic",
    BASEGEOGCRS["NAD83",
        DATUM["North American Datum 1983",
            ELLIPSOID["GRS 1980",6378137,298.257222101,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4269]],
    CONVERSION["North_America_Albers_Equal_Area_Conic",
        METHOD["Albers Equal Area",
            ID["EPSG",9822]],
        PARAMETER["Latitude of false origin",', lat, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8821]],
        PARAMETER["Longitude of false origin",', long, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8822]],
        PARAMETER["Latitude of 1st standard parallel",', lat - 20, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8823]],
        PARAMETER["Latitude of 2nd standard parallel",', lat + 20, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8824]],
        PARAMETER["Easting at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8826]],
        PARAMETER["Northing at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8827]]],
    CS[Cartesian,2],
        AXIS["(E)",east,
            ORDER[1],
            LENGTHUNIT["metre",1]],
        AXIS["(N)",north,
            ORDER[2],
            LENGTHUNIT["metre",1]]]'
	)
	
	out
	
}

#' @describeIn customAlbers Custom coordinate reference system WKT2 string
#' @export
customLambert <- function(x) {

	cent <- .getCentroid(x)
	long <- cent$long
	lat <- cent$lat

	out <- paste0(
		'PROJCRS["Lambert_Conformal_Conic",
    BASEGEOGCRS["NAD83",
        DATUM["North American Datum 1983",
            ELLIPSOID["GRS 1980",6378137,298.257222101,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4269]],
    CONVERSION["Lambert_Conformal_Conic",
        METHOD["Lambert Conic Conformal (2SP)",
            ID["EPSG",9802]],
        PARAMETER["Latitude of false origin",', lat, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8821]],
        PARAMETER["Longitude of false origin",', long, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8822]],
        PARAMETER["Latitude of 1st standard parallel",', lat - 20, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8823]],
        PARAMETER["Latitude of 2nd standard parallel",', lat + 20, ',
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8824]],
        PARAMETER["Easting at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8826]],
        PARAMETER["Northing at false origin",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8827]]],
    CS[Cartesian,2],
        AXIS["(E)",east,
            ORDER[1],
            LENGTHUNIT["metre",1]],
        AXIS["(N)",north,
            ORDER[2],
            LENGTHUNIT["metre",1]]]'
	)
	
	out
	
}

#' @describeIn customAlbers Custom coordinate reference system WKT2 string
#' @param alt	Altitude in meters of the viewpoint in km. The default (35800 km) is geosynchronous orbit.
#' @export
customVNS <- function(x, alt = 35800) {

	alt <- 1000 * alt

	cent <- .getCentroid(x)
	long <- cent$long
	lat <- cent$lat

	out <- paste0(
		'PROJCRS["World_Vertical_Perspective",
			BASEGEOGCRS["WGS 84",
			DATUM["World Geodetic System 1984",
				ELLIPSOID["WGS 84",6378137,298.257223563,
					LENGTHUNIT["metre",1]]],
			PRIMEM["Greenwich",0,
				ANGLEUNIT["Degree",0.0174532925199433]]],
		CONVERSION["World_Vertical_Perspective",
			METHOD["Vertical Perspective",
				ID["EPSG",9838]],
			PARAMETER["Latitude of topocentric origin",', lat, ',
				ANGLEUNIT["Degree",0.0174532925199433],
				ID["EPSG",8834]],
			PARAMETER["Longitude of topocentric origin",', long, ',
				ANGLEUNIT["Degree",0.0174532925199433],
				ID["EPSG",8835]],
			PARAMETER["Viewpoint height",', alt, ',
				LENGTHUNIT["metre",1],
				ID["EPSG",8840]]],
		CS[Cartesian,2],
			AXIS["(E)",east,
				ORDER[1],
				LENGTHUNIT["metre",1]],
			AXIS["(N)",north,
				ORDER[2],
				LENGTHUNIT["metre",1]],
		USAGE[
			SCOPE["Not known."],
			AREA["World."],
			BBOX[-90,-180,90,180]]]'
	)
	
	out
	
}

customVNS <- function(x, alt = 35800) {

	alt <- 1000 * alt

	cent <- .getCentroid(x)
	long <- cent$long
	lat <- cent$lat

	out <- paste0(
		'PROJCRS["World_Vertical_Perspective",
			BASEGEOGCRS["WGS 84",
			DATUM["World Geodetic System 1984",
				ELLIPSOID["WGS 84",6378137,298.257223563,
					LENGTHUNIT["metre",1]]],
			PRIMEM["Greenwich",0,
				ANGLEUNIT["Degree",0.0174532925199433]]],
		CONVERSION["World_Vertical_Perspective",
			METHOD["Vertical Perspective",
				ID["EPSG",9838]],
			PARAMETER["Latitude of topocentric origin",', lat, ',
				ANGLEUNIT["Degree",0.0174532925199433],
				ID["EPSG",8834]],
			PARAMETER["Longitude of topocentric origin",', long, ',
				ANGLEUNIT["Degree",0.0174532925199433],
				ID["EPSG",8835]],
			PARAMETER["Viewpoint height",', alt, ',
				LENGTHUNIT["metre",1],
				ID["EPSG",8840]]],
		CS[Cartesian,2],
			AXIS["(E)",east,
				ORDER[1],
				LENGTHUNIT["metre",1]],
			AXIS["(N)",north,
				ORDER[2],
				LENGTHUNIT["metre",1]],
		USAGE[
			SCOPE["Not known."],
			AREA["World."],
			BBOX[-90,-180,90,180]]]'
	)
	
	out
	
}


### get centroid of raster, vector, sf, etc. object
###################################################
.getCentroid <- function(x) {

	# x		raster, vector, data frame, matrix, etc.

	if (inherits(x, c('SpatRaster', 'SpatVector', 'sf', 'Spatial'))) {
	
		x <- if (inherits(x, 'SpatRaster')) {
			extentToVect(x)
		} else if (inherits(x, 'SpatVector')) {
			terra::centroids(x)
		} else if (inherits(x, 'sf')) {
			sf::st_centroid(sf::st_geometry(x)) 
		}
		
		x <- sf::st_as_sf(x)
		cent <- sf::st_centroid(sf::st_geometry(x))
		cent <- sf::st_coordinates(cent)
		
		long <- cent[1L, 1L]
		lat <- cent[1L, 2L]
	
	} else if (inherits(x, c('matrix', 'data.frame'))) {
		if (ncol(x) != 2L | nrow(x) != 1L) stop('Argument "x" must be a 2-column data frame/matrix with one row, a two-element numeric vector, or a spatial object.')
		long <- x[1L, 1L]
		lat <- x[1L, 2L]
	
	} else if (inherits(x, 'numeric')) {
		long <- x[1L]
		lat <- x[2L]
	} else {
		stop('Argument "x" must be a 2-column data frame/matrix with one row, a two-element numeric vector, or a spatial object.')
	}

	list(long = long, lat = lat)

}

