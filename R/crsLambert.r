#' Generate a custom Lambert azimuthal equal-area projection WKT2 string
#'
#' This function generates a WKT2 (well-known text) coordinate reference string for a Lambert equal-area azimuthal projection of a geographic object or coordinate pair.
#'
#' @param x 	Either an object of class \code{SpatRaster}, \code{SpatVector}, or \code{sf}, \emph{or} a numeric string with two values (longitude and latitude of the center of the projection), \emph{or} a two-column matrix/data frame with the centroid of the projection.
#'
#' @seealso \code{\link{crsGet}}
#' @examples
#'
#' library(sf)
#' data(mad0)
#' wkt2 <- crsLambert(mad0)
#' mad0Lamb <- st_transform(mad0, wkt2)
#'
#' par(mfrow=c(1, 2))
#' plot(st_geometry(mad0), main='Unprojected (WGS84)')
#' plot(st_geometry(mad0Lamb), main='Lambert')
#'
#' # The effect is more noticable when plotting large areas,
#' # especially near the poles.
#' if (FALSE) {
#' 
#' library(geodata)
#' library(terra)
#'
#' can <- gadm('CAN', level=0, path=getwd()) # outline of Canada
#'
#' wktVert <- crsVertical(can)
#' wktLamb <- crsLambert(can)
#'
#' canVert <- project(can, wktVert)
#' canLamb <- project(can, wktLamb)
#'
#' par(mfrow=c(1, 3))
#' plot(can, main='WGS84')
#' plot(canVert, main='Vertical')
#' plot(canLamb, main='Lambert')
#'
#' }
#'
#' @export

crsLambert <- function(x) {

	if (inherits(x, c('SpatRaster', 'SpatVector', 'sf'))) {
	
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
	
	} else {
		stop('Argument "x" must be a coordinate pair or a spatial object of class "SpatRaster", "SpatVector", or "sf".')
	}

	out <- paste0(
		'PROJCRS["Lambert_Conformal_Conic",
    BASEGEOGCRS["NAD83",
        DATUM["North American Datum 1983",
            ELLIPSOID["GRS 1980",6378137,298.257222101,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4269]],
    CONVERSION["North_America_Lambert_Conformal_Conic",
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
