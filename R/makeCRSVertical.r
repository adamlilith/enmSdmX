#' Generate a "Vertical Near-Side" projection WKT2 string
#'
#' This function generates a WKT2 (well-known text) coordinate reference string for a "vertical near-side" projection of a geographic object or coordinate pair. When plotted, such an object will appear as if it would from a satellite in geosyncronous orbit directly above the centroid of the region of interest (showing about 1/3 of the Earth).
#'
#' @param x 	Either an object of class \code{SpatRaster}, \code{SpatVector}, or \code{sf}, \emph{or} a numeric string with two values (longitude and latitude of the center of the projection), \emph{or} a two-column matrix/data frame with the centroid of the projection.
#' @param alt	Altitude in meters of the viewpoint in km. The default (35800 km) is geosynchronous orbit.
#'
#' @seealso \link{getCRS}
#' @examples
#'
#' data(mad0)
#' wkt2 <- makeCRSVertical(mad0)
#' mad0vert <- st_transform(mad0, wkt2)
#'
#' par(mfrow=c(1, 2))
#' plot(mad0, main='Unprojected (WGS84)')
#' plot(mad0vert, main='Vertical')
#'
#' @export

makeCRSVertical <- function(x, alt = 35800) {

	alt <- 1000 * alt

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
