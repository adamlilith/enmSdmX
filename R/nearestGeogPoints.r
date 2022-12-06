#' Minimum convex polygon from a set of spatial polygons and/or points
#'
#' This function implements the "nearest geographic point" method (Smith et al. 2023) to enable the use of occurrence records geolocated only to a general place (e.g., a country or province), along with occurrences georeferenced with little error. The function returns a minimum convex polygon (MCP) constructed from a set of spatial polygons and/or points.
#'
#' @param pts Either \code{NULL} (default) or a set of spatial points. This can be either a \code{SpatVector} (\code{terra} package) or \code{POINTS} or \code{MULTIPOINTS} \code{sf} object (\code{sf} package). \emph{These must be in an equal-area projection!} This can also be a \code{Spatial} object (e.g., \code{SpatialPoints} or \code{SpatialPointsDataFrame}) from the \pkg{sp} package, but \emph{the \code{sp} package will be deprecated in 2023}.
#'
#' @param polys	Either \code{NULL} (default), or an object representing spatial polygons of (for example) counties in which a species is known to reside. \emph{This must be in an equal-area projection!}. This object can be either a \code{SpatVector} (\pkg{terra} package)), or \code{POLYGON}, \code{MULTIPOLYGON}, \code{LINESTRING}, or \code{MULTILINESTRING} \code{sf} object (\pkg{sf} package). This can also be a \code{Spatial} object (e.g., \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}) from the \pkg{sp} package, \emph{the \code{sp} package will be deprecated in 2023}.
#'
#' @param centerFrom Indicates how to locate the "reference" centroid used to identify points on each polygon. This is only relevant if both \code{pts} and \code{polys} are not \code{NULL}.
#' \itemize{
#' 	\item \code{'pts'}: The default is to use the centroid of \code{pts}, which finds the centroid of \code{pts}, then finds the location on the border of each polygon closest to this centroid.
#' 	\item \code{'polys'}: This option will first calculate the centroid of each polygon, then the centroid of these points, and then find the location on the border of each polygon closest to this point.
#' 	\item \code{'both'}: This option first calculates the centroid of each polygon, then finds the joint centroid of these points plus of \code{pts}, and lastly locates on the border of each polygon the point closest to this grand centroid.
#' }
#'
#' @param return Determines what is returned:
#' \itemize{
#' 		\item \code{'mcp'} (default): The minimum convex polygon
#'		\item \code{'mcpPoints'}: Points of the vertices of the minimum convex polygon
#'		\item \code{'polyPoints'}: The point on each \code{poly} polygon closest to the given center
#'	}
#'
#' @param terra If \code{TRUE} (default), the return an object of class \code{SpatVector}. Otherwise, return an object of class \code{sf}.
#
#' @details This function constructs a minimum convex polygon (MCP) from a set of spatial points and/or spatial polygons. The manner in which this is done depends on whether \code{polys} and/or \code{pts} are specified:
#' \itemize{
#' \item Only \code{pts} is supplied: The MCP is constructed directly from the points.
#' \item Only \code{polys} is supplied: The MCP is constructed from the point on each polygon closest to the centroid of the centroids of the polygons.
#' \item Both \code{pts} and \code{polys} are supplied: The MCP is constructed from the combined set of \code{pts} \emph{and} from the point on each polygon closest to the centroid of \code{pts}. By default, the function uses the centroid of the precise occurrences in step (1), but this can be changed to the centroid of the centroids of the polygons or the centroid of the points defined by the union of precise occurrence points plus the centroids of the polygons.
#' }
#'
#' The function can alternatively return the points on the vertices of the MCP, or points on the input polygons closest to the reference centroid.
#'
#' @return \code{SpatVector}, or \code{sf POLYGON} representing a minimum convex polygon.
#'
#' @seealso \code{\link{nearestEnvPoints}} for the "nearest environmental point" method, a related application for environmental space.
#'
#' @references
#' Smith, A.B., Murphy, S.J., Henderson, D., and Erickson, K.D. 2023. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth.  \emph{Global Ecology and Biogeography} In press. Open access pre-print: \doi{10.1101/2021.06.10.447988}
#'
#' @example man/examples/nearestGeogPoints_example.r
#'
#' @export
nearestGeogPoints <- function(
	pts = NULL,
	polys = NULL,
	centerFrom = 'pts',
	return = 'mcp',
	terra = TRUE
) {

	if (!is.null(polys)) {
		if (inherits(polys, c('SpatVector', 'Spatial'))) polys <- sf::st_as_sf(polys)
		polys <- sf::st_geometry(polys)
	}

	if (!is.null(pts)) {
		if (inherits(pts, c('Spatial', 'SpatVector'))) pts <- sf::st_as_sf(pts)
		pts <- sf::st_geometry(pts)
	}

	if (is.null(polys) & is.null(pts)) {
		stop('Either "polys", "pts", or both must be specified.')
	} else if (is.null(polys) & !is.null(pts)) {

		allPoints <- pts
		nearestPolyPoints <- NA

	# just polygons
	} else if (!is.null(polys) & is.null(pts)) {

		### focal centroid
		polyCents <- sf::st_centroid(polys)
		polyCents <- sf::st_union(polyCents)
		center <- sf::st_centroid(polyCents)

		### find closest points to center

		vectors <- sf::st_nearest_points(polys, center)
		nearestPolyPoints <- sf::st_cast(vectors, 'POINT')

		# remove centroid
		disjunct <- sf::st_disjoint(center, nearestPolyPoints)
		disjunct <- disjunct[[1]]
		nearestPolyPoints <- nearestPolyPoints[disjunct]

		# add centroid back in (once) if it lands in a poly
		intersect <- sf::st_intersects(center, polys)
		intersect <- length(intersect) > 0
		if (intersect) nearestPolyPoints <- c(nearestPolyPoints, center)

		nearestPolyPoints <- allPoints <- sf::st_union(nearestPolyPoints)

	} else {

		### focal centroid
		pts <- sf::st_union(pts)

		if (centerFrom == 'pts') {
			center <- sf::st_centroid(pts)
		} else if (centerFrom == 'polys') {
			centers <- sf::st_centroid(polys)
			centers <- sf::st_union(centers)
			center <- sf::st_centroid(centers)
		} else {
			centers <- sf::st_centroid(polys)
			ptsPolyCents <- c(pts, centers)
			ptsPolyCents <- sf::st_union(ptsPolyCents)
			center <- sf::st_centroid(ptsPolyCents)
		}

		### find closest points to center

		vectors <- sf::st_nearest_points(polys, center)
		nearestPolyPoints <- sf::st_cast(vectors, 'POINT')

		# remove centroid
		disjunct <- sf::st_disjoint(center, nearestPolyPoints)
		disjunct <- disjunct[[1]]
		nearestPolyPoints <- nearestPolyPoints[disjunct]

		# add centroid back in (once) if it lands in a poly
		intersect <- sf::st_intersects(center, polys)
		intersect <- length(intersect) > 0
		if (intersect) nearestPolyPoints <- c(nearestPolyPoints, center)

		nearestPolyPoints <- sf::st_union(nearestPolyPoints)
		allPoints <- sf::st_union(nearestPolyPoints, pts)

	}



	if (return %in% c('mcp', 'mcpPoints')) {
		out <- sf::st_union(allPoints)
		out <- sf::st_geometry(out)
		out <- sf::st_convex_hull(out)
		if (return == 'mcpPoints') out <- sf::st_cast(out, 'MULTIPOINT')
	} else if (return == 'polyPoints') {
		out <- nearestPolyPoints
	}
	if (terra) out <- terra::vect(out)
	out

}
