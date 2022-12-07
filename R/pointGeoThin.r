#' Thin geographic points deterministically or randomly
#'
#' This function thins geographic points such that none have nearest neighbors closer than some user-specified distance. For a given set of points that fall within this distance, thinning can be conducted in two ways.  Both begin by first calculating all pairwise distances between points. Then, clusters of points are found based on proximity using the "single-linkage" method (i.e., based on minimum distance between groups). Then, either a deterministic or random method is used to select the retained points:
#' \itemize{
#'	\item Deterministic: For each cluster, distances between each point in the cluster and all points outside of the cluster are calculated. The point retained in each cluster is the one with the greatet minimum pairwise distance to any points in any other cluster. This point will this be maximally isolated from any other point.
#'  \item Random: For each cluster, a random point is chosen.
#' }
#'
#' @param x A "spatial points" object of class \code{SpatVector}, \code{sf}, \code{data.frame}, or \code{matrix}. If \code{x} is a \code{data.frame} or \code{matrix}, then the points will be assumed to have the WGS84 coordinate system (i.e., unprojected).
#' @param minDist Minimum distance (in meters) needed between points to retain them. Points falling closer than this distance will be candidates for being discarded.
#' @param random If \code{FALSE} (default), then use the deterministic method for thinning. If \code{TRUE}, then use the random method.
#' @param longLat This is ignored if \code{x} is a \code{Spaytvector} or \code{sf} object. However, if \code{x} is a \code{data.frame} or \code{matrix}, then this should be a character or integer vector specifiying the columns in \code{x} corresponding to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. The default is to assume that the first two columns in \code{x} represent coordinates.
#' @param ... Additional arguments. Not used.
#'
#' @return Object of class \code{x}.
#'
#' @examples
#' 
#' library(sf)
#' 
#' # lemur occurrence data
#' data(mad0)
#' data(lemurs)
#' crs <- crsGet('WGS84')
#' occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
#' ll <- c('longitude', 'latitude')
#' occs <- st_as_sf(occs, coords = ll, crs = crsGet('WGS84'))
#' 
#' # deterministically thin
#' det <- pointGeoThin(x = occs, minDist = 30000)
#' 
#' # randomly thin
#' set.seed(123)
#' rand <- pointGeoThin(x = occs, minDist = 30000, random = TRUE)
#' 
#' # map
#' par(mfrow = c(1, 2))
#' plot(st_geometry(occs), cex = 1.4, main = 'Deterministic')
#' plot(st_geometry(det), pch = 21, cex = 1.4, bg = 1:nrow(det), add = TRUE)
#' plot(st_geometry(mad0), add = TRUE)
#' 
#' plot(st_geometry(occs), cex = 1.4, main = 'Random')
#' plot(st_geometry(rand), pch = 21, cex = 1.4, bg = 1:nrow(rand), add = TRUE)
#' plot(st_geometry(mad0), add = TRUE)
#' 
#' @export
pointGeoThin <- function(
	x,
	minDist,
	random = FALSE,
	longLat = 1:2,
	...
) {

	if (!inherits(x, c('SpatVector', 'sf'))) {
		x <- sf::st_as_sf(x, coords = longLat, crs = crsGet('WGS84'))
		input <- 'data.frame'
	} else if (inherits(x, 'SpatVector')) {
		x <- sf::st_as_sf(x)
		input <- 'SpatVector'
	} else if (inherits(x, 'sf')) {
		x <- sf::st_cast(x, 'POINT')
		input <- 'sf'
	} else {
		x <- sf::st_as_sf(x)
		input <- 'unknown'
	}

	# cluster based on distances
	dists <- sf::st_distance(x)
	diag(dists) <- NA
	distsDists <- stats::as.dist(dists)
	clust <- stats::hclust(distsDists, method = 'single')

	# define groups
	groups <- stats::cutree(clust, h = minDist)
	uniqueGroups <- sort(unique(groups))

	if (length(uniqueGroups) == 1L) {
	
		warning(paste('All points are less than', minDist, 'm from one another.\nOnly the point farthest from its nearest neighbor\nwillbe returned. Try increasing "minDist".'))
		
		minDist <- apply(dists, 1, min, na.rm = TRUE)
		out <- x[which.max(minDist), ]
		
	# points in > 1 group
	} else {

		# within each group, find point that is least close to point in any other group
		if (exists('out', inherits = FALSE)) rm(out)
		for (i in seq_along(uniqueGroups)) {

			group <- uniqueGroups[i]
			inPoints <- x[groups == group, , drop=FALSE]

			if (random) {
			
				keep <- sample(1L:nrow(inPoints), 1)
				keep <- inPoints[keep, ]
			
			} else {

				outPoints <- x[groups != group, , drop=FALSE]

				distsToOtherPoints <- sf::st_distance(inPoints, outPoints)
				minDistToOtherPoints <- apply(distsToOtherPoints, 1, min)
				whichMaxDistToOtherPoints <- which.max(minDistToOtherPoints)
				
				keep <- inPoints[whichMaxDistToOtherPoints, ]
				
			}
				
			out <- if (exists('out', inherits = FALSE)) {
				rbind(out, keep)
			} else {
				keep
			}

		}
		
	}
		
	if (input == 'data.frame') {
		out <- as.data.frame(x)
		out$geometry <- NULL
	} else if (input == 'SpatVector') {
		out <- terra::vect(out)
	}
	
	out

}
