#' Thin geographic points (mostly) deterministically
#'
#' This function thins geographic points such that none have nearest neighbors closer than some user-specified distance. The results are deterministic. \cr \cr
#' Thinning is conducted by first calculating all pairwise distances between points. Then, clusters of points are found based on proximity using the "single-linkage" method (i.e., based on minimum distance between groups). For eqach cluster, the distance between each point in the cluster and outside of the cluster is calculated.  The point retained in each cluster is the one with the largest pairwise distance to any points in any other cluster.
#'
#' @param x A "spatial points" object of class \code{SpatVector}, \code{sf}, \code{data.frame}, or \code{matrix}. If \code{x} is a \code{data.frame} or \code{matrix}, then the points will be assumed to have the WGS84 coordinate system (i.e., unprojected).
#' @param minDist Minimum distance (in meters) needed between points to retain them. Points falling closer than this distance will be candidates for being discarded.
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
#' # thin
#' thinned <- pointGeoThin(x = occs, minDist = 50000)
#' 
#' # map
#' plot(st_geometry(occs), cex = 2, main = 'Selected Points')
#' plot(st_geometry(thinned), pch = 21, cex = 2, bg = 1:nrow(out), add = TRUE)
#' plot(st_geometry(mad0), add = TRUE)
#' 
#' @export
pointGeoThin <- function(
	x,
	minDist,
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
		x <- sf::st_cast(x, 'POINTS')
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
			outPoints <- x[groups != group, , drop=FALSE]

			distsToOtherPoints <- sf::st_distance(inPoints, outPoints)
			minDistToOtherPoints <- apply(distsToOtherPoints, 1, min)
			whichMaxDistToOtherPoints <- which.max(minDistToOtherPoints)
			
			keep <- inPoints[whichMaxDistToOtherPoints, ]
			
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
