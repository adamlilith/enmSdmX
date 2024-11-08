#' Assign geographically-distinct k-folds
#'
#' This function generates geographically-distinct cross-validation folds, or "geo-folds" ("g-folds" for short). Points are grouped by proximity to one another. Folds can be forced to have at least a minimum number of points in them. Results are deterministic (i.e., the same every time for the same data). \cr \cr
#' More specifically, g-folds are created using this process:
#' \itemize{
#'	\item To start, all pairwise distances between points are calculated. These are used in a clustering algorithm to create a dendrogram of relationships by distance. The dendrogram is then "cut" so it has \code{k} groups (folds). If each fold has at least the minimum desired number of points (\code{minIn}), then the process stops and fold assignments are returned.
#'	\item However, if at least one fold has fewer than the desired number of points, a series of steps is executed.
#'	\itemize{
#'		\item First, the fold with a centroid that is farthest from all others is selected. If it has sufficient points, then the next-most distant fold is selected, and so on.
#'		\item Once a fold is identified that has fewer than the desired number of points, it is grown by adding to it the points closest to its centroid, one at a time. Each time a point is added, the fold centroid is calculated again. The fold is grown until it has the desired number of points. Call this "fold #1". From hereafter, these points are considered "assigned" and not eligible for re-assignment.
#'		\item The remaining "unassigned" points are then clustered again, but this time into \code{k - 1} folds. And again, the most-distant group found that has fewer than the desired number of points is found. This fold is then grown as before, using only unassigned points. This fold then becomes "fold #2."
#'		\item The process repeats iteratively until there are \code{k} folds assigned, each with at least the desired number of points. 
#' 	}
#' }
#' The potential downside of this approach is that the last fold is assigned the remainder of points, so will be the largest. One way to avoid gross imbalance is to select the value of \code{minIn} such that it divides the points into nearly equally-sized groups.\cr\cr
#'
#' Note that in general it is probably mathematically impossible to cluster points in 2-dimensional space into \code{k} groups, each with at least \code{minIn} points, in a manner that seems "reasonable" to the eye in all cases. In experimentation, "unreasonable" results often appear when the number of groups is high.
#'
#' @param x 		A "spatial points" object of class \code{SpatVector}, \code{sf}, \code{data.frame}, or \code{matrix}. If \code{x} is a \code{data.frame} or \code{matrix}, then the points will be assumed to have the WGS84 coordinate system (i.e., unprojected).
#' @param k 		Numeric: Number of folds to create.
#' @param minIn 	Numeric: Minimum number of points required to be in a fold.
#' @param longLat 	Character or integer vector: This is ignored if \code{x} is a \code{SpatVector} or \code{sf} object. However, if \code{x} is a \code{data.frame} or \code{matrix}, then this should be a character or integer vector specifying the columns in \code{x} corresponding to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. The default is to assume that the first two columns in \code{x} represent coordinates.
#' @param method Character: Method used by \code{\link[stats]{hclust}} to cluster points. By default, this is \code{'complete'}, but other methods may give more reasonable results, depending on the case.
#' @param ... Additional arguments (unused)
#'
#' @return A vector of integers the same length as the number of points in \code{x}. Each integer indicates which fold a point in \code{x} belongs to.
#' @seealso \code{\link{geoFoldContrast}}
#'
#' @example man/examples/geoFold_examples.r
#' 
#' @export
geoFold <- function(
	x,
	k,
	minIn = 1,
	longLat = 1:2,
	method = 'complete',
	...
) {

	### initialize folds
	n <- if (inherits(x, 'sf')) {
		countPoints(x)
	} else {
		nrow(x)
	}

	if (n < k * minIn) stop(paste('Not possible to divide points into', k, 'folds each with at least', minIn, 'points each. Reduce the value of minIn or increase k.'))

	# distances
	if (inherits(x, 'SpatVector')) x <- sf::st_as_sf(x)
	if (!inherits(x, 'sf')) x <- terra::vect(x, geom=longLat, crs=getCRS('wgs84'))
	
	dists <- sf::st_distance(x)
	dists <- methods::as(dists, 'matrix')
	diag(dists) <- NA
	dists <- stats::as.dist(dists)	

	clust <- stats::hclust(dists, method = method)

	# define folds
	folds <- stats::cutree(clust, k = k)

	folds <- .renumFolds(folds, decreasing = FALSE)
	nPerFold <- table(folds)
	nAssigned <- sum(nPerFold)

	focalFold <- 1
	while (any(nPerFold < minIn) | nAssigned < n) {

		### calculate centroids of folds
		if (exists('cents', inherits = FALSE)) rm(cents)
		for (kk in 1L:k) {
		
			xInFold <- x[folds == kk, ]
			# if (inherits(x, 'sf')) {
				xInFold <- sf::st_union(xInFold)
				cent <- sf::st_centroid(xInFold)
			# } else {
				# cent <- terra::centroids(xInFold)
			# }
			cents <- if (exists('cents', inherits = FALSE)) {
				c(cents, cent)
			} else {
				cent
			}
		
		}
		
		# find centroid with LARGEST distance to any other centroid
		# if (inherits(x, 'sf')) {
			centToCentDists <- sf::st_distance(cents)
			centToCentDists <- methods::as(centToCentDists, 'matrix')
		# } else {
			# centToCentDists <- terra::distance(cents)
		# }
		centToCentDists <- methods::as(centToCentDists, 'matrix')
		
		diag(centToCentDists) <- -Inf
		maxCentToCentDists <- apply(centToCentDists, 1, max, na.rm = TRUE)
		maxDist <- max(maxCentToCentDists)
		
		# find centroid with SECOND-LARGEST distance to any other centroid
		# centToCentDists <- sweep(centToCentDists, 2L, maxCentToCentDists)
		centToCentDists[centToCentDists == maxDist] <- -Inf
		# centToCentDists[centToCentDists == 0] <- NA
		maxCentToCentDists <- apply(centToCentDists, 1, max, na.rm = TRUE)
		maxDist <- max(maxCentToCentDists)
		# maxDistIndex <- which(maxCentToCentDists == maxDist)
		maxDistIndex <- which.max(maxCentToCentDists == maxDist)
		
		# assign points closest to this centroid to this fold
		cent <- cents[maxDistIndex]
		
		newFolds <- folds
		newFolds[newFolds != focalFold] <- NA
		# newFolds[folds == maxDistIndex] <- focalFold

		nThisFold <- sum(newFolds == focalFold, na.rm = TRUE)

		while (nThisFold < minIn) {
		
			# assign next-closest point to centroid
			xOutFold <- x[which(is.na(newFolds)), ]
			# if (inherits(x, 'sf')) {
				distToCent <- sf::st_distance(xOutFold, cent)[ , 1L]
			# } else {
				# distToCent <- terra::distance(xOutFold, cent) # !!!
			# }
			closest <- which.min(distToCent)
			newFolds[which(is.na(newFolds))[closest]] <- focalFold
			
			# re-calculate group centroid
			xInFold <- x[which(newFolds == focalFold), ]
			# if (inherits(x, 'sf')) {
				xInFold <- sf::st_union(xInFold)
				cent <- sf::st_centroid(xInFold)
			# } else {
				# cent <- terra::centroids(xInFold)
			# }
			
			nThisFold <- nThisFold + 1
			
		} # if focal fold has too few points
		
		xOutFold <- x[which(is.na(newFolds)), ]
		subFolds <- geoFold(x = xOutFold, k = k - focalFold, minIn = minIn, longLat = longLat)
		# subFolds <- subFolds + max(folds, na.rm = TRUE) + 1
		subFolds <- subFolds + focalFold
		newFolds[is.na(newFolds)] <- subFolds
		
		newFolds <- newFolds
		folds <- folds
		
		folds <- newFolds
		folds <- .renumFolds(folds)
		
		nPerFold <- table(folds)
		nAssigned <- sum(nPerFold)

		focalFold <- focalFold + 1
		
	} # if any fold has too few points
	
	folds
	
}

### renumber folds from smallest to largest or vice versa
.renumFolds <- function(folds, decreasing = FALSE) {

	# folds			vector of fold assignments
	# decreasing	if \code{FALSE} (default), renumber so smallest fold is 1, most abundant last; if \code{TRUE}, vice versa

	nPerFold <- table(folds)
	foldOrder <- order(nPerFold, decreasing = decreasing)

	uniqueFolds <- sort(unique(folds))
	uniqueFoldsSorted <- uniqueFolds[foldOrder]
	nFolds <- length(uniqueFolds)

	if (!all(foldOrder == 1:nFolds)) {

		newFolds <- folds
		for (i in seq_along(uniqueFoldsSorted)) newFolds[folds == uniqueFoldsSorted[i]] <- i
		folds <- newFolds
	
	}
	
	folds

}
