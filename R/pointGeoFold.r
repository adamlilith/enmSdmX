#' Assign geographically-distinct k-folds
#'
#' This function assigns geographically-divided k-folds ("g-folds") to a set of points based on clustering by proximity. Folds are created by grouping nearby points together. Users can force folds to have at least a minimum number of points by:
#' \itemize{
#'	\item Swapping points between folds. This ensures the desired number of folds is produced, but can sometimes result in folds that are "squashed" between others.
#' 	\item Collapsing smaller folds into nearby larger ones. This reduces the total number of folds, sometimes even down to a single fold.
#' 	\item A hybrid approach in which folds can be collapsed until a threshold number of folds is reached, then the swapping rountine takes over. This helps ensure at least a given number of folds is produced, but is less likely to yield "squashed" folds.
#' }
#' Regardless of which method is used, the results are deterministic.
#'
#' @param x 		A "spatial points" object of class \code{SpatVector}, \code{sf}, \code{data.frame}, or \code{matrix}. If \code{x} is a \code{data.frame} or \code{matrix}, then the points will be assumed to have the WGS84 coordinate system (i.e., unprojected).
#' @param k 		Number of folds to create.
#' @param minIn 	Positive integer or \code{NULL}. Minimum number of sites required to be in a fold. If left \code{NULL} (default), it is possible to have just one site in a fold.
#' @param collapseFolds		Either:
#' \itemize{
#' 	\item If \code{TRUE} and there are fewer than \code{minIn} points in any fold, collapse folds with fewer than \code{minIn} points into neighboring folds to create larger (but fewer folds). This can sometimes create a single fold.  If this is undesirable, then try setting \code{kMin} to a number >1 but < \code{k}. This will turn the swapping routine on once there are \code{kMin} folds left.
#' 	\item if \code{FALSE} (default) and there are fewer than \code{minIn} points in any fold, then force a swapping routine so that there will be at least \code{minIn} points in each fold. This can sometimes create folds that are "squished" between others.
#' }
#' @param kMin		Minimum number of folds desired. This is used only if \code{collapseFolds} is \code{FALSE}. If >1, then the routine will collapse folds into one another to achieve the desired number of points per fold, but then switch to the swapping routine once there are \code{kMin} folds remaining.  The default is 1, meaning that the function would return a single fold if this was the only way to create a fold with sufficient number of points.
#' @param longLat 	This is ignored if \code{x} is a \code{Spaytvector} or \code{sf} object. However, if \code{x} is a \code{data.frame} or \code{matrix}, then this should be a character or integer vector specifiying the columns in \code{x} corresponding to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. The default is to assume that the first two columns in \code{x} represent coordinates.
#' @param ... Additional arguments. Not used.
#'
#' @return A vector of integers the same length as the dimension of \code{x}. Each integer indicates which fold a point in \code{x} belongs to.
#'
#' @example man/examples/geoFold_examples.r
#' 
#' @export
pointGeoFold <- function(
	x,
	k,
	kMin = 1,
	minIn = 1,
	collapseFolds = FALSE,
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

	n <- nrow(x)
	if (n < k * minIn) stop(paste('Not possible to divide points into', k, 'folds each with at least', minIn, 'points each. Reduce the value of minIn or increase k.'))

	par(mfrow=c(2, 3))

	# cluster based on distances
	allDists <- sf::st_distance(x)
	allDists <- as(allDists, 'matrix')
	diag(allDists) <- NA
	dists <- stats::as.dist(allDists)
	clust <- stats::hclust(dists, method = 'single')

	# define groups
	folds <- stats::cutree(clust, k = k)
	nPerFold <- table(folds)
	
	plot(st_geometry(occs), col=folds, pch=16)
	
	if (any(nPerFold < minIn)) {
	
		if (!collapseFolds) okToSwap <- rep(TRUE, n)

		# if at least one group too small
		while (any(nPerFold < minIn)) {

			### combine folds
			#################
			if (collapseFolds & k > kMin) {
			
				k <- k - 1
				folds <- stats::cutree(clust, k = k)
				if (k == kMin) collapseFolds <- FALSE

			### swap points between folds
			#############################
			} else {

				# define ingroup (group to grow)
				smallestIngroup <- which.min(nPerFold)
				ingroupIndices <- which(folds == smallestIngroup)

				okToSwap[ingroupIndices] <- FALSE
				
				ingroupFold <- unique(folds[ingroupIndices])

				# agglomerate ingroup with closest points in outgroups

				nIngroup <- length(ingroupIndices)
				nNeeded <- minIn - nIngroup
				
				distsToIngroup <- allDists[ , ingroupIndices, drop = FALSE]
				distsToIngroup[!okToSwap, ] <- NA

				distToIngroup <- apply(distsToIngroup, 1, min)
				closestToIngroup <- rank(distToIngroup)
			
				assignToIngroup <- which(closestToIngroup <= nNeeded)
				folds[assignToIngroup] <- ingroupFold
				okToSwap[assignToIngroup] <- FALSE

				folds <- .renumSeq(folds)
				
				nPerFold <- table(folds)
				nFolds <- length(unique(folds))
				
				# if we lost folds, recreate them from unassigned points
				if (nFolds < k) {
				
					outgroupDists <- allDists[okToSwap, okToSwap]
					outgroupDists <- stats::as.dist(outgroupDists)
					clust <- stats::hclust(outgroupDists, method = 'single')
					newFolds <- stats::cutree(clust, k = k - nFolds + 1)
					newFolds <- newFolds + 1E6
					
					newFolds <- .insert(x = folds[!okToSwap], into = newFolds, at = which(!okToSwap), warn = FALSE)
					folds <- .renumSeq(newFolds)
				
				}
				
			} # if swapping
			
			nPerFold <- table(folds)
			plot(st_geometry(occs), col=folds, pch=16)
			
		} # while too few points in at least one fold
		
	} # if too few points in a fold

	folds

}

# insert values into a vector
.insert <- function(x, into, at, warn = TRUE) {

	if (length(x) > length(at)) stop('Length of x is longer than the number of indices.')
	if (any(at > length(at) + length(into))) stop(paste('At least one index is too high. The new vector will be', length(x) + length(into), 'elements long.'))

	out <- rep(NA, length(into) + length(at))
	if (length(x) < length(at)) {
	
		if (warn) warning('Length of x is shorter than the length of at. Recycling x.')
		x <- rep(x, length.out = length(at))
	
	}
	
	out[at] <- x
	out[is.na(out)] <- into
	
	out

}

### renumber folds sequentially
.renumFolds <- function(folds) {

	folds <- folds - min(folds) + 1
	uniqueFolds <- sort(unique(folds))
	nFolds <- length(uniqueFolds)
	if (!all(uniqueFolds == 1:nFolds)) {

		newFolds <- folds
		for (i in seq_along(uniqueFolds)) newFolds[folds == uniqueFolds[i]] <- i
		folds <- newFolds
	
	}
	
	folds

}

.renumSeq <- function(x) {

	# x <- x - min(x, na.rm = TRUE) + 1
	xUnique <- sort(unique(x))
	n <- length(xUnique)
	if (!all(xUnique == 1:n)) {

		xNew <- rep(NA, length(x))
		for (i in seq_along(xUnique)) xNew[x == xUnique[i]] <- i
		x <- xNew
	
	}
	
	x

}
