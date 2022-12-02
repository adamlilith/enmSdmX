#' Proximity-based weighting for occurrences to correct for spatial bias
#'
#' This function calculates weights for points based on proximity to other points and the distance of spatial autocorrelation.\cr \cr Weights can be used, for example, to account for spatial bias in the manner in which the points were observed. Weighting is calculated on the assumption that if two points fell exactly on top of one another, they should each have a weight of 1/2. If theree points had the exact same coordinates, then their weights should be 1/3, and so on.  Increasing distance between points should increase their weight, up to the distance at which there is no "significant" spatial autocorrelation, beyond which a point should have a weight of 1. This distance needs to be supplied by the user, as it will depend on the intended use of the weights. The distance can be calculated from "standard" metrics of spatial autocorrelation (e.g., a variogram), or on the basis of knowledge of the system (e.g., maximum dispersal distance of an organism). \cr \cr For a given point \eqn{i}, the weight is defined as \deqn{w_i = 1 / (1 + \epsilon)} where \deqn{\epsilon = \sum_{n=1}^{N}((1 - d_n)/d_sac)^\alpha} in which \eqn{N} is the total number of points closer than the maximum distance (\eqn{d_sac}) of point \eqn{i}, and \eqn{d_n} the distance between focal point \eqn{i} and point \eqn{n}.  \eqn{\alpha} is a weighting factor. By default, this is set to 1, but can be changed by the user to augment or diminish the effect that neighboring points have on the weight of a focal cell. When \eqn{\alpha} is <1, neighboring points will reduce the weight of the focal point relative to the default, and when \eqn{\alpha} is >1, they will have less effect relative to the default. When all neighboring points are at or beyond the maximum distance of spatial autocorrelation, then the focal point gets a weight \eqn{w_i} of 1. When at least neighboring one point is less than this distance away, the weight of the focal point will be >0 but <1.
#'
#' @param x	A spatial points object of class \code{SpatVector} or \code{sf}.
#' @param maxDist Maximum distance beyond which a two neighboring points are assumed to have no effect on one another for calculaton of weights.
#' @param alpha Scaling parameter (see equations above).
#'
#' @return A numeric vector of weights.
#'
#' @examples
#'
#' # lemur occurrence data
#' data(lemurs)
#' wgs84 <- getCRS('WGS84')
#' occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
#' occs <- sf::st_as_sf(occs, coords=c('longitude', 'latitude'), crs=wgs84)
#' 
#' # weights
#' maxDist <- 30000 # in meters, for this example
#' w <- pointDistWeights(occs, maxDist)
#' 
#' # plot
#' plot(sf::st_geometry(occs), cex=5 * w, main='point size ~ weight')
#' plot(st_geometry(mad0), col='gainsboro', border='gray70', add=TRUE)
#' plot(sf::st_geometry(occs), cex=5 * w, add=TRUE)
#' 
#' @export

pointDistWeights <- function(x, maxDist, alpha = 1) {

	if (inherits(x, 'SpatVector')) x <- sf::st_as_sf(x)
	
	dists <- sf::st_distance(x)
	diag(dists) <- NA
	dists <- as(dists, 'matrix')
	
	w <- rep(1, nrow(x))
	for (i in 1:nrow(x)) {
	
		neighs <- which(dists[i, ] < maxDist)
		if (length(neighs) > 0) {
		
			theseDists <- dists[i, neighs]
			epsilon <- sum((1 - (theseDists / maxDist))^alpha)
			w[i] <- 1 / (1 + epsilon)
		
		}
	
	}

	w

}
