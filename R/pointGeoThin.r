#' Thin geographic points (mostly) deterministically
#'
#' This function thins geographic points such that none have nearest neighbors closer than some user-specified distance. The results are almost deterministic (see Details).
#'
#' @param x A "points" object of class \code{SpatVector}, \code{sf}, \code{data.frame}, or \code{matrix}. If \code{x} is a \code{data.frame} or \code{matrix}, then the points will be assumed to have the WGS84 coordinate system (i.e., unprojected).
#' @param minDist Minimum distance (in meters) needed between points to retain them. Points falling closer than this distance will be candidates for being discarded.
#' @param longLat This is ignored if \code{x} is a \code{Spaytvector} or \code{sf} object. However, if \code{x} is a \code{data.frame} or \code{matrix}, then this should be a character or integer vector specifiying the columns in \code{x} corresponding to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. The default is to assume that the first two columns in \code{x} represent coordinates.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param ... Additional arguments. Not used.
#' @details
#' The procedure for removing points is as follows:
#' \itemize{
#' 	\item Calculate all pairwise distances between points.
#' 	\item Find point(s) with largest number of neighbors (<\code{minDist} away). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the points with the closest neighbor within \code{minDist}. If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the point that is closest to the centroid of all non-removed points. If just one such point exists, remove it, but if there is more than one...
#' 	\item Of these find the point that has the smallest median distance to all points (even if > \code{minDist}). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these randomly select a point and remove it.
#' 	\item Repeat.
#' }
#' Thus the results are deterministic up to the last tie-breaking step.
#' @return Object of class \code{x}.
#' @seealso \code{\link{geoThinApprox}}
#' @examples

# lemur occurrence data

library(sf)

data(lemurs)
crs <- crsGet('WGS84')
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
occs <- sf::st_as_sf(occs, coords=c('longitude', 'latitude'), crs=crs)

dists <- sf::st_distance(occs)
diag(dists) <- NA
distsDists <- stats::as.dist(dists)

clust <- stats::hclust(distsDists, method = 'single')

minDist <- 10000 # in meters


plot(cutree(clust, h = minDist))



#' @export

geoThin <- function(
	x,
	minDist,
	longLat = 1:2,
	verbose = FALSE,
	...
) {

	dists <- 

}
