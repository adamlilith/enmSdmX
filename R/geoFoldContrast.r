#' Assign geographically-distinct k-folds to background/absence sites
#'
#' This function generates geographically-distinct cross-validation folds, or "geo-folds" of background or absence sites (i.e., "contrast" sites). Each contrast site is assigned to a fold based on the fold of the presence site that is closest. Typically, this function is run after \code{\link{geoFold}} is run to assign presences to folds.
#'
#' @param contrast	A "spatial points" object representing contrast sites:
#' \itemize{
#'	\itemize{
#'		\item A \code{SpatVector} or \code{sf} vector with points
#'  	\item A \code{data.frame} or \code{matrix}: Points will be assumed to have the WGS84 coordinate system (i.e., unprojected), and \code{contrastLongLat} should denote the columns with coordinates.
#' }
#' }
#' @param pres 		A "spatial points" object representing presence sites:
#'	\itemize{
#'		\item A \code{SpatVector} or \code{sf} vector with points
#'  	\item A \code{data.frame} or \code{matrix}: Points will be assumed to have the WGS84 coordinate system (i.e., unprojected), and \code{presLongLat} should denote the columns with coordinates.
#' }
#' @param presFolds		Numeric vector: These provide the folds to which \code{pres} are assigned. There must be one value per point in \code{pres}.
#' @param contrastLongLat,presLongLat 	Character or integer vector: A character or integer vector specifying the columns in \code{contrast} and \code{pres} corresponding to longitude and latitude (in that order). The default is to assume that the first two columns in \code{contrast} represent coordinates. These are ignored if \code{contrast} or \code{pres} are a \code{SpatVector} or an \code{sf} object.
#' @param ... Additional arguments (unused)
#'
#' @return A vector of integers the same length as the number of points in \code{contrast}. Each integer indicates which fold a point in \code{contrast} belongs to.
#'
#' @seealso \code{\link{geoFold}}
#'
#' @example man/examples/geoFold_examples.r
#' 
#' @export
geoFoldContrast <- function(
	contrast,
	pres,
	presFolds,
	contrastLongLat = 1:2,
	presLongLat = 1:2,
	...
) {

	if (!inherits(contrast, c('SpatVector', 'sf'))) contrast <- sf::st_as_sf(contrast, coords=contrastLongLat, crs=getCRS('wgs84'))
	if (!inherits(pres, c('SpatVector', 'sf'))) pres <- sf::st_as_sf(pres, coords=presLongLat, crs=getCRS('wgs84'))
	
	if (!inherits(contrast, 'sf')) contrast <- sf::st_as_sf(contrast)
	if (!inherits(pres, 'sf')) pres <- sf::st_as_sf(pres)
	
	dists <- sf::st_distance(contrast, pres)

	closest <- apply(dists, 1L, which.min)

	presFolds[closest]
	
}
