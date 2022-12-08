#' Number of points in a "spatial points" object
#'
#' Returns the number of points in a \code{sf} or \code{SpatVector} object. This is typically done using either \code{length(x)} or \code{nrow(x)}, depending on whether the object in question has rows or not. This function helps in ambiguous cases, so users need not care if \code{nrow} or \code{length} is needed.

pointCount <- function(x, byFeature = FALSE) {
	x <- sf::st_as_sf(x)
	if (byFeature) {
		.countVertices(sf::st_geometry(x))
	} else {
		sum(.countVertices(sf::st_geometry(x)))
	}
}

.countVertices <- function(x) {
	out <- if (is.list(x)) {
		sapply(sapply(x, .countVertices), sum)
	} else if (is.matrix(x)) {
		nrow(x)
	} else if (sf::st_is_empty(x)) {
		0
	} else {
		1
	}
	names(out) <- NULL
	out
}
