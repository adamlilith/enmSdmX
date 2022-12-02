#' WKT string for a named coordinate reference system or a spatial object
#'
#' Retrieve the Well-Known text string (WKT2) for a coordinate reference system (CRS) by name or from a spatial object. The most common usage of the function is to return the WKT2 string using an easy-to-remember name. For example, \code{getCRS('wgs84')} returns the WKT2 string for the WGS84 datum. To get a table of strings, just use \code{getsCRS()}.
#'
#' @param x This can be any of:
#' \itemize{
#'		\item Name of CRS: Each CRS has one "long" name and at least one "short" name, which appear in the table returned by \code{getCRS()}. You can use the "long" name of the CRS, or either of the two "short" names.  Spaces, case, and dashes are ignored, but to make the codes more memorable, they are shown as having them.
#' 		\item \code{NULL} (default): This returns a table of projections with their "long" and "short" names (nearly the same as \code{data(crss)}).
#'		\item An object of class \code{SpatVector}, \code{SpatRaster}, or \code{sf}. If this is a "\code{Spat}" object, then a character vector with the CRS in WKT form is returned. If a \code{sf} is supplied, then a \code{crs} object is returned in WKT format.
#' }
#'
#' @param warn If \code{TRUE} (default), then print a warning if the name of the CRS cannot be found.
#'
#' @return A string representung WKT2 (well-known text) object or a \code{data.frame}.
#' @examples
#'
#' # NB Using cat() make the CRS display nicely for human eyes.
#' cat(getCRS('WGS 84'))
#' cat(getCRS('Mollweide'))
#' cat(getCRS('WorldClim'))
#'
#' getCRS()
#' 
#' data(mad0)
#' getCRS(mad0)
#'
#' @export

getCRS <- function(
	x = NULL,
	warn = TRUE
) {

	data('crss', envir=environment(), package='enmSdmX')
	
	# return table
	if (is.null(x)) {
		out <- crss[ , c('long', 'short1', 'short2', 'region', 'projected', 'projectionGeometry', 'datum', 'type', 'notes')]
	} else if (inherits(x, c('SpatVector', 'SpatRaster'))) {
		out <- terra::crs(x)
	} else if (inherits(x, 'sf')) {
		out <- sf::st_crs(x)
	# return WKT2
	} else {

		x <- tolower(x)
		x <- gsub(x, pattern=' ', replacement='')
		x <- gsub(x, pattern='-', replacement='')
		
		crsLong <- crss$long
		crsShort1 <- crss$short1
		crsShort2 <- crss$short2
		
		crsLong <- tolower(crsLong)
		crsShort1 <- tolower(crsShort1)
		crsShort2 <- tolower(crsShort2)
		
		crsLong <- gsub(crsLong, pattern=' ', replacement='')
		crsShort1 <- gsub(crsShort1, pattern=' ', replacement='')
		crsShort2 <- gsub(crsShort2, pattern=' ', replacement='')
		
		crsLong <- gsub(crsLong, pattern='-', replacement='')
		crsShort1 <- gsub(crsShort1, pattern='-', replacement='')
		crsShort2 <- gsub(crsShort2, pattern='-', replacement='')
		
		this <- which(x == crsLong)
		if (length(this) == 0L) this <- which(x == crsShort1)
		if (length(this) == 0L) this <- which(x == crsShort2)
		if (length(this) != 0L) {
			out <- crss$wkt2[this]
		} else {
			out <- NA
			if (warn) warning('CRS name cannot be found. It must match one of "long", "short1", or "short2" in the table retuned by data(crss).')
		}
		
	}

	out

}
