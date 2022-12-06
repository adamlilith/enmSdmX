#' Calculate the precision of a geographic coordinate
#'
#' This function calculates the imprecision of geographic coordinates due to rounded coordinate values. See \emph{Details} for an explanation of how this is calculated.
#'
#' @param x Spatial points represented as a \code{SpatVector} or \code{sf} object. Alternatively, this can also be a data frame or matrix, in which the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84).
#'
#' @param dms Logical, if \code{FALSE} (default), it is assumed that the original format in which coordinate were reported is in decimal notation. If \code{TRUE}, then it will be calculated assuming the coordinate were originally in degrees-minutes-seconds format.  If you do not know the original format, the less presumptive approach is to calculate coordinate imprecision twice with or without \code{dm = TRUE}, and use the larger of the two values.
#'
#' @param epsilon Zero or positive integer, number of digits to which to round seconds value if \code{dms} is \code{TRUE}. Ignored if \code{dms} is \code{FALSE}. This is used to accommodate inexact integer values when converting from DMS to decimal. For example, -108.932222 converted to DMS format is 108deg 55min 7.9992sec, but if \code{epsilon} = 2 then it would be converted to 108deg 55min 08sec.
#'
#' @return Numeric values (by default in units of meters).
#'
#' @details
#' For coordinates originally reported in decimal notation, coordinate imprecision is \emph{half} the distance between the two opposing corners on a bounding box whose size is based on the number of significant digits in the coordinates. This box is defined by 1) finding the maximum number of significant digits after the decimal in the longitude/latitude pair; 2) adding/subtracting 5 to the decimal place that falls just after this; and 3) calculating the distance between these points then dividing by 2. For example, if longitude is 82.37 and latitude 45.8 then the number of significant digits after the decimal place is 2 and 1, respectively so 2 is used on the assumption that latitude is measured to the nearest 100th degree. The precision is then the distance between the point pairs (82.37 - 0.05 = 82.365, 45.8 - 0.05 = 45.795) and (82.37 + 0.05 = 82.375, 45.8 + 0.05 = 45.805). \cr \cr
#'
#' For coordinates originally reported in degree-minus-second (DMS) format, the bounding box is defined by adding/subtracting 0.5 units (degrees, minutes, or seconds, depending on the smallest non-zero unit reported) from the coordinate. For example, if longitude is 90deg 00min 00sec and latitude is 37deg 37min 37sec, then the bounding box will be defined by adding/subtracting 0.5 arcsec to the coordinates.
#'
#' @examples
#'
#' # coarse-precision cases
#' long <-	c(45, 45.1, 45.1)
#' lat <-  c(45, 45.1, 45)
#' ll <- cbind(long, lat)
#' precision_m <- coordImprecision(ll)
#' cbind(ll, precision_m)
#'
#' # fine-precision cases
#' long <-	rep(45, 8)
#' lat <-  c(45, 45.1, 45.11, 45.111, 45.1111, 45.11111, 45.111111, 45.1111111)
#' ll <- cbind(long, lat)
#' precision_m <- coordImprecision(ll)
#' cbind(ll, precision_m)
#'
#' # precision varies with latitude
#' long <- rep(45, 181)
#' lat <- seq(-90, 90)
#' ll <- cbind(long, lat)
#' precision_m <- coordImprecision(ll)
#' cbind(ll, precision_m)
#' plot(lat, precision_m / 1000, xlab='Latitude', ylab='Precision (km)')
#' 
#' # dateline/polar cases
#' long <-	c(0, 180, 45, 45)
#' lat <-  c(45, 45, 90, -90)
#' ll <- cbind(long, lat)
#' precision_m <- coordImprecision(ll)
#' cbind(ll, precision_m)
#'
#' # original coordinates in degrees-minutes-seconds format
#' longDD <- c(90, 90, 90, 90, 90, 90)
#' longMM <- c(0, 0, 0, 11, 11, 0)
#' longSS <- c(0, 0, 0, 0, 52, 52)
#' latDD <- c(38, 38, 38, 38, 38, 38)
#' latMM <- c(0, 37, 37, 37, 37, 0)
#' latSS <- c(0, 0, 38, 38, 38, 0)
#' longHemis <- rep('W', 6)
#' latHemis <- rep('N', 6)
#' longDec <- dmsToDecimal(longDD, longMM, longSS, longHemis)
#' latDec <- dmsToDecimal(latDD, latMM, latSS, latHemis)
#' decimal <- cbind(longDec, latDec)
#' (decImp <- coordImprecision(decimal))
#' (dmsImp <- coordImprecision(decimal, dms=TRUE))
#'
#' # What if we do not know if coordinates were originally reported in
#' # decimal or degrees-minutes-seconds format? Most conservative option
#' # is to use maximum:
#' pmax(decImp, dmsImp)
#'
#' if (FALSE) {
#'   # known error when longitude is negative and latitude is -90
#'   long <- -45
#'   lat <- -90
#'   ll <- cbind(long, lat)
#'   coordImprecision(ll)
#' }
#'
#' @export
coordImprecision <- function(
	x,
	dms = FALSE,
	epsilon = 2
) {

	if (inherits(x, 'SpatVector')) {
		crs <- terra::crs(x)
		x <- terra::geom(x)[ , c('x', 'y')]
	} else if (inherits(x, 'sf')) {
		crs <- sf::st_crs(x)
		x <- sf::st_coordinates(x)[ , c('x', 'y')]
	} else if (inherits(x, c('matrix', 'data.frame'))) {
		crs <- getCRS('wgs84')
	} else {
		stop('Argument "x" must be a SpatVector or sf object, or a matrix or data.frame.')
	}
	
	## correct for original conversion from DDMMSS format
	if (dms) {
	
		longDms <- decimalToDms(x[ , 1L])
		latDms <- decimalToDms(x[ , 2L])
	
		longDms[ , 'ss'] <- round(longDms[ , 'ss'], epsilon)
		latDms[ , 'ss'] <- round(latDms[ , 'ss'], epsilon)
		
		# locate bounding box corners for each coordinate pair
		longDmsPlus <- longDmsMinus <- longDms
		latDmsPlus <- latDmsMinus <- latDms
		
		for (i in 1L:nrow(x)) {
		
			# recorded to nearest degree
			if (longDms[i, 'mm'] == 0 & longDms[i, 'ss'] == 0 &
				latDms[i, 'mm'] == 0 & latDms[i, 'ss'] == 0) {
				
				longDmsPlus[i, 'dd'] <- longDmsPlus[i, 'dd'] + 0.5
				latDmsPlus[i, 'dd'] <- latDmsPlus[i, 'dd'] + 0.5
		
				longDmsMinus[i, 'dd'] <- longDmsMinus[i, 'dd'] - 0.5
				latDmsMinus[i, 'dd'] <- latDmsMinus[i, 'dd'] - 0.5
		
			# recorded to nearest arcmin
			} else if (longDms[i, 'ss'] == 0 & latDms[i, 'ss'] == 0) {
			
				longDmsPlus[i, 'mm'] <- longDmsPlus[i, 'mm'] + 0.5
				latDmsPlus[i, 'mm'] <- latDmsPlus[i, 'mm'] + 0.5
			
				longDmsMinus[i, 'mm'] <- longDmsMinus[i, 'mm'] - 0.5
				latDmsMinus[i, 'mm'] <- latDmsMinus[i, 'mm'] - 0.5
			
			# recorded to nearest arcsec
			} else if (trunc(longDms[i, 'ss']) == longDms[i, 'ss'] & trunc(latDms[i, 'ss']) == latDms[i, 'ss']) {
			
				longDmsPlus[i, 'ss'] <- longDmsPlus[i, 'ss'] + 0.5
				latDmsPlus[i, 'ss'] <- latDmsPlus[i, 'ss'] + 0.5
			
				longDmsMinus[i, 'ss'] <- longDmsMinus[i, 'ss'] - 0.5
				latDmsMinus[i, 'ss'] <- latDmsMinus[i, 'ss'] - 0.5
			
			# recorded using decimal seconds
			} else {
			
				longDigits <- omnibus::countDecDigits(longDms[i, 'ss'])
				latDigits <- omnibus::countDecDigits(latDms[i, 'ss'])
				
				digits <- max(longDigits, latDigits)
			
				delta <- 5 * 10^-(1 + digits)
				
				longDmsPlus[i, 'ss'] <- longDmsPlus[i, 'ss'] + delta
				latDmsPlus[i, 'ss'] <- latDmsPlus[i, 'ss'] + delta
			
				longDmsMinus[i, 'ss'] <- longDmsMinus[i, 'ss'] - delta
				latDmsMinus[i, 'ss'] <- latDmsMinus[i, 'ss'] - delta
			
			}
		
		} # next coordinate pair
		
		longHemis <- ifelse(x[ , 1] >= 0, 'E', 'W')
		latHemis <- ifelse(x[ , 2] >= 0, 'N', 'S')
		
		longDecPlus <- dmsToDecimal(longDmsPlus[ , 'dd'], longDmsPlus[ , 'mm'], longDmsPlus[ , 'ss'], hemis=longHemis)
		latDecPlus <- dmsToDecimal(latDmsPlus[ , 'dd'], latDmsPlus[ , 'mm'], latDmsPlus[ , 'ss'], hemis=latHemis)

		longDecMinus <- dmsToDecimal(longDmsMinus[ , 'dd'], longDmsMinus[ , 'mm'], longDmsMinus[ , 'ss'], hemis=longHemis)
		latDecMinus <- dmsToDecimal(latDmsMinus[ , 'dd'], latDmsMinus[ , 'mm'], latDmsMinus[ , 'ss'], hemis=latHemis)

		plusPlus <- cbind(longDecPlus, latDecPlus)
		minusMinus <- cbind(longDecMinus, latDecMinus)
		
	# if decimal format
	} else {

		plusPlus <- minusMinus <- x
	
		for (i in 1L:nrow(x)) {

			## get coordinates and add deltas
			xy <- c(x[i, ])
			
			# number of significant digits based on plain count of digits
			digitsCount <- omnibus::countDecDigits(xy)
			digits <- max(digitsCount)

			# get bounding box of imprecision
			delta <- 5 * 10^-(1 + digits)

			plusPlus[i, ] <- plusPlus[i, ] + delta
			minusMinus[i, ] <- minusMinus[i, ] - delta

		}
		
	}

	
	## calculate coordinate uncertainty due to imprecision
	out <- rep(NA, nrow(x))
	
	for (i in 1L:nrow(x)) {

		## correct datelines/polar cases

		# crosses dateline but not at pole
		if (plusPlus[i, 1] > 180 & (minusMinus[i, 2] > -90 & plusPlus[i, 2] < 90)) {

			# shift east 10 deg
			plusPlus[i, 1] <- plusPlus[i, 1] - 10
			minusMinus[i, 1] <- minusMinus[i, 1] - 10

		# crosses date line but not at pole
		} else if (minusMinus[i, 1] < -180 & (minusMinus[i, 2] > -90 | plusPlus[i, 2] < 90)) {

			# shift east 10 deg
			plusPlus[i, 1] <- plusPlus[i, 1] + 10
			minusMinus[i, 1] <- minusMinus[i, 1] + 10

		# crosses dateline and over a pole
		} else if (plusPlus[i, 1] > 180 & (minusMinus[i, 2] < -90 | plusPlus[i, 2] > 90)) {

			delta <- abs(plusPlus[i, 2] - minusMinus[i, 2])
			
			# north pole: shift east 10 deg and "down"
			if (plusPlus[i, 2] > 90) {

				plusPlus[i, ] <- c(plusPlus[i, 1] - 10, 90)
				minusMinus[i, ] <- c(minusMinus[i, 1] - 10, 90 - delta)

			# south pole: shift east 10 deg and "up"
			} else if (minusMinus[i, 2] < -90) {

				plusPlus[i, ] <- c(plusPlus[i, 1] - 10, -90)
				minusMinus[i, ] <- c(minusMinus[i, 1] - 10, -90 + delta)

			}

		# crosses dateline and over a pole
		} else if (plusPlus[i, 1] < -180 & (minusMinus[i, 2] < -90 | plusPlus[i, 2] > 90)) {

			delta <- abs(plusPlus[i, 2] - minusMinus[i, 2])
			
			# north pole: shift east 10 deg and "down"
			if (plusPlus[i, 2] > 90) {

				plusPlus[i, ] <- c(plusPlus[i, 1] + 10, 90)
				minusMinus[i, ] <- c(minusMinus[i, 1] + 10, 90 - delta)

			# south pole: shift east 10 deg and "up"
			} else if (minusMinus[i, 2] < -90) {

				plusPlus[i, ] <- c(plusPlus[i, 1] + 10, -90)
				minusMinus[i, ] <- c(minusMinus[i, 1] + 10, -90 + delta)

			}

		# at north pole but does not cross dateline
		} else if (plusPlus[i, 2] > 90) {

			delta <- abs(plusPlus[i, 2] - minusMinus[i, 2])
			plusPlus[i, ] <- c(plusPlus[i, 1], 90)
			minusMinus[i, ] <- c(minusMinus[i, 1], 90 - delta)

		# at south pole but does not cross dateline
		} else if (minusMinus[i, 2] < -90) {

			delta <- abs(plusPlus[i, 2] - minusMinus[i, 2])
			plusPlus[i, ] <- c(plusPlus[i, 1], -90)
			minusMinus[i, ] <- c(minusMinus[i, 1], -90 + delta)

		}

	}
	
	# calculate precision
	plusPlus <- sf::st_as_sf(as.data.frame(plusPlus), coords=1:2, crs=crs)
	minusMinus <- sf::st_as_sf(as.data.frame(minusMinus), coords=1:2, crs=crs)
	
	out <- 0.5 * sf::st_distance(plusPlus, minusMinus, by_element = TRUE)
	out

}
