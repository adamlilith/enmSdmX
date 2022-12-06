#' Weighted Tjur's R2
#'
#' This function calculates Tjur's R2 metric of model discrimination accuracy. Unweighted R2 is simply the difference between the mean predicted value at presence sites and the mean predicted value at absence/background sites. The weighted version allows for differing weights between presences and between absences/contrast values (i.e., the difference between the weighted mean of predictions at presences and weighted mean predictions at absences/contrast locations).
#' @param pres Predictions at presence sites.
#' @param contrast Predictions at absence/background sites.
#' @param presWeight Weights of presence cases. The default is to assign each presence case a weight of 1.
#' @param contrastWeight Weights of absence/background cases. The default is to assign each case a weight of 1.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param ... Other arguments (unused).
#'
#' @return Numeric value.
#'
#' @references Tjur, T. 2009. Coefficients of determination in logistic regression models—A new proposal: The coefficient of discrimination. \emph{The American Statistician} 63:366–372. \doi{10.1198/tast.2009.08210}.
#'
#' @seealso \code{\link[dismo]{evaluate}}, \code{\link{evalAUC}}, \code{\link{evalMultiAUC}}, \code{\link{evalContBoyce}}, \code{\link{evalThreshold}}, \code{\link{evalThresholdStats}}, \code{\link{evalTSS}}
#'
#' @examples
#' pres <- seq(0.5, 1, by=0.1)
#' contrast <- seq(0, 1, by=0.01)
#'
#' # unweighted
#' evalTjursR2(pres, contrast)
#'
#' # weighted (weight presences with low predictions more)
#' presWeight <- c(1, 1, 1, 0.5, 0.5, 0.5)
#' evalTjursR2(pres, contrast, presWeight=presWeight)
#'
#' # weighted (weight presences with high predictions more)
#' presWeight <- c(0.5, 0.5, 0.5, 1, 1, 1)
#' evalTjursR2(pres, contrast, presWeight=presWeight)
#'
#' # weight presences and absences
#' contrastWeight <- sqrt(contrast)
#' evalTjursR2(pres, contrast, presWeight=presWeight, contrastWeight=contrastWeight)
#' @export
evalTjursR2 <- function(
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	na.rm = FALSE,
	...
) {

	# if all NAs
	if (all(is.na(pres)) | all(is.na(contrast)) | all(is.na(presWeight)) | all(is.na(contrastWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
	if (length(contrastWeight) != length(contrast)) stop('You must have the same number of absence/background predictions and absence/background weights ("contrast" and "contrastWeight").')
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(contrast, contrastWeight)
		contrast <- cleanedBg[[1]]
		contrastWeight <- cleanedBg[[2]]

	}

	# Tjur's R2
	sum(pres * presWeight) / sum(presWeight) - sum(contrast * contrastWeight) / sum(contrastWeight)
	

}
