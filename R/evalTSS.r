#' Weighted True Skill Statistic (TSS)
#'
#' This function calculates the True Skill Statistic (TSS).
#' @param pres Numeric vector. Predicted values at test presences
#' @param contrast Numeric vector. Predicted values at background/absence sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param thresholds Numeric vector. Thresholds at which to calculate the sum of sensitivity and specificity. The default evaluates all values from 0 to 1 in steps of 0.01.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param ... Other arguments (unused).
#' @return Numeric value.
#' @details This function calculates the maximum value of the True Skill Statistic (i.e., across all thresholds, the values that maximizes sensitivity plus specificity).
#' @references  See Allouche, O., Tsoar, A., and Kadmon, R. 2006. Assessing the accuracy of species distribution models: Prevalence, kappa and the true skill statistic (TSS). \emph{Journal of Applied Ecology} 43:1223-1232. \doi{10.1111/j.1365-2664.2006.01214.x}
#'
#' @seealso \code{\link[dismo]{evaluate}}, \code{\link{evalAUC}}, \code{\link{evalMultiAUC}}, \code{\link{evalContBoyce}}, \code{\link{evalThreshold}}, \code{\link{evalThresholdStats}}, \code{\link{evalTjursR2}}
#'
#' @examples
#' set.seed(123)
#' 
#' # set of bad and good predictions at presences
#' bad <- runif(30)^2
#' good <- runif(30)^0.1
#' hist(good, breaks=seq(0, 1, by=0.1), border='green', main='Presences')
#' hist(bad, breaks=seq(0, 1, by=0.1), border='red', add=TRUE)
#' pres <- c(bad, good)
#' contrast <- runif(1000)
#' evalTSS(pres, contrast)
#' 
#' # upweight bad predictions
#' presWeight <- c(rep(1, 100), rep(0.1, 30))
#' evalTSS(pres, contrast, presWeight=presWeight)
#' 
#' # upweight good predictions
#' presWeight <- c(rep(0.1, 100), rep(1, 30))
#' evalTSS(pres, contrast, presWeight=presWeight)
#' 
#' e <- dismo::evaluate(pres, contrast)
#' max(e@TPR + e@TNR) - 1
#' 
#' # Why different values from dismo's evaluate() function?
#' # Because dismo's function uses thresholds based on presence/non-presence
#' # values, whereas evalTSS uses equall-spaced thresholds.
#' head(e@t)
#'
#' @export
evalTSS <- function(
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	thresholds = seq(0, 1, by=0.001),
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

		cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
		contrast <- cleanedContrast[[1]]
		contrastWeight <- cleanedContrast[[2]]

	}

	# stats
	sumPresWeights <- sum(presWeight)
	sumContrastWeights <- sum(contrastWeight)
	
	numPres <- length(pres)
	numContrast <- length(contrast)
	
	# TSS
	tss <- rep(NA, length(thresholds))
	
	# for each threshold
	for (i in seq_along(thresholds)) {
		
		thisThresh <- thresholds[i]
	
		# which presences/contrast sites are correctly predicted at this threshold
		whichCorrectPres <- which(pres >= thisThresh)
		whichCorrectContrast <- which(contrast < thisThresh)
		
		numCorrectPres <- length(whichCorrectPres)
		numCorrectContrast <- length(whichCorrectContrast)
		
		anyCorrectPres <- (numCorrectPres > 0)
		anyCorrectContrast <- (numCorrectContrast > 0)
		
		# weights of correctly predicted predictions
		correctPresWeights <- if (anyCorrectPres) {
			sum(presWeight[whichCorrectPres])
		} else {
			0
		}
		
		correctContrastWeights <- if (anyCorrectContrast) {
			sum(contrastWeight[whichCorrectContrast])
		} else {
			0
		}
		
		# true positive/negative rates
		tpr <- correctPresWeights / sumPresWeights
		tnr <- correctContrastWeights / sumContrastWeights
		tss[i] <- tpr + tnr - 1
	
	}
	
	thresholdMaxTss <- thresholds[which.max(tss)]
	tss <- max(tss)
	attr(tss, 'thresholdMaxTss') <- thresholdMaxTss
	tss

}
