#' Scales predictors
#' 
#' Scales predictors.
#'
#' @param scale Either \code{NA} (default), or \code{TRUE} or \code{FALSE}. If \code{TRUE}, the predictors will be centered and scaled by dividing by subtracting their means then dividing by their standard deviations. The means and standard deviations will be returned in the model object under an element named "\code{scales}". For example, if you do something like \code{model <- trainGLM(data, scale=TRUE)}, then you can get the means and standard deviations using \code{model$scales$means} and \code{model$scales$sds}. If \code{FALSE}, no scaling is done. If \code{NA} (default), then the function will check to see if non-factor predictors have means ~0 and standard deviations ~1. If not, then a warning will be printed, but the function will continue to do it's operations.
#'
#' @param preds A character vector with names of the predictors in \code{data}.
#'
#' @param data A data frame.
#'
#' @noRd
.scalePredictors <- function(scale, preds, data) {

	numPreds <- length(preds)
	predIsNotFactor <- rep(FALSE, numPreds)
	means <- sds <- mags <- rep(NA_real_, numPreds)
	names(means) <- names(sds) <- names(mags) <- names(predIsNotFactor) <- preds
	for (pred in preds) {
		
		if (!inherits(data[ , pred], 'factor')) {
			means[[pred]] <- mean(data[ , pred])
			sds[[pred]] <- sd(data[ , pred])
			mags[[pred]] <- abs(max(data[ , pred]))
			predIsNotFactor[[pred]] <- TRUE
		}
	
	}
	
	if (any(predIsNotFactor)) {
	
		if (is.na(scale)) {
			mags <- mags * 1E-6
			if (any(abs(means) > mags) | any(abs(sds - 1) > 1E-6)) {
				warning('Predictors do not seem to be centered and scaled. Model may be unstable.')
			}
		} else {
			for (pred in preds[predIsNotFactor]) {
				data[ , pred] <- (data[ , pred] - means[[pred]]) / sds[[pred]]
			}
		}
	}
	
	list(mean=means, sd=sds)
	
}
