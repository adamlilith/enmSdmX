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
#' @import data.table
#' @noRd
.scalePredictors <- function(scale, preds, data) {

	# which predictors are not factors?
	predIsNotFactor <- rep(FALSE, length(preds))
	names(predIsNotFactor) <- preds
	for (pred in preds) {
		predIsNotFactor[[pred]] <- !inherits(data[ , (pred)], 'factor')
	}
	
	if (any(predIsNotFactor)) {
	
		# get centers and scales
		predsNotFactors <- preds[predIsNotFactor]
		if (inherits(data, 'data.table')) {
			.SD <- .SDcols <- NULL
			means <- unlist(data[ , lapply(.SD, mean), .SDcols=predsNotFactors])
			sds <- unlist(data[ , lapply(.SD, stats::sd), .SDcols=predsNotFactors])
			mags <- unlist(abs(data[ , lapply(.SD, max), .SDcols=predsNotFactors]))
		} else {
			means <- colMeans(data[ , predsNotFactors, drop=FALSE])
			sds <- apply(data[ , predsNotFactors, drop=FALSE], 2, stats::sd)
			mags <- abs(apply(data[ , predsNotFactors, drop=FALSE], 2, max))
		}

		# scale
		if (is.na(scale)) {
			mags <- mags * 1E-6
			if (any(abs(means) > mags) | any(abs(sds - 1) > 1E-6)) {
				warning('Predictors do not seem to be centered and scaled. Model may be unstable.')
			}
		} else {
			for (pred in predsNotFactors) {
				if (inherits(data, 'data.table')) {
					data[ , (pred)] <- (data[[pred]] - means[[pred]]) / sds[[pred]]
				} else {
					data[ , pred] <- (data[ , pred] - means[[pred]]) / sds[[pred]]
				}
			}
		}
		
	}
	
	list(data = data, scales = list(mean=means, sd=sds))
	
}
