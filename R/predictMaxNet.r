#' Predictions from a MaxNet model
#'
#' This function is the same as the \code{predict} function in the \pkg{maxnet} package, except that:
#' \itemize{
#'		\item	If the input is a data frame, the output is a vector as output (not a single-column matrix);
#'		\item	If the input is a \code{SpatRaster}, the output is a \code{SpatRaster};
#'		\item	The default output is on the cloglog scale;
#'		\item   The function can be explicitly called (versus doing, say, \code{maxnet:::predict.maxnet}, which does not work even when that would be really useful...).
#' }
#'
#' @param model		Object of class \code{maxnet}.
#' @param newdata	Object of class \code{data.frame} or \code{SpatRaster} (\pkg{terra} package).
#' @param clamp		If \code{TRUE} (default), predict outside the range of training data by 'clamping' values to the last value.
#' @param type		One of:
#' \itemize{
#'		\item		\code{cloglog} (default): Predictions are on a complementary log-log scale.
#'		\item		\code{logistic}: Predictions are on a logistic scale (and thus technically the same to several decimal places as predictions from MaxEnt <=3.3.3k, except for differences in default features).
#'		\item		\code{link}: Predictions are on the scale of the predictors.
#'		\item		\code{exponential}: Predictions are on an exponential ('raw') scale.
#'	}
#' @param ... Other arguments (unused).
#' @return Numeric vector or \code{SpatRaster}
#'
#' @example man/examples/trainXYZ_examples.R
#'
#' @seealso \code{\link[raster]{predict}} from the \pkg{raster} package, \code{\link[terra]{predict}} from the \pkg{terra} package, and \code{\link[maxnet]{maxnet}} (see the \code{predict} function therein)
#'
#' @export
predictMaxNet <- function(
	model,
	newdata,
	clamp = TRUE,
	type = 'cloglog',
	...
) {
	
	if (inherits(newdata, 'SpatRaster')) {

		nd <- terra::as.data.frame(newdata, na.rm=FALSE)
		notNas <- which(stats::complete.cases(nd))
		nd <- nd[notNas, , drop=FALSE]
		preds <- do.call(predictMaxNet, args = list(model = model, newdata = nd, clamp = clamp, type = type, ...))

		if (inherits(newdata, 'SpatRaster')) {
			out <- newdata[[1L]]
			out[] <- NA
			out <- setValueByCell(out, preds, cell=notNas, format='raster')
		} else {
			out <- rep(NA, nrow(newdata))
			out[notNas] <- preds
		}

	} else {

		if (clamp) {
			for (v in intersect(names(model$varmax), names(newdata))) {
				newdata[ , v] <- pmin(pmax(newdata[ , v], model$varmin[v]), model$varmax[v])
			}
		}
		terms <- sub('hinge\\((.*)\\):(.*):(.*)$', 'maxnet:::hingeval(\\1,\\2,\\3)', names(model$betas))
		terms <- sub('categorical\\((.*)\\):(.*)$', 'maxnet:::categoricalval(\\1,\'\\2\')', terms)
		terms <- sub('thresholds\\((.*)\\):(.*)$', 'maxnet:::thresholdval(\\1,\\2)', terms)
		f <- stats::formula(paste('~', paste(terms, collapse=' + '), '-1'))
		mm <- stats::model.matrix(f, data.frame(newdata))
		if (clamp) mm <- t(pmin(pmax(t(mm), model$featuremins[names(model$betas)]), 
			model$featuremaxs[names(model$betas)]))
		link <- (mm %*% model$betas) + model$alpha
		if (type == 'cloglog') {
			out <- 1 - exp(0 - exp(model$entropy + link))
		} else if (type == 'logistic') {
			out <- 1 / (1 + exp(-model$entropy - link))
		} else if (type == 'exponential') {
			out <- exp(link)
		} else if (type == 'link') {
			out <- link
		} else {
			stop('Argument "type" must be one of "cloglog", "logistic", "exponential", or "link".')
		}
		out <- out[ , 1L]
	}
	out
}
