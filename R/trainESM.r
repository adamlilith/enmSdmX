#' Calibrate an ensemble of small models
#'
#' This function calibrates a set of "ensembles of small models" (ESM), which are designed for modeling species with few occurrence records. In the original formulation, each model has two covariates interacting additively. Models are calibrated using all possible combinations of covariates. By default, this function does the same, but can also include univariate models, models with two covariates plus their interaction term, and models with quadratic and corresponding linear terms. This function will \emph{only} train generalized linear models. Extending the types of algorithms is planned!
#'
#' @param data Data frame or matrix. Response variable and environmental predictors (and no other fields) for presences and non-presence sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character vector or integer vector. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data} as predictors.
#' @param scale Either \code{NA} (default), or \code{TRUE} or \code{FALSE}. If \code{TRUE}, the predictors will be centered and scaled by dividing by subtracting their means then dividing by their standard deviations. The means and standard deviations will be returned in the model object under an element named "\code{scales}". For example, if you do something like \code{model <- trainGLM(data, scale=TRUE)}, then you can get the means and standard deviations using \code{model$scales$mean} and \code{model$scales$sd}. If \code{FALSE}, no scaling is done. If \code{NA} (default), then the function will check to see if non-factor predictors have means ~0 and standard deviations ~1. If not, then a warning will be printed, but the function will continue to do its operations.
#' @param univariate,quadratic,interaction \code{TRUE} or \code{FALSE}: Whether or not to include univariate models, quadratic models, and/or models with 2-way interactions (default is \code{FALSE}).
#' @param interceptOnly If \code{TRUE}, include an intercept-only model (default is \code{FALSE}).
#' @param method Character: Name of function used to solve the GLM. For "normal" GLMs, this can be \code{'glm.fit'} (default), \code{'brglmFit'} (from the \pkg{brglm2} package), or another function.
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param family Character or function. Name of family for data error structure (see \code{\link[stats]{family}}). Default is to use the 'binomial' family.
#' @param ... Arguments to pass to \code{\link[stats]{glm}}
#' @param verbose Logical. If \code{TRUE} then display progress.
#'
#' @return A list object with several named elements:
#' \itemize{
#' 	 \item \code{models}: A list with each ESM model.
#' 	 \item \code{tuning}: A \code{data.frame} with one row per model, in the order as they appear in \code{$models}.
#' }
#'
#' @references
#' Breiner, F.T., Guisan, A., Bergamini, A., and Nobis, M.P.  2015.  Overcoming limitations of modeling rare species by using ensembles of small models.  \emph{Methods in Ecology and Evolution} 6:1210-1218.. \doi{10.1111/2041-210X.12403}
#' Lomba, A., L. Pellissier, C. Randin, J. Vicente, J. Horondo, and A. Guisan.  2010.  Overcoming the rare species modeling complex: A novel hierarchical framework applied to an Iberian endemic plant. \emph{Biological Conservation} 143:2647-2657. \doi{10.1016/j.biocon.2010.07.007}
#'
#' @seealso \code{\link[enmSdmX]{trainBRT}}, \code{\link[enmSdmX]{trainGAM}}, \code{\link[enmSdmX]{trainGLM}}, \code{\link[enmSdmX]{trainMaxEnt}}, \code{\link[enmSdmX]{trainMaxNet}}, \code{\link[enmSdmX]{trainNS}}, \code{\link[enmSdmX]{trainRF}}, \code{\link[enmSdmX]{trainByCrossValid}}
#'
#' @example man/examples/trainESM_examples.r
#' 
#' @export
trainESM <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	univariate = FALSE,
	quadratic = FALSE,
	interaction = FALSE,
	interceptOnly = FALSE,
	method = 'glm.fit',
	scale = NA,
	w = TRUE,
	family = stats::binomial(),
	...,
	verbose = FALSE
) {

	# response and predictors
	if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
	if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

	# weights and scaling
	w <- .calcWeights(w, data = data, resp = resp, family = family)
	if (is.na(scale) || scale) {
		scaleds <- .scalePredictors(scale, preds, data)
		data <- scaleds$data
		scales <- scaleds$scales
	}

	# intercept-only
	if (interceptOnly) {
		tuning <- data.frame(intercept = 1, pred1 = NA, pred2 = NA)
	}

	# bivariate
	if (!exists('tuning', inherits = FALSE)) tuning <- data.frame()
	for (i in 1:(length(preds) - 1)) {
		for (j in (i + 1):length(preds)) {
			tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds[i], pred2 = preds[j]))
		}
	}

	# univariate
	if (univariate) {
		tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds, pred2 = NA))
	}

	# interaction
	if (interaction) {

		tuning$pred3 <- NA
		for (i in 1:(length(preds) - 1)) {
			for (j in (i + 1):length(preds)) {
				tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds[i], pred2 = preds[j], pred3 = paste0(preds[i], ':', preds[j])))
			}
		}

	}

	# quadratic
	if (quadratic) {

		if (!any(names(tuning) == 'pred3')) tuning$pred3 <- NA
		tuning$pred4 <- NA

		if (univariate) tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds, pred2 = paste0('I(', preds, '^2)'), pred3 = NA, pred4 = NA))
			
		for (i in 1:(length(preds) - 1)) {
			for (j in (i + 1):length(preds)) {
				tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds[i], pred2 = preds[j], pred3 = paste0('I(', preds[i], '^2)'), pred4 = NA))
			}
		}
		for (i in 1:(length(preds) - 1)) {
			for (j in (i + 1):length(preds)) {
				tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds[i], pred2 = preds[j], pred3 = paste0('I(', preds[j], '^2)'), pred4 = NA))
			}
		}
		for (i in 1:(length(preds) - 1)) {
			for (j in (i + 1):length(preds)) {
				tuning <- rbind(tuning, data.frame(intercept = 1, pred1 = preds[i], pred2 = preds[j], pred3 = paste0('I(', preds[i], '^2)'), pred4 = paste0('I(', preds[i], '^2)')))
			}
		}
	}

	models <- list()
	tuning$model <- NA
	for (i in 1:nrow(tuning)) {
		
		form <- paste0(resp, ' ~')
		for (j in 1:(ncol(tuning) - 1)) {
			if (j == 1) {
				if (!is.na(tuning[i, j])) form <- paste(form, tuning[i, j])
			} else {
				if (!is.na(tuning[i, j])) form <- paste0(form, ' + ', tuning[i, j])
			}
		}
		tuning$model[i] <- form

		if (verbose) omnibus::say(form)

		form <- stats::as.formula(form)
		args <- list(formula = form, data = data, weights = w, method = method, family = family, ...)
		models[[length(models) + 1]] <- do.call(stats::glm, args = args)

	}

	out <- list(models = models, tuning = tuning)
	if (is.na(scale) || scale) out$scale <- scales
	out

}
