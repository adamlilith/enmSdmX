#' Plot response curves for one or more models
#'
#' This function creates plots of response curves for one or more models. Response curves show how model predictions change as a particular predictor is changed, while holding all other predictors constant. The function can plot response curves for a single model or for multiple models, and can plot response curves for multiple predictors across multiple models.
#'
#' @param models Either a model object (like \code{glm}, \code{gam}, etc. -- any model produced by a \code{trainXYZ()} function), or a \code{list} of such model objects. If a single model object is passed, the function will create response curves for that model. If a list of model objects is passed, the function will create response curves for each model in the list. The models must be able to be passed to \code{\link{predictEnmSdm}}.
#'
#' @param env A \code{data.frame} containing environmental data. Typically, this data represents the range of environmental conditions over which predictions will be made. The data must contain one column for each predictor variable in the model(s) for which response curves will be plotted. Values do not need to be in a particular order.
#'
#' @param ref Either \code{NULL} (default), a \code{data.frame}, or a named vector of number ic values. These are used for setting the value of non-focal variables.
#' \itemize{
#'   \item \code{NULL}: The function will use the mean value of each predictor in \code{env} as the "constant" value for each predictor.
#'   \item \code{data.frame}: Mean values across each predictor are used as their "constant" values.
#'   \item \code{named vector}: The names of the vector should correspond to the names of the predictors in \code{env}. The values of the vector should be the "constant" values for each predictor.
#' }
#'
#' @param vars Either \code{NULL} (default) or a character vector of predictor names. If \code{NULL}, response curves will be plotted for all predictors in \code{env}. If a character vector, response curves will be plotted for the predictors specified in \code{vars}.
#'
#' @param modelNames Either \code{NULL} or a character vector of model names. If \code{NULL} (default), the names of the models will be set to 'Model 1', 'Model 2', etc. If a character vector, the names of the models will be set to the values in \code{modelNames}. The length of \code{modelNames} must be equal to the number of models passed to \code{models}.
#'
#' @param constantFx Function used to calculate the constant value for each predictor. The default is \code{mean}. The function must take a numeric vector as input and return a single numeric value. The function is applied to each column of \code{ref} to "calculate" the constant value for each predictor.
#'
#' @param combine Logical. If \code{TRUE} (default), graphs will be combined using \code{\link[cowplot]{plot_grid}}. If \code{FALSE}, the output will be a \code{list} object, with one plot per element.
#'
#' @param ncol Number of columns to use when combining plots. If \code{combine = FALSE}, this argument is ignored. Default is \code{NULL}, in which case the function will use the square root of the number of plots as the number of columns.
#'
#' @param n Number of points to use for each variable in the response curve. Default is 200.
#'
#' @returns Either a \pkg{ggplot} \code{grob} object or a \code{list} of \pkg{ggplot} \code{grob} objects.
#'
#' @example man/examples/trainXYZ_examples.r
#'
#' @export 
responseCurves <- function(
	models,
	env,
	ref = NULL,
	vars = NULL,
	modelNames = NULL,
	constantFx = mean,
	combine = TRUE,
	ncol = NULL,
	n = 200
) {

	xDUMMYDUMMY <- NULL
	yDUMMYDUMMY <- NULL
	modelDUMMYDUMMY <- NULL

	if (is.null(vars)) vars <- names(env)

	if (!all(vars %in% names(env))) stop('The names of the variables in `vars` must match the names of the variables in `env`.')

	if (!is.null(ref)) if (!all(vars %in% names(ref))) stop('The names of the variables in `vars` must match the names of the variables in `ref`.')

	# calculate constant values
	if (is.null(ref)) ref <- env
	constants <- apply(ref[ , vars, drop = FALSE], 2, constantFx, na.rm = TRUE)

	# models
	if (!is.list(models)) models <- list(models)

	# create data frame with increasing values for each variable
	minsEnv <- sapply(env[ , vars], min, na.rm = TRUE)
	maxsEnv <- sapply(env[ , vars], max, na.rm = TRUE)

	minsRef <- sapply(ref[ , vars], min, na.rm = TRUE)
	maxsRef <- sapply(ref[ , vars], max, na.rm = TRUE)

	mins <- pmin(minsEnv, minsRef, na.rm = TRUE)
	maxs <- pmax(maxsEnv, maxsRef, na.rm = TRUE)

	gradients <- data.frame(n = 1:n)
	for (variable in vars) {

		# get range of values for the variable
		vals <- seq(mins[variable], maxs[variable], length.out = n)

		# create data frame
		gradients[[variable]] <- seq(mins[variable], maxs[variable], length.out = n)

	}

	# create plots
	out <- list()
	for (i in seq_along(vars)) {
	
		variable <- vars[i]

		# create data frame with focal variable varying but all others constant
		constantGrad <- gradients
		varsSansFocal <- vars[vars != variable]
		constantGrad[ , varsSansFocal] <- constants[[varsSansFocal]]

		# predict for each model
		df <- data.frame()
		for (m in seq_along(models)) {

			model <- models[[m]]

			# get predictions
			preds <- predictEnmSdm(model, newdata = constantGrad)

			constantGrad$yDUMMYDUMMY <- preds
			modelName <- ifelse(is.null(modelNames), paste('Model', m), modelNames[m])
			constantGrad$modelDUMMYDUMMY <- modelName

			df <- rbind(df, constantGrad)

		}

		names(df)[names(df) == variable] <- 'xDUMMYDUMMY'

		# plot
		graph <- ggplot2::ggplot() +
			ggplot2::geom_line(
				data = df,
				ggplot2::aes(
					x = xDUMMYDUMMY,
					y = yDUMMYDUMMY,
					color = modelDUMMYDUMMY
					# linetype = modelDUMMYDUMMY
				)
			) +
			# ggplot2::scale_color_manual(name = 'Model') +
			# ggplot2::scale_linetype_manual(name = 'Model') +
			ggplot2::ylim(0, 1) +
			ggplot2::labs(x = variable, y = 'Prediction', color = 'Model') +
			ggplot2::ggtitle(variable)

		# if (rug) {
		
		# 	graph <- graph + 
		# 		ggplot2::geom_rug(
		# 			data = env,
		# 			ggplot2::aes(x = ggplot2::.data[[variable]], y = ggplot2::.data[[variable]]),
		# 			sides = 't'
		# 		)
		
		# }

		out[[i]] <- graph
	
	} # next variable

	if (combine) {
		if (is.null(ncol)) ncol <- ceiling(sqrt(length(out)))
		out <- cowplot::plot_grid(plotlist = out, ncol = ncol)
	}
	out

}
