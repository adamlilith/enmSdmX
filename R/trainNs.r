#' Calibrate a natural splines model
#'
#' This function constructs a natural-spline model by evaluating all possible models given the available predictors and constraints. "Constraints" in this case include the degrees of freedom for a spline, whether or not interaction terms are included, minimum number of presence sites per model term, and maximum number of terms to include in the model. Its output is any or all of: a table with AICc for all evaluated models; all models evaluated; and/or the single model with the lowest AICc.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response variable. The default is to use the first column in \code{data} as the response.
#' @param preds Character vector or integer vector. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param scale Either \code{NA} (default), or \code{TRUE} or \code{FALSE}. If \code{TRUE}, the predictors will be centered and scaled by dividing by subtracting their means then dividing by their standard deviations. The means and standard deviations will be returned in the model object under an element named "\code{scales}". For example, if you do something like \code{model <- trainGLM(data, scale=TRUE)}, then you can get the means and standard deviations using \code{model$scales$means} and \code{model$scales$sds}. If \code{FALSE}, no scaling is done. If \code{NA} (default), then the function will check to see if non-factor predictors have means ~0 and standard deviations ~1. If not, then a warning will be printed, but the function will continue to do it's operations.
#' @param method Character, name of function used to solve. This can be \code{'glm.fit'} (default), \code{'brglmFit'} (from the \pkg{brglm2} package), or another function.
#' @param df A vector of integers > 0 or \code{NULL}. Sets flexibility of model fit. See documentation for \code{\link[splines]{ns}}.
#' @param interaction If \code{TRUE} (default), include two-way interaction terms.
#' @param interceptOnly If \code{TRUE} (default) and model selection is enabled, then include an intercept-only model.
#' @param presPerTermFinal Minimum number of presence sites per term in initial starting model.
#' @param maxTerms Maximum number of terms to be used in any model, not including the intercept (default is 8). Used only if \code{construct} is \code{TRUE}.
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by AICc.
#' }
#' @param cores Number of cores to use. Default is 1. If you have issues when \code{cores} > 1, please see the \code{\link{troubleshooting_parallel_operations}} guide.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device. Default is \code{FALSE}.
#' @param ... Arguments to send to \code{\link[stats]{glm}}.
#'
#' @returns The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these. If \code{scale} is \code{TRUE}, any model object will also have an element named \code{$scale}, which contains the means and standard deviations for predictors that are not factors.
#'
#' @seealso \code{\link[splines]{ns}}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export
trainNS <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	scale = NA,
	df = 1:4,
	interaction = TRUE,
	interceptOnly = TRUE,
	method = 'glm.fit',
	presPerTermFinal = 10,
	maxTerms = 8,
	w = TRUE,
	family = 'binomial',
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########

		# degrees of freedom
		if (is.null(df)) df <- 'NULL'

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		w <- .calcWeights(w, data = data, resp = resp, family = family)
		
		if (is.na(scale) || scale) {
			scaleds <- .scalePredictors(scale, preds, data)
			data <- scaleds$data
			scales <- scaleds$scales
		}

	### parallelization
	###################
			
		cores <- min(cores, parallel::detectCores(logical = FALSE))

		if (cores > 1L) {

			`%makeWork%` <- foreach::`%dopar%`
			# cl <- parallel::makeCluster(cores, setup_strategy = 'sequential')
			cl <- parallel::makeCluster(cores)
			parallel::clusterEvalQ(cl, requireNamespace('parallel', quietly=TRUE))
			doParallel::registerDoParallel(cl)
			on.exit(parallel::stopCluster(cl), add=TRUE)
			
			# `%makeWork%` <- doRNG::`%dorng%`
			# doFuture::registerDoFuture()
			# future::plan(future::multisession(workers = cores))
			# on.exit(future:::ClusterRegistry('stop'), add=TRUE)

		} else {
			`%makeWork%` <- foreach::`%do%`
		}

		paths <- .libPaths() # need to pass this to avoid "object '.doSnowGlobals' not found" error!!!
		mcOptions <- list(preschedule = TRUE, set.seed = TRUE, silent = verbose)

	### make list of candidate model terms
	######################################

		n <- if (family %in% c('binomial', 'quasibinomial')) {
			sum(data[ , resp, drop=TRUE])
		} else {
			nrow(data)
		}

		### create vector of terms
		terms <- character()
		for (thisPred in preds) {
		
			for (thisDf in df) {
			
				term <- if (!inherits(data[ , thisPred], 'factor')) {
					paste0('splines::ns(', thisPred, ', df=', thisDf, ')')
				} else {
					thisPred
				}
				
				terms <- c(terms, term)
			}
		
		}

		# interaction terms
		if (interaction & length(preds) > 1L & n >= 2 * presPerTermFinal) {
		
			predCombos <- utils::combn(preds, m = 2, simplify = FALSE)

			for (i in seq_along(predCombos)) {

				pred1 <- predCombos[[i]][1]
				pred2 <- predCombos[[i]][2]

				terms <- c(terms, paste0('splines::ns(', pred1, ' * ', pred2, ', df=', df, ')'))
					
			}
			
		} # if more than one term
		
		terms <- sort(terms)

	### model selection
	###################

		### make table of all possible terms
		
		# single-predictor terms
		candidates <- list(x = paste0('splines::ns(', preds[1L], ', df=', df, ')'))
		candidates$x <- c(candidates$x, NA)
		
		if (length(preds) > 1L) {
			for (i in 2L:length(preds)) {
				candidates <- c(
					candidates,
					list(x = c(paste0('splines::ns(', preds[i], ', df=', df, ')'), NA))
				)
			}
		
		}
		names(candidates) <- paste0('pred', seq_along(preds))
		
		# add interaction terms
		if (interaction & length(preds) > 1L & n >= 2 * presPerTermFinal) {
		
			for (countPred1 in 1L:(length(preds) - 1L)) { # for each predictor test two-variable terms

				pred1 <- preds[countPred1]

				for (countPred2 in 2L:length(preds)) { # for each second predictor test two-variable terms

					pred2 <- preds[countPred2]
					candidates$x <- c(paste0('splines::ns(', preds[countPred1], ' * ', preds[countPred2], ', df=', df, ')'), NA)
					names(candidates)[length(candidates)] <- paste0('pred', length(candidates))
					
				} # next second term
				
			} # next first term
			
		} # if more than one term

		candidates <- expand.grid(candidates, stringsAsFactors = FALSE)
		
		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[numTerms > 0, , drop=FALSE]

		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[numTerms <= maxTerms, , drop=FALSE]

		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[n >= numTerms * presPerTermFinal, , drop=FALSE]
		
		forms <- rep(NA, nrow(candidates))
		for (i in 1L:nrow(candidates)) {
		
			terms <- which(!is.na(unlist(candidates[i , ])))
			terms <- unlist(candidates[i, terms])
			if (length(terms) > 1L) terms <- paste(terms, collapse = ' + ')
			forms[i] <- terms
		
		}
		
		forms <- paste0('1 + ', forms)
		if (interceptOnly) forms <- c(forms, '1')
		
		work <- foreach::foreach(
			i = seq_along(forms),
			.options.multicore = mcOptions,
			.combine = 'c',
			.inorder = FALSE,
			.export = c('.trainNsWorker')
		) %makeWork% {
			.trainNsWorker(
				i = i,
				forms = forms,
				data = data,
				resp = resp,
				family = family,
				method = method,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = ('models' %in% out | 'model' %in% out),
				...
			)
		}
	
		# if (cores > 1L) parallel::stopCluster(cl)

		# tuning table
		tuning <- data.frame(
			model = work[[1L]]$formula,
			AICc = work[[1L]]$AICc
		)
		
		if (length(work) > 1L) {
			for (i in 2L:length(work)) {
				
				tuning <- rbind(
					tuning,
					data.frame(
						model = work[[i]]$formula,
						AICc = work[[i]]$AICc
					)
				)
				
			}
		}
	
		if (!is.na(scale)) {
			if (scale) for (i in seq_along(work)) work[[i]]$model$scale <- scales
		}
	
		bestOrder <- order(tuning$AICc)
		if ('model' %in% out) model <- work[[bestOrder[1L]]]$model
		if ('models' %in% out) {
			models <- list()
			models[[1]] <- work[[1L]]$model
			for (i in 2L:length(work)) models[[i]] <- work[[i]]$model
			models <- models[bestOrder]
		}
		tuning <- tuning[bestOrder, , drop = FALSE]
		rownames(tuning) <- NULL
	
		if (verbose) {
		
			omnibus::say('Model selection:', level=2)
			print(tuning)
			utils::flush.console()
		
		}

	### return
	##########

	if (length(out) > 1L) {
		output <- list()
		if ('models' %in% out) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuning
		output
	} else if ('models' %in% out) {
		models
	} else if ('model' %in% out) {
		model
	} else if ('tuning' %in% out) {
		tuning
	}
		
}

#################
### train GAM ###
#################

.trainNsWorker <- function(
	i,
	forms, # formulae (without LHS and maybe without intercept)
	data,
	resp,
	family,
	method,
	w,
	paths,
	insertIntercept, # if TRUE, add "1 +" to the RHS side
	modelOut, # if TRUE, return model *and* data frame with model
	...
) {

	 # need to call this to avoid "object '.doSnowGlobals' not found" error!!!
	.libPaths(paths)

	form <- forms[i]
	if (insertIntercept) {
		form <- if (form == '') {
			'1'
		} else {
			paste('1', form, sep=' + ')
		}
	}
	thisForm <- paste0(resp, ' ~ ', form)

	model <- stats::glm(
		formula = stats::as.formula(thisForm),
		data = data,
		family = family,
		method = method,
		weights = w,
		...
	)
	
	AICc <- AICcmodavg::AICc(model)
	
	# out
	out <- if (modelOut) {
		
		list(
			list(
				model = model,
				formula = form,
				AICc = AICc
			)
		)
		
	} else {
	
		data.frame(
			formula = form,
			AICc = AICc
		)
	
	}
		
	out

}
