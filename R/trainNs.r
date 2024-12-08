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
#' @param maxTerms Maximum number of terms to be used in any model, not including the intercept (default is 8).
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param removeInvalid Logical. If \code{TRUE} (default), remove models that either did not converge or have parameter estimates near the boundaries (usually negative or positive infinity).
#' @param failIfInvalid Logical. If \code{TRUE} (default), and the "full" model either does not converge or has parameters near the boundary, then the function will fail. If \code{FALSE}, then return \code{NULL} in this case.
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
	removeInvalid = TRUE,
	failIfInvalid = TRUE,
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {

	if (FALSE) {
	
		resp <- 'presBg'

		scale <- TRUE
		df <- 1:4
		interaction <- TRUE
		interceptOnly <- TRUE
		method <- 'glm.fit'
		presPerTermFinal <- 10
		maxTerms <- 8
		w <- TRUE
		family <- 'binomial'
		removeInvalid <- TRUE
		failIfInvalid <- TRUE
		out <- 'model'
		cores <- 1
		verbose <- TRUE
	
	}

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

		n <- if (family %in% c('binomial', 'quasibinomial')) {
			sum(data[ , resp, drop=TRUE])
		} else {
			nrow(data)
		}

	### model selection
	###################

		### make list of candidate model terms

		terms <- character()
		factors <- sapply(data[ , preds, drop = FALSE], is.factor)
		names(factors) <- preds
		componentTerms <- character()

		# univariate terms
		for (thisPred in preds) {
		
			if (factors[thisPred]) {
				terms <- c(terms, thisPred)
				componentTerms <- c(componentTerms, thisPred)
			} else {

				for (thisDf in df) {
				
					term <- paste0('splines::ns(', thisPred, ', df=', thisDf, ')')
					terms <- c(terms, term)
					componentTerms <- c(componentTerms, thisPred)

				}

			}
		
		}

		numTerms <- rep(1, length(terms))

		# interaction terms
		if (interaction & length(preds) > 1L & n >= 2 * presPerTermFinal) {
		
			predCombos <- utils::combn(preds, m = 2, simplify = FALSE)

			for (i in seq_along(predCombos)) {

				pred1 <- predCombos[[i]][1]
				pred2 <- predCombos[[i]][2]

				if (factors[[pred1]] & !factors[[pred2]]) {
					newTerm <- paste0(pred1, '* splines::ns(', pred2, ', df=', df, ')')
				} else if (!factors[[pred1]] & factors[[pred2]]) {
					newTerm <- paste0(pred2, '* splines::ns(', pred1, ', df=', df, ')')
				} else {
					newTerm <- paste0('splines::ns(', pred1, ' * ', pred2, ', df=', df, ')')
				}

				componentTerms <- c(componentTerms, rep(paste(pred1, pred2, collapse = ' '), length(df)))
				numTerms <- c(numTerms, rep(2, length(df)))
				terms <- c(terms, newTerm)
					
			}
			
		} # if more than one term

		### select best set of terms for full model
		###########################################

		forms <- terms
		# if (interceptOnly) forms <- c(forms, '1')
		
		if (verbose) omnibus::say('Evaluating simple models with each candidate term...')
		work <- foreach::foreach(
			i = seq_along(forms),
			.options.multicore = mcOptions,
			.combine = 'rbind',
			.inorder = FALSE,
			.export = c('.trainNsWorker')
		) %makeWork% {
			.trainNsWorker(
				i = i,
				forms = forms,
				numTerms = numTerms,
				componentTerms = componentTerms,
				data = data,
				resp = resp,
				family = family,
				method = method,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = FALSE,
				...
			)
		}

		if (removeInvalid) {

			bads <- which(!work$converged | work$boundary)
			if (length(bads) > 0) {
				work <- work[-bads, , drop = FALSE]
				
				if (nrow(work) == 0) {
					msg <- 'No single-term models converged or all models had parameter estimates near the boundary.'
					if (failIfInvalid) {
						stop(msg)
					} else {
						warning(msg)
						return(NULL)
					}
				}
			}
		}

		work <- work[order(work$AICc), , drop = FALSE]
		rownames(work) <- NULL

		if (verbose) {
			omnibus::say('Term-by-term evaluation:', pre = 1)
			print(work[ , c('term', 'converged', 'boundary', 'AICc')])
			utils::flush.console()
		}

	### model construction
	######################

		topTerms <- work$term[1L]
		totalTerms <- work$numTerms[1L]

		termsSoFar <- work$componentTerms[1L]
			
		i <- 2L
		while (i <= nrow(work) && (n >= (totalTerms + numTerms[i]) * presPerTermFinal & totalTerms < maxTerms + numTerms[i])) {

			if (!(work$componentTerms[i] %in% termsSoFar)) {
				topTerms <- c(topTerms, work$term[i])
				totalTerms <- totalTerms + work$numTerms[i]
				termsSoFar <- c(termsSoFar, work$componentTerms[i])
			}

			i <- i + 1L

		}

		formsNA <- list()
		for (i in seq_along(topTerms)) {
			formsNA[[i]] <- c(topTerms[i], NA_character_)
		}

		formGrid <- expand.grid(formsNA, stringsAsFactors = FALSE)

		forms <- character()
		for (i in 1:nrow(formGrid)) {
		
			terms <- unlist(formGrid[i, , drop = TRUE])
			terms <- terms[!is.na(terms)]
			if (length(terms) > 0) {
				form <- paste(terms, collapse = ' + ')
				forms[i] <- form
			}
		
		}


		if (verbose) omnibus::say('Evaluating all possible models constructed from full model...', pre = 1)

		wantModels <- any(c('models', 'model') %in% out)
		combine <- if (wantModels) { 'c' } else { 'rbind' }

		work <- foreach::foreach(
			i = seq_along(forms),
			.options.multicore = mcOptions,
			.combine = combine,
			.inorder = FALSE,
			.export = c('.trainNsWorker')
		) %makeWork% {
			.trainNsWorker(
				i = i,
				forms = forms,
				numTerms = numTerms,
				componentTerms = componentTerms,
				data = data,
				resp = resp,
				family = family,
				method = method,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = wantModels,
				...
			)
		}

		if (wantModels) models <- list()
		tuning <- data.frame()
		for (i in seq_along(work)) {
		
			if (wantModels) {
				models[[i]] <- work[[i]]$model
				models[[i]]$scales <- scales
			}

			tuning <- rbind(
				tuning,
				data.frame(
					model = work[[i]]$formula,
					converged = work[[i]]$model$converged,
					boundary = work[[i]]$model$boundary,
					AICc = work[[i]]$AICc
				)
			)
		
		}

		if (removeInvalid) {

			bads <- which(!tuning$converged | tuning$boundary)
			if (length(bads) > 0) {
				
				tuning <- work[-bads, , drop = FALSE]
				
				if (nrow(work) == 0) {
					msg <- 'No single-term models converged or all models had parameter estimates near the boundary.'
					if (failIfInvalid) {
						stop(msg)
					} else {
						warning(msg)
						return(NULL)
					}
				}

				if (wantModels) models <- models[-bads]

			}

		}

		bestOrder <- order(tuning$AICc)
		tuning <- tuning[bestOrder, , drop = FALSE]
		rownames(tuning) <- NULL
		if (wantModels) {
			models <- models[bestOrder]
			model <- models[[1L]]
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

### train NS
############

.trainNsWorker <- function(
	i,
	forms, # formulae (without LHS and intercept)
	numTerms, # vector of number of predictors in each term
	componentTerms, # character vector of predictors in each term
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

	# form <- forms[i]
	# if (insertIntercept) {
	# 	form <- if (form == '') {
	# 		'1'
	# 	} else {
	# 		paste('1', form, sep=' + ')
	# 	}
	# }
	thisForm <- paste0(resp, ' ~ ', forms[i])
	thisForm <- stats::as.formula(thisForm)

	model <- suppressWarnings(stats::glm(
		formula = thisForm,
		data = data,
		family = family,
		method = method,
		weights = w,
		...
	))
	
	AICc <- AICcmodavg::AICc(model)
	
	# out
	out <- if (modelOut) {
		
		list(
			list(
				model = model,
				formula = forms[i],
				numTerms = numTerms[i],
				componentTerms = componentTerms[i],
				converged = model$converged,
				boundary = model$boundary,
				AICc = AICc
			)
		)
		
	} else {
	
		data.frame(
			term = forms[i],
			componentTerms = componentTerms[i],
			numTerms = numTerms[i],
			converged = model$converged,
			boundary = model$boundary,
			AICc = AICc
		)
	
	}
	out

}
