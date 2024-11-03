#' Calibrate a generalized additive model (GAM)
#'
#' This function constructs a generalized additive model. By default, the model is constructed in a two-stage process.  First, the "construct" phase generates a series of simple models with univariate and bivariate interaction terms. These simple models are then ranked based on their AICc. Second, the "select" phase creates a "full" model from the simple models such that there is at least \code{presPerTermInitial} presences (if the response is binary) or data rows (if not) for each smooth term to be estimated (not counting the intercept). Finally, it selects the best model using AICc from all possible subsets of this "full" model. Its output is any or all of: a table with AICc for all evaluated models; all models evaluated in the "selection" phase; and/or the single model with the lowest AICc.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response variable. The default is to use the first column in \code{data} as the response.
#' @param preds Character vector or integer vector. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}).
#' @param gamma Initial penalty to degrees of freedom to use (larger ==> smoother fits).
#' @param scale A numeric value indicating the "scale" parameter (see argument \code{scale} in \code{\link[mgcv]{gam}}). The default is 0 (which allows a single smoother for Poisson and binomial error families and unknown scale for all others.)
#' @param construct If \code{TRUE} (default), then construct the model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{maxTerms}.
#' @param select If \code{TRUE} (default), then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if \code{construct} is \code{TRUE}).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only if \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is \code{TRUE}.
#' @param maxTerms Maximum number of terms to be used in any model, not including the intercept (default is 8). Used only if \code{construct} is \code{TRUE}.
#' @param interaction Character or \code{NULL}. Type of interaction term to use (\code{te}, \code{ts}, \code{s}, etc.). See \code{?te} (for example) for help on any one of these. If \code{NULL}, then interactions are not used.
#' @param interceptOnly If \code{TRUE} (default) and model selection is enabled, then include an intercept-only model.
#' @param smoothingBasis Character. Indicates the type of smoothing basis. The default is \code{'cs'} (cubic splines), but see \code{\link[mgcv]{smooth.terms}} for other options. This is the value of argument \code{bs} in a \code{\link[mgcv]{s}} function.
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by AICc.
#' }
#' @param cores Integer >= 1. Number of cores to use when calculating multiple models. Default is 1. If you have issues when \code{cores} > 1, please see the \code{\link{troubleshooting_parallel_operations}} guide.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device.
#' @param ... Extra arguments (not used).
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these.
#'
#' @seealso \code{\link[mgcv]{gam}}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export


trainGAM <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	gamma = 1,
	scale = 0,
	smoothingBasis = 'cs',
	interaction = 'te',
	interceptOnly = TRUE,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 10,
	maxTerms = 8,
	w = TRUE,
	family = 'binomial',
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {

	### setup
	#########

		# ellipses <- list(...)

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		w <- .calcWeights(w, data = data, resp = resp, family = family)

	### parallelization
	###################
			
		cores <- if (!construct) {
			1L
		} else {
			min(cores, parallel::detectCores(logical = FALSE))
		}

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
		
		# single-predictor terms
		for (thisPred in preds) {
		
			if (!inherits(data[ , thisPred], 'factor')) {
				terms <- c(terms, paste0('s(', thisPred, ', bs=\'', smoothingBasis, '\')'))
			} else {
				thisPred
			}
		}
		
		# interaction terms
		if (!is.null(interaction) & length(preds) > 1L & n >= 2 * presPerTermInitial) {
		
			predCombos <- utils::combn(preds, m = 2, simplify = FALSE)

			for (i in seq_along(predCombos)) {
				
				pred1 <- predCombos[[i]][1]
				pred2 <- predCombos[[i]][2]

				# create term
				if (!inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

					terms <- c(terms, paste0(interaction, '(', pred1, ', ', pred2, ', bs=\'', smoothingBasis, '\')'))

				} else if (inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

					terms <- c(terms, paste0(interaction, '(', pred2, ', by=', pred1, ', bs=\'', smoothingBasis,'\')'))

				} else if (!inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

					terms <- c(terms, paste0(interaction, '(', pred1, ', by=', pred2, ', bs=\'', smoothingBasis, '\')'))

				} else if (inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

					terms <- c(terms, paste0(pred1, ' * ', pred2))

				}
				
			}
							
			# for (countPred1 in 1L:(length(preds)-1L)) { # for each predictor test two-variable terms

				# pred1 <- preds[countPred1]

				# for (countPred2 in 2:length(preds)) { # for each second predictor test two-variable terms

					# pred2 <- preds[countPred2]

					# # create term
					# if (!inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

						# terms <- c(terms, paste0(interaction, '(', pred1, ', ', pred2, ', bs=\'', smoothingBasis, '\')'))

					# } else if (inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

						# terms <- c(terms, paste0(interaction, '(', pred2, ', by=', pred1, ', bs=\'', smoothingBasis,'\')'))

					# } else if (!inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

						# terms <- c(terms, paste0(interaction, '(', pred1, ', by=', pred2, ', bs=\'', smoothingBasis, '\')'))

					# } else if (inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

						# terms <- c(terms, paste0(pred1, ' * ', pred2))

					# }
					
				# } # next second term
				
			# } # next first term
			
		} # if more than one term and wanting interactions
			
	## term-by-term model construction
	##################################
	if (construct) {

		assess <- foreach::foreach(
			i = seq_along(terms),
			.options.multicore = mcOptions,
			.combine = 'rbind',
			.inorder = FALSE,
			.export = c('.trainGamWorker')
		) %makeWork% {
			.trainGamWorker(
				i = i,
				forms = terms,
				data = data,
				resp = resp,
				family = family,
				gamma = gamma,
				scale = scale,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = FALSE,
				...
			)
				
		}
		
		assess <- assess[order(assess$AICc), , drop=FALSE]
		rownames(assess) <- NULL

		if (verbose) {
			omnibus::say('Term-by-term evaluation', level=2)
			print(assess)
			utils::flush.console()
		}

		### no selection
		################
		
		# just return model with "best" terms without further selection
		
		if (!select) {
			
			requiredPresPerTerm <- presPerTermFinal * 1L:nrow(assess)
			
			these <- which(requiredPresPerTerm <= n)
			theseTerms <- assess$formula[these]
			form <- paste(theseTerms, collapse=' + ')
			thisForm <- paste0(resp, ' ~ 1 + ', form)
			
			model <- mgcv::gam(
				formula = stats::as.formula(thisForm),
				family = family,
				data = data,
				method = 'ML',
				optimizer = c('outer', 'newton'),
				scale = scale,
				select = TRUE,
				gamma = gamma,
				weights = w,
				...
			)
			
			AICc <- AICcmodavg::AICc(model)
			
			tuning <- data.frame(
				model = form,
				AICc = AICc
			)

			models <- NULL

			if (verbose) {

				omnibus::say('Final model (construction from best terms, but no selection)', level=2)
				print(summary(model))
				utils::flush.console()
				
			}

		### model selection
		###################
		
		# select best model from all possible subsets of "full" model with best terms
		
		} else {

			# make formulae for all possible models
			candidates <- list(x = c(TRUE, FALSE))
			if (length(terms) > 1L) for (i in 2L:length(terms)) candidates <- c(candidates, list(x = c(TRUE, FALSE)))
			candidates <- expand.grid(candidates, stringsAsFactors = FALSE)
			names(candidates) <- terms
			candidates <- candidates[rowSums(candidates) != 0, , drop=FALSE]

			numTerms <- rowSums(candidates)
			requiredPresPerTerm <- numTerms * presPerTermFinal
			candidates <- candidates[requiredPresPerTerm <= n, , drop = FALSE]
			
			numTerms <- rowSums(candidates)
			candidates <- candidates[numTerms <= maxTerms, , drop=FALSE]

			forms <- character()
			for (i in 1L:nrow(candidates)) {
				thisForm <- paste0('1 + ', paste(terms[unlist(candidates[i, , drop = TRUE])], collapse = ' + '))
				forms <- c(forms, thisForm)
			}
			if (interceptOnly) forms <- c(forms, '1')

			work <- foreach::foreach(
				i = seq_along(forms),
				.options.multicore = mcOptions,
				.combine = 'c',
				.inorder = FALSE,
				.export = c('.trainGamWorker')
			) %makeWork% {
				.trainGamWorker(
					i = i,
					forms = forms,
					data = data,
					resp = resp,
					family = family,
					gamma = gamma,
					scale = scale,
					w = w,
					insertIntercept = FALSE,
					paths = paths,
					modelOut = ('models' %in% out | 'model' %in% out),
					...
				)
			}
			
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
			
				omnibus::say('Model selection', level=2)
				print(tuning)
				utils::flush.console()
			
			}
		
		} # if selecting best model from subsets of "full" model

	### if not constructing model term-by-term (selection not possible)
	###################################################################
	} else {
	
		form <- paste(terms, collapse = ' + ')
		thisForm <- paste0(resp, ' ~ 1 + ', form)
	
		model <- mgcv::gam(
			formula = stats::as.formula(thisForm),
			family = family,
			data = data,
			method = 'ML',
			optimizer = c('outer', 'newton'),
			scale = scale,
			select = TRUE,
			gamma = gamma,
			weights = w,
			...
		)
		
		AICc <- AICcmodavg::AICc(model)
		
		tuning <- data.frame(
			model = form,
			AICc = AICc
		)
		
		models <- NULL
		
		if (select) warning('Model selection is not performed when argument "construct" is FALSE.')
		
		if (verbose) {
		
			omnibus::say('Model (no construction or selection)', level=2)
			print(summary(model))
			utils::flush.console()
		
		}
		
	} # if not constructing model term-by-term

	# if (cores > 1L) parallel::stopCluster(cl)

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

.trainGamWorker <- function(
	i,
	forms, # formulae (without LHS and maybe without intercept)
	data,
	resp,
	family,
	gamma,
	scale,
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

	model <- mgcv::gam(
		formula = stats::as.formula(thisForm),
		family = family,
		data = data,
		method = 'ML',
		optimizer = c('outer', 'newton'),
		scale = scale,
		select = TRUE,
		gamma = gamma,
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
