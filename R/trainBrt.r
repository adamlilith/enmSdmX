#' Calibrate a boosted regression tree (generalized boosting machine) model
#'
#' This function calibrates a boosted regression tree (or gradient boosting machine) model, and is a wrapper for \code{\link[dismo]{gbm.step}}. The function uses a grid search to assess the best combination of learning rate, tree depth, and bag fraction based on cross-validated deviance. If a particular combination of parameters leads to an unconverged model, the script attempts again using slightly different parameters. Its output is any or all of: a table with deviance of evaluated models; all evaluated models; and/or the single model with the lowest deviance.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response variable. The default is to use the first column in \code{data} as the response.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param family Character. Name of error family.  See \code{\link[dismo]{gbm.step}}.
#' @param learningRate Numeric. Learning rate at which model learns from successive trees (Elith et al. 2008 recommend 0.0001 to 0.1).
#' @param treeComplexity Positive integer. Tree complexity: depth of branches in a single tree (1 to 16).
#' @param bagFraction Numeric in the range [0, 1]. Bag fraction: proportion of data used for training in cross-validation (Elith et al. 2008 recommend 0.5 to 0.7).
#' @param minTrees Positive integer. Minimum number of trees to be scored as a "usable" model (Elith et al. 2008 recommend at least 1000). Default is 1000.
#' @param maxTrees Positive integer. Maximum number of trees in model set (same as parameter \code{max.trees} in \code{\link[dismo]{gbm.step}}).
#' @param tries Integer > 0. Number of times to try to train a model with a particular set of tuning parameters. The function will stop training the first time a model converges (usually on the first attempt). Non-convergence seems to be related to the number of trees tried in each step.  So if non-convergence occurs then the function automatically increases the number of trees in the step size until \code{tries} is reached.
#' @param tryBy Character list. A list that contains one or more of \code{'learningRate'}, \code{'treeComplexity'}, \code{numTrees}, and/or \code{'stepSize'}. If a given combination of \code{learningRate}, \code{treeComplexity}, \code{numTrees}, \code{stepSize}, and \code{bagFraction} do not allow model convergence then then the function tries again but with alterations to any of the arguments named in \code{tryBy}:
#' * \code{learningRate}: Decrease the learning rate by a factor of 10.
#' * \code{treeComplexity}: Randomly increase/decrease tree complexity by 1 (minimum of 1).
#' * \code{maxTrees}: Increase number of trees by 20%.
#' * \code{stepSize}: Increase step size (argument \code{n.trees} in \code{gbm.step()}) by 50%.
#' If \code{tryBy} is NULL then the function attempts to train the model with the same parameters up to \code{tries} times.
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param anyway Logical. If \code{FALSE} (default), it is possible for no models to be returned if none converge and/or none had a number of trees is >= \code{minTrees}). If \code{TRUE} then all models are returned but with a warning.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest deviance.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest deviance.
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by deviance.
#' }
#' @param cores Integer >= 1. Number of cores to use when calculating multiple models. Default is 1.
#' @param parallelType Either \code{'doParallel'} (default) or \code{'doSNOW'}. Issues with parallelization might be solved by trying the non-default option.
#' @param verbose Logical. If \code{TRUE} display progress.
#' @param ... Arguments to pass to \code{\link[dismo]{gbm.step}}.
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of two or more of these.
#'
#' @seealso \code{\link[dismo]{gbm.step}}
#'
#' @references
#' Elith, J., J.R. Leathwick, & T. Hastie. 2008. A working guide to boosted regression trees. \emph{Journal of Animal Ecology} 77:802-813. \doi{10.1111/j.1365-2656.2008.01390.x}
#'
#' @example man/examples/trainXYZ_examples.R
#'
#' @export
trainBRT <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	learningRate = c(0.0001, 0.001, 0.01),
	treeComplexity = c(5, 3, 1),
	bagFraction = 0.6,
	minTrees = 1000,
	maxTrees = 8000,
	tries = 5,
	tryBy = c('learningRate', 'treeComplexity', 'maxTrees', 'stepSize'),
	w = TRUE,
	anyway = FALSE,
	family = 'bernoulli',
	out = 'model',
	cores = 1,
	parallelType = 'doParallel',
	verbose = FALSE,
	...
) {
	
	### setup 
	#########

		# # add dummy variable if using univariate model to avoid errors
		# if (ncol(data)==2) {
			# data$DUMMY <- 1
			# preds <- c(preds, 'DUMMY')
		# }

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		# model weights
		w <- .calcWeights(w, data = data, resp = resp, family = family)

	### generate table of parameterizations
	#######################################
		
		params <- expand.grid(
			learningRate = learningRate,
			treeComplexity = treeComplexity,
			bagFraction = bagFraction,
			maxTrees = maxTrees,
			stringsAsFactors = FALSE
		)
		
	### MAIN
	########

		cores <- min(cores, nrow(params), parallel::detectCores(logical = FALSE))

		if (cores > 1L) {

			`%makeWork%` <- foreach::`%dopar%`
			cl <- parallel::makeCluster(cores, setup_strategy = 'sequential')
			doParallel::registerDoParallel(cl)
			# on.exit(parallel::stopCluster(cl))
			
		} else {
			`%makeWork%` <- foreach::`%do%`
		}

		paths <- .libPaths() # need to pass this to avoid "object '.doSnowGlobals' not found" error!!!
		mcOptions <- list(preschedule = TRUE, set.seed = TRUE, silent = verbose)

		work <- foreach::foreach(
			i = 1L:nrow(params),
			.options.multicore = mcOptions,
			.combine='c',
			.inorder = FALSE,
			.export = c('.trainBrtWorker')
		) %makeWork% {
			.trainBrtWorker(
				i = i,
				params = params,
				data = data,
				preds = preds,
				resp = resp,
				family = family,
				learningRate = learningRate,
				treeComplexity = treeComplexity,
				bagFraction = bagFraction,
				minTrees = minTrees,
				maxTrees = maxTrees,
				tries = tries,
				tryBy = tryBy,
				w = w,
				paths = paths,
				...
			)
				
		}
						
		if (cores > 1L) parallel::stopCluster(cl)

	### collate models
	##################
		
		models <- list()
		tuning <- data.frame()

		for (i in seq_along(work)) {
		
			models[[i]] <- work[[i]]$model
			tuning <- rbind(tuning, work[[i]]$workerTuning)
		
		}
		
	### process models
	##################
	
		if (anyway) {
			origModels <- models
			origTuning <- tuning
		}
	
		# remove non-converged models
		keeps <- which(tuning$converged)
		tuning <- tuning[keeps, , drop=FALSE]
		models <- models[keeps]

		if (length(models) > 0L) {
		
			# remove models with fewer trees than required
			keeps <- which(omnibus::naCompare('>=', tuning$nTrees, minTrees))
			tuning <- tuning[keeps, , drop=FALSE]
			models <- models[keeps]
			
			if (length(models) > 0) {
			
				# sort from best to worst model
				modelOrder <- order(tuning$dev, decreasing=FALSE)
				tuning <- tuning[modelOrder, , drop=FALSE]
				models <- models[modelOrder]
				
				rownames(tuning) <- NULL
				
			}
				
		}
		
		if (anyway & length(models) == 0) {
			models <- origModels
			tuning <- origTuning
			warning('No models converged and/or had sufficient trees.')
		}
		
	### return
	##########
		
		if (verbose) {
			omnibus::say('')
			print(tuning, digits=4)
			omnibus::say('')
		}

		if (length(out) > 1) {
			output <- list()
			if ('models' %in% out) output$models <- models
			if ('model' %in% out) output$model <- if (length(models) > 0) { models[[1]] } else { NA }
			if ('tuning' %in% out) output$tuning <- tuning
			output
		} else if (out == 'models') {
			if (length(models) > 0) {
				models
			} else {
				NULL
			}
		} else if (out == 'model') {
			if (length(models) > 0) {
				models[[1]]
			} else {
				NULL
			}
		} else if (out == 'tuning') {
			tuning
		}

}


#######################
### worker function ###
#######################

.trainBrtWorker <- function(
	i,								# iterator
	params,							# parameterizations
	data,							# data frame
	preds,							# character
	resp,							# character
	family,							# character
	learningRate,					# learning rate
	treeComplexity,					# tree depth
	bagFraction,					# bag fraction
	minTrees,						# minimum number of trees in a model
	maxTrees,						# maximum number of trees in a model
	tries,							# number of times to try if non-convergence
	tryBy,							# one or more of c('learningRate', 'treeComplexity', 'maxTrees', 'stepSize')
	w,								# weights (numeric vector),
	paths,							# .libPaths() output
	...								# other (to pass to step.gbm)
) {

	 # need to call this to avoid "object '.doSnowGlobals' not found" error!!!
	.libPaths(paths)

	# flag to indicate if model converged or not
	converged <- FALSE

	# starter values
	tempLr <- params$learningRate[i]
	tempTc <- params$treeComplexity[i]
	tempBf <- params$bagFraction[i]
	tempMaxTrees <- params$maxTrees[i]
	
	tempStepSize <- 50 # default for n.trees in gbm.step

	# tuning table
	workerTuning <- data.frame()
	
	# by TRY
	numTries <- 0L
	while (numTries <= tries & !converged) {

		numTries <- numTries + 1L

		# try with different parameter combinations
		if (numTries > 1L && !is.null(tryBy)) {

			if ('learningRate' %in% tryBy) tempLr <- tempLr / 10
			if ('treeComplexity' %in% tryBy) tempTc <- max(1, tempTc + ifelse(stats::runif(1) > 0.5, 1, -1))
			if ('maxTrees' %in% tryBy) tempMaxTrees <- round(1.2 * tempMaxTrees)
			if ('stepSize' %in% tryBy) tempStepSize <- round(0.8 * tempStepSize)

		}

		# train model... using tryCatch because model may not converge
		model <- tryCatch(
			model <- dismo::gbm.step(
				data = data,
				gbm.x = preds,
				gbm.y = resp,
				family = family,
				tree.complexity = tempTc,
				learning.rate = tempLr,
				bag.fraction = tempBf,
				max.trees = tempMaxTrees,
				n.trees = tempStepSize,
				plot.main = FALSE,
				plot.folds = FALSE,
				silent = TRUE,
				verbose = TRUE,
				site.weights = w,
				...
			),
			error=function(err) return(NULL)
		)

		# if model training succeeded (model will be gbm object if training succeeded)
		if (!is.null(model)) {

			converged <- TRUE

			# tuning table
			workerTuning <- rbind(
				workerTuning,
				data.frame(
					learningRate = tempLr,
					treeComplexity = tempTc,
					bagFraction = tempBf,
					maxTrees = tempMaxTrees,
					stepSize = tempStepSize,
					nTrees = model$gbm.call$best.trees,
					converged = TRUE,
					deviance = model$cv.statistics$deviance.mean
				)
			)

		} else {
		
			model <- NA
		
			# tuning table
			workerTuning <- rbind(
				workerTuning,
				data.frame(
					learningRate = tempLr,
					treeComplexity = tempTc,
					bagFraction = tempBf,
					maxTrees = tempMaxTrees,
					stepSize = tempStepSize,
					nTrees = NA,
					converged = FALSE,
					deviance = NA
				)
			)
			
		}

	} # while trying to train model

	workerOut <- list(
		list(
			model=model,
			workerTuning=workerTuning
		)
	)
	
	workerOut
	
}
