#' Calibrate a random forest model
#'
#' @description This function trains a random forest model. It identifies the optimal number of trees and value for \code{mtry} (number of variables sampled as candidates at each split) using out-of-bag error (OOB). The number of trees in each candidate model is set by the user with argument \code{numTrees}. The number of predictors to test per split, \code{mtry}, is found by exploring a range of values. If the response (\code{y}) is a factor, the starting value for \code{mtry} is \code{max(1, floor(p / 3))}, where \code{p} is the number of predictors. If the response is not a factor, the starting value is \code{max(1, floor(sqrt(p)))}. Values y\code{mtryIncrement} argument until the total number of predictors is used.  See \code{\link[randomForest]{randomForest}} for more details.
#'
#' The output of the function is any or all of: a table with out-of-bag (OOB) error of evaluated models; all evaluated models; and/or the single model with the lowest OOB error.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response variable. The default is to use the first column in \code{data} as the response.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param family Character. If "\code{binomial}" then the response is converted to a binary factor with levels 0 and 1. Otherwise, this argument has no effect.
#' @param numTrees Vector of number of trees to grow. All possible combinations of \code{mtry} and \code{numTrees} will be assessed.
#' @param mtryIncrement Positive integer (default is 2).  Number of predictors to add to \code{mtry} until all predictors are in each tree.
#' @param w Weights. For random forests, weights are simply used as relative probabilities of selecting a row in \code{data} to be used in a particular tree. This argument takes any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest out-of-bag (OOB) error rate.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest OOB.
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by OOB error rate.
#' }
#' @param cores Number of cores to use. Default is 1. If you have issues when \code{cores} > 1, please see the \code{\link{troubleshooting_parallel_operations}} guide.
#' @param verbose Logical. If \code{TRUE} then display progress for finding optimal value of \code{mtry}.
#' @param ... Arguments to pass to \code{\link[randomForest]{randomForest}}.
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these.
#'
#' @seealso \code{\link[randomForest]{randomForest}}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export
trainRF <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	numTrees = c(500, 1000),
	mtryIncrement = 2,
	w = TRUE,
	family = 'binomial',
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {

	# response and predictors
	if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
	if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

	# model weights
	w <- .calcWeights(w, data = data, resp = resp, family = family)
	
	# binomial response
	if (family == 'binomial') data[ , resp] <- factor(data[ , resp], levels=0:1, labels=c('0', '1'))

	x <- data[ , preds, drop = FALSE]
	y <- data[ , resp, drop = TRUE]

	mtryStart <- if (!is.null(y) && !is.factor(y)) max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))
	mtryEnd <- length(preds)
	mtries <- seq(mtryStart, mtryEnd, mtryIncrement)
	if (utils::tail(mtries)[1L] != mtryEnd) mtries <- c(mtries, mtryEnd)

	params <- expand.grid(numTrees = numTrees, mtry = mtries, stringsAsFactors = FALSE)

	### parallelization
	###################
			
		cores <- min(cores, nrow(params), parallel::detectCores(logical = FALSE))

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

	### grow forest
	###############
	
		work <- foreach::foreach(
			i = 1L:nrow(params),
			.options.multicore = mcOptions,
			.combine = 'c',
			.inorder = FALSE,
			.export = c('.trainRfWorker')
		) %makeWork% {
			.trainRfWorker(
				i = i,
				x = x,
				y = y,
				w = w,
				params = params,
				paths = paths,
				modelOut = ('models' %in% out | 'model' %in% out),
				...
			)
		}
	
		# if (cores > 1L) parallel::stopCluster(cl)
	
		# tuning table
		tuning <- data.frame(
			numTrees = work[[1L]]$numTrees,
			mtry = work[[1L]]$mtry,
			oobError = work[[1L]]$oob
		)
		
		if (length(work) > 1L) {
			for (i in 2L:length(work)) {
				
				tuning <- rbind(
					tuning,
					data.frame(
						numTrees = work[[i]]$numTrees,
						mtry = work[[i]]$mtry,
						oobError = work[[i]]$oob
					)
				)
				
			}
		}
	
		bestOrder <- order(tuning$oobError, tuning$numTrees)
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

.trainRfWorker <- function(
	i,
	x,
	y,
	w,
	params,
	paths,
	modelOut,
	...
) {

	.libPaths(paths)

	ntree <- params$numTrees[i]
	mtry <- params$mtry[i]

	model <- randomForest::randomForest(
		x = x,
		y = y,
		mtry = mtry,
		ntree = ntree,
		weights = w,
		strata = y,
		...
	)

	oob <- model$err.rate[ntree, 1L]

	# out
	out <- if (modelOut) {
		
		list(
			list(
				model = model,
				numTrees = ntree,
				mtry = mtry,
				oob = oob
			)
		)
		
	} else {
	
		data.frame(
			numTrees = ntree,
			mtry = mtry,
			oob = oob
		)
	
	}
		
	out
	
}
