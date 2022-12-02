#' Calibrate a boosted regression tree (generalized boosting machine) model
#'
#' This function is a wrapper for \code{\link[dismo]{gbm.step}}. It returns the model with best combination of learning rate, tree depth, and bag fraction based on cross-validated deviance. It can also return a table with deviance of different combinations of tuning parameters that were tested, and all of the models tested. See Elith, J., J.R. Leathwick, and T. Hastie. 2008. A working guide to boosted regression trees. \emph{Journal of Animal Ecology} 77:802-813.
#' @param data data frame with first column being response
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
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
#' @param w Either logical in which case \code{TRUE} (default) causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) \emph{or} a numeric list of weights, one per row in \code{data} \emph{or} the name of the column in \code{data} that contains site weights. If \code{FALSE}, then each datum gets a weight of 1.
#' @param anyway Logical. If \code{FALSE} (default), it is possible for no models to be returned if none converge and/or none had a number of trees is >= \code{minTrees}). If \code{TRUE} then all models are returned but with a warning.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest deviance.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest deviance.
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by deviance.
#' }
#' @param cores Integer >= 1. Number of cores to use when calculating multiple models. Default is 1.
#' @param parallelType Either \code{'doParallel'} or \code{'doSNOW'} (default).
#' @param verbose Logical. If \code{TRUE} display progress.
#' @param ... Arguments to pass to \code{\link[dismo]{gbm.step}}.
#'
#' @return If \code{out = 'model'} this function returns an object of class \code{gbm}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and cross-validation deviance for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gbm} object and the data frame. Note that if a model does not converge or does not meet sufficiency criteria (i.e., the number of optimal trees is < \code{minTrees}, then the model is not returned (a \code{NULL} value is returned for \code{'model'} and models are simply missing from the \code{tuning} and \code{models} output.
#' @seealso \code{\link[dismo]{gbm.step}}
#' @examples
#'
#' # The examples below show a very basic modeling workflow. They have been 
#' # designed to work fast, not produce accurate, defensible models.
#' set.seed(123)
#' 
#' ### setup data
#' 
#' # environmental rasters
#' rastFile <- system.file('extdata/madEnv.tif', package='enmSdmX')
#' madEnv <- rast(rastFile)
#' madEnv <- madEnv / 100 # values were rounded to nearest 100th then * by 100
#' 
#' crs <- sf::st_crs(madEnv)
#' 
#' # lemur occurrence data
#' data(lemurs)
#' occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
#' occs <- sf::st_as_sf(occs, coords=c('longitude', 'latitude'), crs=crs)
#' occEnv <- extract(madEnv, occs, ID=FALSE)
#' occEnv <- occEnv[complete.cases(occEnv), ]
#' 	
#' # create 10000 background sites (or as many as raster can support)
#' bgEnv <- terra::spatSample(madEnv, 20000)
#' bgEnv <- bgEnv[complete.cases(bgEnv), ]
#' bgEnv <- bgEnv[1:min(10000, nrow(bgEnv)), ]
#' 
#' # collate occurrences and background sites
#' presBg <- data.frame(
#' 	presBg = c(
#'    rep(1, nrow(occEnv)),
#'    rep(0, nrow(bgEnv))
#'    )
#' )
#' 
#' env <- rbind(occEnv, bgEnv)
#' env <- cbind(presBg, env)
#' 
#' predictors <- c('bio1', 'bio12')
#' 
#' ## MaxEnt
#' mx <- trainMaxEnt(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	regMult = 1, # too few values for reliable model, but fast
#' 	verbose = TRUE
#' )
#' 
#' ## generalized linear model (GLM)
#' # Normally, we'd center and standardize variables before modeling.
#' gl <- trainGlm(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	verbose = TRUE
#' )
#' 
#' ## generalized additive model (GAM)
#' ga <- trainGam(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	verbose = TRUE
#' )
#' 
#' ## natural splines
#' nat <- trainNs(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	verbose = TRUE
#' )
#' 
#' ## boosted regression trees
#' envSub <- env[1:2000, ] # subsetting data to run faster
#' brt <- trainBrt(
#' 	data = envSub,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	learningRate = 0.001, # too few values for reliable model(?)
#' 	treeComplexity = 2, # too few values for reliable model, but fast
#' 	minTrees = 1200, # minimum trees for reliable model(?), but fast
#' 	maxTrees = 1200, # too small for reliable model(?), but fast
#' 	tryBy = 'treeComplexity',
#' 	anyway = TRUE, # return models that did not converge
#' 	verbose = TRUE
#' )
#' 
#' ## random forests
#' rf <- trainRf(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = predictors,
#' 	verbose = TRUE
#' )
#' 
#' ## make maps of models
#' 
#' mxMap <- predictEnmSdm(mx, madEnv)
#' glMap <- predictEnmSdm(gl, madEnv)
#' gaMap <- predictEnmSdm(ga, madEnv)
#' natMap <- predictEnmSdm(nat, madEnv)
#' brtMap <- predictEnmSdm(brt, madEnv)
#' rfMap <- predictEnmSdm(rf, madEnv)
#' 
#' maps <- c(
#' 	mxMap,
#' 	glMap,
#' 	gaMap,
#' 	natMap,
#' 	brtMap,
#' 	rfMap
#' )
#' 
#' names(maps) <- c('MaxEnt', 'GLM', 'GAM', 'Natural Splines', 'BRTs', 'RFs')
#' fun <- function() plot(occs[1], col='black', add=TRUE)
#' plot(maps, fun=fun)
#' 
#' ## compare model responses to BIO12 (mean annual precipitation)
#' 
#' # make a data frame holding all other variables at mean across occurrences,
#' # varying only BIO12
#' occEnvMeans <- colMeans(occEnv, na.rm=TRUE)
#' occEnvMeans <- rbind(occEnvMeans)
#' occEnvMeans <- as.data.frame(occEnvMeans)
#' climFrame <- occEnvMeans[rep(1, 100), ]
#' rownames(climFrame) <- NULL
#' 
#' minBio12 <- min(env$bio12)
#' maxBio12 <- max(env$bio12)
#' climFrame$bio12 <- seq(minBio12, maxBio12, length.out=100)
#' 
#' predMx <- predictEnmSdm(mx, climFrame)
#' predGl <- predictEnmSdm(gl, climFrame)
#' predGa <- predictEnmSdm(ga, climFrame)
#' predNat <- predictEnmSdm(nat, climFrame)
#' predBrt <- predictEnmSdm(brt, climFrame)
#' predRf <- predictEnmSdm(rf, climFrame)
#' 
#' 
#' plot(climFrame$bio12, predMx,
#' xlab='BIO12', ylab='Prediction', type='l', ylim=c(0, 1))
#' 
#' lines(climFrame$bio12, predGl, lty='dotted', col='blue')
#' lines(climFrame$bio12, predGa, lty='dashed', col='green')
#' lines(climFrame$bio12, predNat, lty=4, col='purple')
#' lines(climFrame$bio12, predBrt, lty=5, col='orange')
#' lines(climFrame$bio12, predRf, lty=6, col='cyan')
#' 
#' legend(
#'    'topleft',
#'    inset = 0.01,
#'    legend = c(
#' 	'MaxEnt',
#' 	'GLM',
#' 	'GAM',
#' 	'NS',
#' 	'BRT',
#' 	'RF'
#'    ),
#'    lty = 1:6,
#'    col = c(
#' 	'black',
#' 	'blue',
#' 	'green',
#' 	'purple',
#' 	'orange',
#' 	'cyan'
#'    ),
#'    bg = 'white'
#' )
#' 
#' 
#' @export
trainBrt_parallel <- function(
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
	parallelType = 'doSnow',
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
		if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
		if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

		# model weights
		if (class(w)[1] == 'logical') {
			w <- if (w) {
				c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
			} else {
				rep(1, nrow(data))
			}
		} else if (class(w) == 'character') {
			w <- data[ , w]
		}

		w <- w / max(w)

	### generate table of parameterizations
	#######################################
		
		params <- expand.grid(learningRate=learningRate, treeComplexity=treeComplexity, bagFraction=bagFraction, maxTrees=maxTrees)
		
	### MAIN
	########

		cores <- min(cores, nrow(params), parallel::detectCores(logical = FALSE))

		if (cores > 1L) {

			`%makeWork%` <- foreach::`%dopar%`
			cl <- parallel::makeCluster(cores, setup_strategy = 'sequential')

			if (parallelType == 'doParallel') {
				doParallel::registerDoParallel(cl)
				# opts <- NULL
				say('doParallel')
			} else if (parallelType == 'doSNOW') {
				doSNOW::registerDoSNOW(cl)
				# opts <- list(progress = verbose)
				say('doSNOW')
			}
			
		} else {
			`%makeWork%` <- foreach::`%do%`
			# opts <- NULL
			say('1 core')
		}

		paths <- .libPaths() # need to pass this to avoid "object '.doSnowGlobals' not found" error!!!
		mcOptions <- list(preschedule = TRUE, set.seed = TRUE, silent = verbose)

		work <- foreach::foreach(
			i = 1L:nrow(params),
			# .options.snow = opts,
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
						
		if (cores > 1) parallel::stopCluster(cl)

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
	numTries <- 0
	while (numTries <= tries & !converged) {

		numTries <- numTries + 1

		# try with different parameter combinations
		if (numTries > 1 && !is.null(tryBy)) {

			if ('learningRate' %in% tryBy) tempLr <- tempLr / 10
			if ('treeComplexity' %in% tryBy) tempTc <- max(1, tempTc + ifelse(stats::runif(1) > 0.5, 1, -1))
			if ('maxTrees' %in% tryBy) tempMaxTrees <- round(1.2 * tempMaxTrees)
			if ('stepSize' %in% tryBy) tempStepSize <- round(0.8 * tempStepSize)

		}

		# train model... using tryCatch because model may not converge
		model <- tryCatch(
			model <- dismo::gbm.step(
				data=data,
				gbm.x=preds,
				gbm.y=resp,
				family=family,
				tree.complexity=tempTc,
				learning.rate=tempLr,
				bag.fraction=tempBf,
				max.trees=tempMaxTrees,
				n.trees=tempStepSize,
				plot.main=FALSE,
				plot.folds=FALSE,
				silent=TRUE,
				verbose=TRUE,
				site.weights=w,
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
