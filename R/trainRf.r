#' Calibrate a random forest model
#'
#' This function trains a random forest model. It identifies the optimal value for \code{mtry} (number of variables sampled as candidates at each split) using out-of-bag error (OOB). See \code{\link[randomForest]{randomForest}} for more details.
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response varoable. The default is to use the first column in \code{data} as the response.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param family Character. If "\code{binomial}" then the response is converted to a binary factor with levels 0 and 1. Otherwise, this argument has no effect.
#' @param numTrees Vector of number of trees to grow. All possible combinations of \code{mtry} and \code{numTrees} will be assessed.
#' @param mtryFactor Positive integer (default is 2).  Number of predictors to add to \code{mtry} until all predictors are in each tree.
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
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by OOB error rate.
#' }
#' @param cores Number of cores to use. Default is 1.
#' @param parallelType Either \code{'doParallel'} (default) or \code{'doSNOW'}. Issues with parallelization might be solved by trying the non-default option.
#' @param verbose Logical. If \code{TRUE} then display progress for finding optimal value of \code{mtry}.
#' @param ... Arguments to pass to \code{\link[randomForest]{randomForest}}.
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these.
#'
#' @seealso \code{\link[randomForest]{randomForest}}
#'
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
trainRf <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	numTrees = c(500, 1000),
	mtryFactor = 2,
	w = TRUE,
	cores = 1,
	parallelType = 'doParallel',
	out = 'model',
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
	mtries <- seq(mtryStart, mtryEnd, mtryFactor)
	if (tail(mtries)[1L] != mtryEnd) mtries <- c(mtries, mtryEnd)

	params <- expand.grid(numTrees = numTrees, mtry = mtries, stringsAsFactors = FALSE)

	### parallelization
	###################
			
		cores <- min(cores, nrow(params), parallel::detectCores(logical = FALSE))

		if (cores > 1L) {

			`%makeWork%` <- foreach::`%dopar%`
			cl <- parallel::makeCluster(cores, setup_strategy = 'sequential')

			if (tolower(parallelType) == 'doparallel') {
				doParallel::registerDoParallel(cl)
			} else if (tolower(parallelType) == 'dosnow') {
				doSNOW::registerDoSNOW(cl)
			} else {
				stop('Argument "parallelType" must be either "doParallel" or "doSNOW".')
			}
			
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
	
		if (cores > 1L) parallel::stopCluster(cl)
	
		# tuning table
		tuning <- data.frame(
			numTrees = work[[1L]]$numTrees,
			oobError = work[[1L]]$oob
		)
		
		if (length(work) > 1L) {
			for (i in 2L:length(work)) {
				
				tuning <- rbind(
					tuning,
					data.frame(
						numTrees = work[[i]]$numTrees,
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
