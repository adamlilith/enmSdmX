#' Generic predict function for SDMs/ENMs
#'
#' This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, MaxEnt, and more.
#'
#' @param model		Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{maxnet}, \code{prcomp}, \code{kde}, \code{gbm}, and possibly others (worth a try!).
#'
#' @param newdata	Data frame or matrix, or \code{SpatRaster} with data to which to predict.
#'
#' @param maxentFun	This argument is only used if the \code{model} object is a MaxEnt model; otherwise, it is ignored. It takes a value of either \code{'terra'}, in which case a MaxEnt model is predicted using the default \code{predict} function from the \pkg{terra} package, or \code{'enmSdmX'} in which case the function \code{\link[enmSdmX]{predictMaxEnt}} function from the \pkg{enmSdmX} package (this package) is used.
#'
#' @param cores	 Integer >= 1. Number of cores to use when calculating multiple models. Default is 1. This is forced to 1 if \code{newdata} is a \code{SpatRaster} (i.e., as of now, there is no parallelization when predicting to a raster... sorry!).  If you have issues when \code{cores} > 1, please see the \code{\link{troubleshooting_parallel_operations}} guide.
#'
#' @param nrows		Number of rows of \code{newdata} to predict at a time. This is only used if \code{newdata} is a \code{data.frame} or \code{matrix}. The default is to predict all rows at once, but for very large data frames/matrices this can lead to memory issues in some cases. By setting the number of rows, \code{newdata} is divided into chunks, and predictions made to each chunk, which may ease memory limitations. This can be combined with multi-coring (which will increase memory requirements). In this case, all cores combined will get \code{nrows} of data. How many rows are too many? You will have to decide depending on the capabilities of your system. For example, predicting the outcome of a GLM on data with 10E6 rows may be fine, but predicting a PCA (with multiple axes) to the data data may require too much memory.
#'
#' @param paths Locations where packages are stored. This is typically not useful to the general user, and is only supplied for when the function is called as a functional.
#'
#' @param ... Arguments to pass to the algorithm-specific \code{predict} function.
#'
#' @return Numeric or \code{SpatRaster}.
#'
#' @seealso \code{\link[stats]{predict}} from the \pkg{stats} package, \code{\link[terra]{predict}} from the \pkg{terra} package, \code{\link[enmSdmX]{predictMaxEnt}}, \code{\link[enmSdmX]{predictMaxNet}}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export
predictEnmSdm <- function(
	model,
	newdata,
	maxentFun = 'terra',
	cores = 1,
	nrows = nrow(newdata),
	paths = .libPaths(),
	...
) {
	
	cores <- min(cores, parallel::detectCores(logical = FALSE))
	i <- NULL
	
	# doing this for cases where this function is called as a functional in a multi-core state and packages are in a custom location
	ells <- list(...)
	.libPaths(paths)

	### predict to each subset of rows
	##################################
	n <- nrow(newdata)
	if (!inherits(newdata, c('SpatRaster')) & nrows < n) {

		# predict to each set of rows
		rowJobs <- ceiling(n / nrows)
		for (countRowJob in 1L:rowJobs) {

			thisNewData <- newdata[(1 + (countRowJob - 1) * nrows):min(n, countRowJob * nrows), , drop=FALSE]
			thisOut <- predictEnmSdm(model = model, newdata = thisNewData, maxentFun = maxentFun, cores = cores, nrows = nrow(thisNewData), ...)
			
			if (exists('out', inherits = FALSE)) {
				# if (!inherits(newdata, c('SpatRaster'))) {
					# out <- rbind(out, thisOut)
				# } else {
					out <- c(out, thisOut)
				# }
			} else {
				out <- thisOut
			}
		
		}
	
	### parallelization
	###################
	} else if (cores > 1L & nrows > n) {

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

		paths <- .libPaths() # need to pass this to avoid "object '.doSnowGlobals' not found" error!!!
		mcOptions <- list(preschedule = TRUE, set.seed = TRUE, silent = TRUE)

		nPerJob <- floor(nrow(newdata) / cores)
		jobs <- rep(1:cores, each=nPerJob)
		if (length(jobs) < nrow(newdata)) jobs <- c(jobs, rep(cores, nrow(newdata) - length(jobs)))

		combine <- if (inherits(model, 'prcomp')) {
			'rbind'
		} else {
			'c'
		}
		
		out <- foreach::foreach(
			i = 1L:cores,
			.options.multicore = mcOptions,
			.combine = combine,
			.multicombine = TRUE,
			.inorder = TRUE,
			.export=c('predictEnmSdm', 'predictMaxEnt', 'predictMaxNet'),
			.packages=c('enmSdmX')
		) %makeWork% {
			enmSdmX::predictEnmSdm(
				i = i,
				model = model,
				newdata = newdata[jobs == i, , drop = FALSE],
				maxentFun = maxentFun,
				cores = 1L,
				paths = paths,
				...
			)
		}

		# if (cores > 1L) parallel::stopCluster(cl)

	### single-core/raster
	######################
	} else {
		
		if (inherits(newdata, 'data.table')) newdata <- as.data.frame(newdata)

		# GAM
		if (inherits(model, 'gam')) {

			out <- if (inherits(newdata, c('SpatRaster'))) {
				terra::predict(newdata, model, type='response', ...)
			} else {
				mgcv::predict.gam(model, newdata, type='response', ...)
			}

		# GLM
		} else if (inherits(model, c('glm'))) {

			out <- if (inherits(newdata, c('SpatRaster'))) {
				terra::predict(newdata, model, type='response', ...)
			} else {
				stats::predict.glm(model, newdata, type='response', ...)
			}

		# LM
		} else if (inherits(model, 'lm')) {

			out <- stats::predict.lm(model, newdata, ...)

		# BRT
		} else if (inherits(model, 'gbm')) {

			if (inherits(newdata, 'SpatRaster')) {
				nd <- terra::as.data.frame(newdata, na.rm=FALSE)
				notNas <- which(stats::complete.cases(nd))
				nd <- nd[notNas, , drop=FALSE]
				preds <- gbm::predict.gbm(model, nd, n.trees=model$gbm.call$n.trees, type='response', ...)
				out <- newdata[[1L]]
				out[] <- NA
				out <- setValueByCell(out, preds, cell=notNas, format='raster')
			} else {
				out <- gbm::predict.gbm(model, newdata, n.trees=model$gbm.call$n.trees, type='response', ...)
			}

			
		# KDE from ks package
		} else if (inherits(model, 'kde')) {

			# hack... not calling ks functions explicitly at least once in package generates warning
			if (FALSE) fhat <- ks::kde(stats::rnorm(100))
			predictKde <- utils::getFromNamespace('predict.kde', 'ks')
			out <- predictKde(model, x=as.matrix(newdata), ...)

		# Maxent
		} else if (inherits(model, c('MaxEnt', 'MaxEnt_model'))) {

			out <- if (maxentFun == 'terra') {
				terra::predict(model, newdata, ...)
			} else if (maxentFun == 'enmSdmX') {
				predictMaxEnt(model, newdata, ...)
			}

		# MaxNet
		} else if (inherits(model, 'maxnet')) {

			out <- predictMaxNet(model = model, newdata = newdata, ...)

		# random forest in randomForest package
		} else if (inherits(model, 'randomForest')) {

			predictRandomForest <- utils::getFromNamespace('predict.randomForest', 'randomForest')

			if (inherits(newdata, 'SpatRaster')) {

				nd <- terra::as.data.frame(newdata, na.rm=FALSE)
				notNas <- which(stats::complete.cases(nd))
				nd <- nd[notNas, , drop=FALSE]
				preds <- predictRandomForest(object = model, newdata = nd, type = 'prob', ...)
				preds <- preds[ , '1']
				out <- newdata[[1L]]
				out[] <- NA
				out <- setValueByCell(out, preds, cell=notNas, format='raster')
				names(out) <- 'randomForest'
				
			} else {
				preds <- predictRandomForest(object = model, newdata = newdata, type = 'prob', ...)
				out <- preds[ , '1']
			}

		# random forest in ranger package
		} else if (inherits(model, 'ranger')) {

			predictRanger <- utils::getFromNamespace('predict.ranger', 'ranger')
			binary <- if (!exists('binary', where=model, inherits=FALSE)) {
				FALSE
			} else {
				model$binary
			}


			if (inherits(newdata, 'SpatRaster')) {

				nd <- as.data.frame(newdata, na.rm=FALSE)
				notNas <- which(stats::complete.cases(nd))
				nd <- nd[notNas, , drop=FALSE]
				
				preds <- predictRanger(model, data=nd, predict.all=binary, type='response', ...)
				
				if (binary) {
					preds <- rowMeans(preds$predictions)
					preds <- preds - 1
				}
				
				out <- newdata[[1L]]
				out[] <- NA
				out <- setValueByCell(out, preds, cell=notNas, format='raster')
				names(out) <- 'ranger'
				
			} else {
			
				out <- predictRanger(model, data=newdata, predict.all=binary, type='response', ...)
				
				if (binary) {
					out <- rowMeans(out$predictions)
					out <- out - 1
				}
			
			}

		# anything else!
		} else {
			out <- do.call('predict', args=list(object=model, newdata=newdata, type='response', ...))
		}
		
	} # single-core

	out

}
