#' Calibrate a MaxEnt (ver 3.3.3+ or "maxent") model using AICc
#'
#' This function calculates the "best" Maxent model using AICc across all possible combinations of a set of master regularization parameters and feature classes. The best model has the lowest AICc, with ties broken by number of features (fewer is better), regularization multiplier (higher better), then finally the number of coefficients (fewer better). The function can return the best model (default), a list of models created using all possible combinations of feature classes and regularization multipliers, and/or a data frame with tuning statistics for each model. Models in the list and in the data frame are sorted from best to worst. The function requires the \code{maxent} jar file (see \emph{Details}).  Its output is any or all of: a table with AICc for all evaluated models; all models evaluated in the "selection" phase; and/or the single model with the lowest AICc.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response variable. The default is to use the first column in \code{data} as the response.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh}) or any combination of \code{lpqht}, where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If \code{TRUE} (default) then test all possible combinations of classes (note that all tested models will at least have linear features). If \code{FALSE} then use the classes provided (these will not vary between models).
#' @param dropOverparam Logical, if \code{TRUE} (default), drop models if they have more coefficients than training occurrences. It is possible for no models to fulfill this criterion, in which case no models will be returned.
#' @param anyway Logical. Same as \code{dropOverparam} (included for backwards compatibility. If \code{NULL} (default), then the value of \code{dropOverparam} will take precedence. If \code{TRUE} or \code{FALSE} then \code{anyway} will override the value of \code{dropOverparam}.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by AICc.
#' }
#' @param forceLinear Logical. If \code{TRUE} (default) then require any tested models to include at least linear features.
#' @param jackknife Logical. If \code{TRUE} (default) the the returned model will be also include jackknife testing of variable importance.
#' @param arguments \code{NULL} (default) or a character list. Options to pass to \code{maxent()}'s \code{args} argument. (Do not include \code{l}, \code{p}, \code{q}, \code{h}, \code{t}, \code{betamultiplier}, or \code{jackknife}!)
#' @param scratchDir Character. Directory to which to write temporary files. Leave as NULL to create a temporary folder in the current working directory.
#' @param cores Number of cores to use. Default is 1.
#' @param verbose Logical. If \code{TRUE} report progress and AICc table.
#' @param ... Extra arguments. Not used.
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these.
#'
#' @details This function is a wrapper for \code{maxent()}. That function relies on a maxent \code{jar} file in the folder \code{./library/dismo/java}. See \code{\link[dismo]{maxent}} for more details. The \code{maxent} function creates a series of files on disk for each model. This function assumes you do not want those files, so deletes most of them. However, there is one that cannot be deleted and the normal ways of changing its permissions in \pkg{R} do not work. So the function simply writes over that file (which is allowed) to make it smaller. Regardless, if you run many models your temporary directory (argument \code{scratchDir}) can fill up and require manual deletion.
#'
#' @seealso \code{\link[dismo]{maxent}}
#' @references
#' Warren, D.L. and S.N. Siefert. 2011. Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. \emph{Ecological Applications} 21:335-342. \doi{10.1890/10-1171.1}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export
trainMaxEnt <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10),
	classes = 'default',
	testClasses = TRUE,
	dropOverparam = TRUE,
	anyway = TRUE,
	forceLinear = TRUE,
	jackknife = FALSE,
	arguments = NULL,
	scratchDir = NULL,
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########
	
		# need this bc on some platforms fails CRAN checks bc rJava isn't imported and if we import it we need to use it
		# so dumb!
		if (FALSE) {
			rJava::clone(data)
		}

		if (!is.null(anyway)) dropOverparam <- anyway
		
		# response and predictors
		if (inherits(data, 'data.table')) data <- as.data.frame(data)
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		# get response and predictors
		presBg <- data[ , resp]
		data <- data[ , preds, drop=FALSE]

		## collate all presences
		allPres <- data[presBg == 1, , drop=FALSE]

		## collate all background sites
		allBg <- data[presBg == 0, , drop=FALSE]

	### generate table of parameterizations
	#######################################
		
		### get combinations of features to test for each regularization multiplier
		classesToTest <- if (classes == 'default') {
			c('l', 'p', 'q', 'h')
		} else {
			unlist(strsplit(classes, ''))
		}

		if (any('p' %in% classesToTest) & ncol(data) == 1) {
			product <- FALSE
			warning('Data has only one variable so forcing product features to FALSE.')
		}

		# create df of 1/0 to indicate each combination of classes to test
		if (testClasses) {
			classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)), stringsAsFactors = FALSE)
			classGrid <- classGrid[-which(rowSums(classGrid) == 0), , drop=FALSE]
		} else {
			classGrid <- data.frame(matrix(rep(1, length(classesToTest)), nrow=1))
		}

		names(classGrid) <- classesToTest

		if (forceLinear & any(classGrid$l == 0)) classGrid <- classGrid[-which(classGrid$l == 0), , drop=FALSE]

		tuning <- data.frame()
		
		# by beta
		for (thisRegMult in regMult) {

			# by feature combination
			for (countParam in 1:nrow(classGrid)) {

				# classes
				classes <- paste0(
					ifelse('l' %in% names(classGrid) && classGrid$l[countParam] == 1, 'l', ''),
					ifelse('p' %in% names(classGrid) && classGrid$p[countParam] == 1, 'p', ''),
					ifelse('q' %in% names(classGrid) && classGrid$q[countParam] == 1, 'q', ''),
					ifelse('h' %in% names(classGrid) && classGrid$h[countParam] == 1, 'h', ''),
					ifelse('t' %in% names(classGrid) && classGrid$t[countParam] == 1, 't', '')
				)
				
				tuning <- rbind(
					tuning,
						data.frame(
						regMult=thisRegMult,
						classes=classes
					)
				)
				
			}
			
		}

	### parallelization
	###################
			
		cores <- min(cores, nrow(tuning), parallel::detectCores(logical = FALSE))

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

	### train models
	################
	
		work <- foreach::foreach(
			i = 1L:nrow(tuning),
			.options.multicore = mcOptions,
			.combine = 'c',
			.inorder = FALSE,
			.export = c('.trainMaxEntWorker'),
			.packages = c('rJava')
		) %makeWork% {
			.trainMaxEntWorker(
				i = i,
				scratchDir = scratchDir,
				tuning = tuning,
				presBg = presBg,
				data = data,
				allPres = allPres,
				allBg = allBg,
				jackknife = jackknife,
				arguments = arguments
			)
		}

		if (cores > 1L) parallel::stopCluster(cl)
	
	### collate results
	###################

		models <- list()
		tuning <- data.frame()

		for (i in seq_along(work)) {
		
			models[[i]] <- work[[i]]$model
			tuning <- rbind(tuning, work[[i]]$thisTuning)
		
		}
		
	### process models
	##################
		
		# sort from best to worst model
		modelOrder <- order(tuning$AICc, tuning$numClasses, tuning$regMult, tuning$numCoeff, decreasing=c(FALSE, FALSE, TRUE, FALSE))
		tuning <- tuning[modelOrder, ]
		models <- models[modelOrder]

		# remove models with more parameters than data points that have more than 0 parameters
		if (dropOverparam) {
		
			topModel <- models[[1]]
			topTuning <- tuning[1, , drop=FALSE]

			overparamModels <- which(tuning$n < tuning$numCoeff)
			if (length(overparamModels) > 0) {
				tuning <- tuning[-overparamModels, ]
				models <- models[-overparamModels]
			}

			if (length(models) == 0 & anyway) {
				warning('No models had fewer coefficients than predictors. Returning best model anyway.', immediate.=TRUE)
				model <- topModel
				tuning <- topTuning
			} else if (length(models) == 0 & !anyway) {
				warning('No models had fewer coefficients than predictors. NA model(s) returned.', immediate.=TRUE)
				model <- NULL
			} else {
				model <- models[[1]]
			}
				
		}
		
		# AICc weights
		if (nrow(tuning) > 0) {

			tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
			tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
			tuning$aicWeight <- tuning$relLike / sum(tuning$relLike)

			rownames(tuning) <- NULL
		
		}

		if (verbose) {

			omnibus::say('')
			print(tuning, digits=3)
			omnibus::say('')

		}

		# remove temporary files... note that "species.lambda" file cannot be removed unless R is closed, so we'll just make it smaller to reduce disk space usage
		unlink(scratchDir, recursive=TRUE, force=TRUE)

	### return
	##########
		
		if (length(out) > 1) {
			output <- list()
			if ('models' %in% out) output$models <- models
			if ('model' %in% out) output$model <- model
			if ('tuning' %in% out) output$tuning <- tuning
			output
		} else if (out == 'models') {
			models
		} else if (out == 'model') {
			model
		} else if (out == 'tuning') {
			tuning
		}

}

#######################
### worker function ###
#######################

.trainMaxEntWorker <- function(
	i,								# iterator
	scratchDir,						# master temp path
	tuning,							# tuning data frame
	presBg,							# vector with 1 (present) / 0 (background)
	data,							# df with all presence/background environmental data
	allPres,						# df with all presence environmental data
	allBg,							# df with all background environmental data
	jackknife,						# logical
	arguments						# string of arguments for maxent
) {
	
	thisRegMult <- tuning$regMult[i]
	thisClasses <- tuning$classes[i]
	
	# create scratch directory
	thisScratchDir <- if (is.null(scratchDir)) {
		base::tempfile(pattern='/_maxentTempFiles/')
	} else {
		base::tempfile(pattern='/_maxentTempFiles/', tmpdir=scratchDir)
	}
	
	omnibus::dirCreate(thisScratchDir, '/plots')

	# get parameters
	params <- c(
		paste0('betamultiplier=', thisRegMult),
		paste0('linear=', ifelse(grepl('l', thisClasses), 'true', 'false')),
		paste0('product=', ifelse(grepl('p', thisClasses), 'true', 'false')),
		paste0('quadratic=', ifelse(grepl('q', thisClasses), 'true', 'false')),
		paste0('hinge=', ifelse(grepl('h', thisClasses), 'true', 'false')),
		paste0('threshold=', ifelse(grepl('t', thisClasses), 'true', 'false')),
		paste0('jackknife=', ifelse(jackknife, 'true', 'false'))
	)

	params <- if (!is.null(arguments)) {
		c(params, arguments)
	} else {
		params
	}

	# train model
	model <- dismo::maxent(
		x=data,
		p=as.vector(presBg),
		path=thisScratchDir,
		args=params,
		silent=TRUE
	)

	## predict to training (and maybe test presences)
	predPres <- dismo::predict(
		object=model,
		x=allPres,
		na.rm=TRUE,
		arguments='outputformat=raw'
	)

	## predict to background
	predBg <- dismo::predict(
		object=model,
		x=allBg,
		na.rm=TRUE,
		arguments='outputformat=raw'
	)

	rawSum <- sum(c(predPres, predBg), na.rm=TRUE)

	## log likelihood
	ll <- sum(log(predPres / rawSum), na.rm=TRUE)

	## number of parameters
	numCoeff <- 0

	for (thisLambda in model@lambdas) { # for each line in lambda object

		# if not a meta-data line
		if (!grepl(thisLambda, pattern='linearPredictorNormalizer') & !grepl(thisLambda, pattern='densityNormalizer') & !grepl(thisLambda, pattern='numBackgroundPoints') & !grepl(thisLambda, pattern='entropy')) {

			split <- strsplit(thisLambda, ', ')
			paramValue <- as.numeric(split[[1]][2])
			if (paramValue !=0) numCoeff <- numCoeff + 1 # increment number of parameters

		}

	}

	# AICc
	AICc <- -2 * ll + 2 * numCoeff + (2 * numCoeff * (numCoeff + 1)) / (sum(presBg) - numCoeff - 1)

	# remove temporary files... note that "species.lambda" file cannot be removed unless R is closed, so we'll just make it smaller to reduce disk space usage
	utils::write.csv(NULL, paste0(thisScratchDir, '/species.lambdas'))
	if (file.exists(paste0(thisScratchDir, '/presences'))) utils::write.csv(NULL, paste0(scratchDir, '/', i, '/presences'))
	if (file.exists(paste0(thisScratchDir, '/absences'))) utils::write.csv(NULL, paste0(scratchDir, '/', i, '/absences'))
	unlink(paste0(thisScratchDir, '/plots'), recursive=TRUE, force=TRUE)
	unlink(thisScratchDir, recursive=TRUE, force=TRUE)

	out <- list(list(
		model=model,
		thisTuning=data.frame(
			regMult=thisRegMult,
			classes=thisClasses,
			numClasses=nchar(thisClasses),
			numCoeff=numCoeff,
			logLik=ll,
			AICc=AICc
		)
	))
	
	out

}
