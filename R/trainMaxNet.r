#' Calibrate a MaxNet (MaxEnt) model using AICc
#'
#' This function calculates the "best" MaxNet model using AICc across all possible combinations of a set of master regularization parameters and feature classes. The "best" model has the lowest AICc, with ties broken by number of features (fewer is better), regularization multiplier (higher better), then finally the number of coefficients (fewer better). The function can return the best model (default), a list of models created using all possible combinations of feature classes and regularization multipliers, and/or a data frame with tuning statistics for each model. Models in the list and in the data frame are sorted from best to worst. Its output is any or all of: a table with AICc for all evaluated models; all models evaluated in the "selection" phase; and/or the single model with the lowest AICc.
#'
#' @param data  Data frame or matrix. Contains a column indicating whether each row is a presence (1) or background (0) site, plus columns for environmental predictors.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh}) or any combination of \code{lpqht}, where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If \code{TRUE} (default) then test all possible combinations of classes (note that all tested models will at least have linear features). If \code{FALSE} then use the classes provided (these will not vary between models).
#' @param dropOverparam Logical, if \code{TRUE} (default), drop models if they have more coefficients than training occurrences. It is possible for no models to fulfill this criterion, in which case no models will be returned.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning parameters, one row per model, sorted by AICc.
#' }
#' @param forceLinear Logical. If \code{TRUE} (default) then require any tested models to include at least linear features.
#' @param cores Number of cores to use. Default is 1.
#' @param parallelType Either \code{'doParallel'} (default) or \code{'doSNOW'}. Issues with parallelization might be solved by trying the non-default option.
#' @param verbose Logical. If \code{TRUE} report the AICc table.
#' @param ... Extra arguments. Not used.
#'
#' @return If \code{out = 'model'} this function returns an object of class \code{MaxEnt}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{MaxEnt} object and the data frame.
#'
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxEnt}}
#'
#' @references
#' Phillips, S.J., Anderson, R.P., Dudík, M. Schapire, R.E., and Blair, M.E.  2017.  Opening the black box: An open-source release of Maxent. \emph{Ecography} 40:887-893. \doi{10.1111/ecog.03049}
#' Warren, D.L. and S.N. Siefert. 2011. Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. \emph{Ecological Applications} 21:335-342. \doi{10.1890/10-1171.1}
#'
#' @example man/examples/trainXYZ_examples.R
#' 
#' @export

trainMaxNet <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10),
	classes = 'default',
	testClasses = TRUE,
	dropOverparam = TRUE,
	forceLinear = TRUE,
	out = 'model',
	cores = 1,
	parallelType = 'doParallel',
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########

		# response and predictors
		if (inherits(data, 'data.table')) data <- as.data.frame(data)
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		# get response and predictors
		presBg <- data[ , resp, drop=TRUE]
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
			classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
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
						regMult = thisRegMult,
						classes = classes
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

	### train models
	################
		
		work <- foreach::foreach(
			i = 1L:nrow(tuning),
			.options.multicore = mcOptions,
			.combine = 'c',
			.inorder = FALSE,
			.export = c('.trainMaxNetWorker', 'predictMaxNet'),
			.packages = c('maxnet')
		) %makeWork% {
			.trainMaxNetWorker(
				i = i,
				tuning = tuning,
				presBg = presBg,
				data = data,
				allPres = allPres,
				allBg = allBg
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
		tuning <- tuning[modelOrder, , drop = FALSE]
		models <- models[modelOrder]

		# remove models with more parameters than data points that have more than 0 parameters
		if (dropOverparam) {
		
			topModel <- models[[1L]]
			topTuning <- tuning[1L, , drop=FALSE]

			overparamModels <- which(tuning$n < tuning$numCoeff)
			if (length(overparamModels) > 0L) {
				tuning <- tuning[-overparamModels, , drop = FALSE]
				models <- models[-overparamModels]
			}

			if (length(models) == 0L & dropOverparam) {
				warning('No models had fewer coefficients than predictors. Returning best model anyway.', immediate.=TRUE)
				model <- topModel
				tuning <- topTuning
			} else if (length(models) == 0L & !dropOverparam) {
				warning('No models had fewer coefficients than predictors. NA model(s) returned.', immediate.=TRUE)
				model <- NULL
			} else {
				model <- models[[1L]]
			}
				
		}
		
		# AICc weights
		if (nrow(tuning) > 0L) {

			tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
			tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
			tuning$aicWeight <- tuning$relLike / sum(tuning$relLike)

			rownames(tuning) <- 1L:nrow(tuning)
		
		}

		if (verbose) {

			omnibus::say('')
			print(tuning, digits=3)
			omnibus::say('')

		}

	### return
	##########
		
		if (length(out) > 1L) {
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

.trainMaxNetWorker <- function(
	i,								# iterator
	tuning,							# tuning data frame
	presBg,							# vector with 1 (present) / 0 (background)
	data,							# df with all presence/background environmental data
	allPres,						# df with all presence environmental data
	allBg							# df with all background environmental data
) {
	
	thisRegMult <- tuning$regMult[i]
	thisClasses <- tuning$classes[i]
	
	# add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
	if (ncol(data) == 1L & thisClasses == 'l') {

		thisData <- data
		thisPres <- allPres
		thisBg <- allBg

		thisData$DUMMY <- rep(1, nrow(thisData))
		thisPres$DUMMY <- rep(1, nrow(allPres))
		thisBg$DUMMY <- rep(1, nrow(allBg))

	} else {

		thisData <- data
		thisPres <- allPres
		thisBg <- allBg

	}
	
	# train model
	model <- maxnet::maxnet(
		p=as.vector(presBg),
		data=thisData,
		f=maxnet::maxnet.formula(p=as.vector(presBg), data=thisData, classes=thisClasses),
		regfun=maxnet::maxnet.default.regularization,
		regmult=thisRegMult
	)

	## predict to training (and maybe test presences)
	predPres <- predictMaxNet(
		model=model,
		newdata=thisPres,
		na.rm=TRUE,
		type = 'exponential'
	)

	## predict to background
	predBg <- predictMaxNet(
		model=model,
		newdata=thisBg,
		na.rm=TRUE,
		type = 'exponential'
	)

	rawSum <- sum(c(predPres, predBg), na.rm=TRUE)

	## log likelihood
	ll <- sum(log(predPres / rawSum), na.rm=TRUE)

	## number of parameters
	numCoeff <- length(model$betas)

	# AICc
	AICc <- -2 * ll + 2 * numCoeff + (2 * numCoeff * (numCoeff + 1)) / (sum(presBg) - numCoeff - 1)

	out <- list(
		list(
			model = model,
			thisTuning = data.frame(
				regMult = thisRegMult,
				classes = thisClasses,
				numClasses = nchar(thisClasses),
				numCoeff = numCoeff,
				logLik = ll,
				AICc = AICc
			)
		)
	)
	
	out

}
