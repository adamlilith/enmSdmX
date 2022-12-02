#' Summarize distribution/niche model cross-validation object
#'
#' This function summarizes models calibrated using the \code{\link[enmSdmX]{trainByCrossValid}} function. It returns aspects of the best models across k-folds (the particular aspects depends on the kind of models used).
#' @param x The output from the \code{\link{trainCrossValid}} function (which is a list). Note that the object \emph{must} include a sublist named \code{tuning}.
#' @param trainFxName Character, name of function used to train the SDM (examples: \code{'trainGlm'}, \code{'trainMaxEnt'}, \code{'trainBrt'}).
#' @param metric Metric by which to select the best model in each k-fold. This can be any of the columns that appear in the data frames in \code{x$tuning} (or any columns added manually), but typically is one of the following \emph{plus} either \code{Train}, \code{Test}, or \code{Delta} (e.g., \code{'logLossTrain'}, \code{'logLossTest'}, or \code{'logLossDelta'}):
#' \itemize{
#' 	\item \code{'logLoss'}: Log loss.
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdmX]{evalContBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdmX]{evalAUC}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdmX]{evalTSS}}.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainSe95'} and/or \code{'trainSe90'}: Sensitivity at the threshold that ensures either 95% or 90% of all training presences are classified as "present" (training sensitivity = 0.95 or 0.9).
#' }
#' @param decreasing Logical, if \code{TRUE} (default), for each k-fold sort models by the value listed in \code{metric} in decreasing order (highest connotes "best", lowest "worst"). If \code{FALSE} use the lowest value of \code{metric}.
#' @return Data frame with statistics on the best set of models across k-folds. Depending on the model algorithm, this could be:
#' \itemize{
#' 	\item BRTs (boosted regression trees): Learning rate, tree complexity, and bag fraction.
#' 	\item GLMs (generalized linear models): Frequency of use of each term in the best models.
#' 	\item Maxent: Frequency of times each specific combination of feature classes was used in the best models plus mean master regularization multiplier for each feature set.
#' 	\item NSs (natural splines): Data frame, one row per fold and one column per predictor, with values representing the maximum degrees of freedom used for each variable in the best model of each fold.
#' 	\item RFs (random forests): Data frame, one row per fold, with values representing the optimal value of \code{mtry} (see \code{\link[randomForest]{randomForest}}).
#' }
#' @seealso \code{\link[enmSdmX]{trainByCrossValid}}, \code{\link[enmSdmX]{trainBrt}}, \code{\link[enmSdmX]{trainGam}}, \code{\link[enmSdmX]{trainGlm}}, \code{\link[enmSdmX]{trainMaxEnt}}, \code{\link[enmSdmX]{trainNs}}, \code{\link[enmSdmX]{trainRf}}
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
#' # create background sites (using just 1000 to speed things up!)
#' bgEnv <- terra::spatSample(madEnv, 3000)
#' bgEnv <- bgEnv[complete.cases(bgEnv), ]
#' bgEnv <- bgEnv[sample(nrow(bgEnv), 1000), ]
#' 
#' # collate occurrences and background sites
#' presBg <- data.frame(
#'    presBg = c(
#'       rep(1, nrow(occEnv)),
#'       rep(0, nrow(bgEnv))
#'    )
#' )
#' 
#' env <- rbind(occEnv, bgEnv)
#' env <- cbind(presBg, env)
#' 
#' predictors <- c('bio1', 'bio12')
#' 
#' # using "vector" form of "folds" argument
#' folds <- dismo::kfold(env, 3) # just 3 folds (for speed)
#' 
#' ## MaxEnt
#' mxx <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainMaxEnt,
#' 	regMult = 1:2 # too few values for valid model, but fast!
#' )
#' 
#' ## generalized linear models
#' glx <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainGlm
#' )
#' 
#' ## generalized linear models
#' gax <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainGam
#' )
#' 
#' ## natural splines
#' natx <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainNs
#' )
#' 
#' ## boosted regression trees
#' brtx <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainBrt,
#' 	learningRate = 0.001, # too few values for reliable model(?)
#' 	treeComplexity = 2, # too few values for reliable model, but fast
#' 	minTrees = 1200, # minimum trees for reliable model(?), but fast
#' 	maxTrees = 1200, # too small for reliable model(?), but fast
#' 	tryBy = 'treeComplexity',
#' 	anyway = TRUE # return models that did not converge
#' )
#' 
#' 
#' ## random forests
#' rfx <- trainByCrossValid(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = c('bio1', 'bio12'),
#' 	folds = folds,
#' 	trainFx = trainRf
#' )
#' 
#' # summarize MaxEnt feature sets and regularization across folds
#' summaryByCrossValid(mxx, trainFxName = 'trainMaxEnt')
#' 
#' # summarize GLM terms across folds
#' summaryByCrossValid(glx, trainFxName = 'trainGlm')
#' 
#' # summarize GAM terms across folds
#' summaryByCrossValid(gax, trainFxName = 'trainGam')
#' 
#' # summarize natural splines terms across folds
#' summaryByCrossValid(natx, trainFxName = 'trainNs')
#' 
#' # summarize BRT parameters across folds
#' # Note that to get BRTs to run fast we allowed no variation
#' # so the summary in this example is fairly boring.
#' summaryByCrossValid(brtx, trainFxName = 'trainBrt')
#' 
#' # summarize random forests 'mtry' parameter across folds
#' summaryByCrossValid(natx, trainFxName = 'trainRf')
#' 
#' @export
summaryByCrossValid <- function(
	x,
	trainFxName = 'trainGlm',
	metric = 'cbiTest',
	decreasing = TRUE
) {

	trainFxName <- tolower(trainFxName)

	tuning <- x$tuning

	### order each tuning attempt by performance
	for (k in seq_along(tuning)) {

		tuningOrder <- order(tuning[[k]][ , metric], decreasing=decreasing)
		tuning[[k]] <- tuning[[k]][tuningOrder, ]

	}

	### BRT
	if (trainFxName == tolower('trainBrt')) {

		# get list of terms in best models
		params <- data.frame()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]
			if (thisTuning$converged[1]) {

				params <- rbind(
					params,
					data.frame(
						learningRate = thisTuning$learningRate[1],
						treeComplexity = thisTuning$treeComplexity[1],
						bagFraction = thisTuning$bagFraction[1],
						nTrees = thisTuning$nTrees[1]
					)
				)

			}

		}

		# summarize best models
		if (nrow(params) > 0) {

			out <- rbind(
				apply(params, 2, min, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.25, na.rm=TRUE),
				apply(params, 2, mean, na.rm=TRUE),
				apply(params, 2, stats::median, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.75, na.rm=TRUE),
				apply(params, 2, max, na.rm=TRUE)
			)

			out <- as.data.frame(out)
			out$treeComplexity <- max(1, round(out$treeComplexity))
			out$nTrees <- round(out$nTrees)

			value <- data.frame(value=c('min', '25th percent quantile', 'mean', 'median', '75th percent quantile', 'max'))
			out <- omnibus::insertCol(value, out, 1)

		} else {
			out <- data.frame()
		}


	} else if (trainFxName == tolower('trainGam')) {

		# get list of terms in best models
		term <- character()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1L]
			thisTerm <- strsplit(thisTerm, ' + ', fixed=TRUE)[[1L]]
			thisTerm <- thisTerm[2L:length(thisTerm)]
			if (any(thisTerm == 'presBg')) thisTerm <- thisTerm[-which(thisTerm == 'presBg')]
			if (any(thisTerm == '~')) thisTerm <- thisTerm[-which(thisTerm == '~')]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			if (length(thisTerm) == 0) thisTerm <- 1

			term <- c(term, thisTerm)

		}
		
		term <- unique(term)
		out <- data.frame(term)
		out$frequencyInBestModels <- 0

		# tally usage of best term(s)
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1]
			thisTerm <- strsplit(thisTerm, ' + ', fixed=TRUE)[[1L]]
			thisTerm <- thisTerm[2L:length(thisTerm)]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			for (countTerm in seq_along(thisTerm)) {
				whichTerm <- which(out$term == thisTerm[countTerm])
				out$frequencyInBestModels[whichTerm] <- out$frequencyInBestModels[whichTerm] + 1
			}

		}

		out <- out[order(out$frequencyInBestModels, decreasing = TRUE), ]
		out$proportionOfModels <- out$frequencyInBestModels / length(x$tuning)

	} else if (trainFxName == tolower('trainGlm')) {

		# get list of terms in best models
		term <- character()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1]
			thisTerm <- strsplit(thisTerm, ' ')[[1]]
			thisTerm <- thisTerm[3:length(thisTerm)]
			if (any(thisTerm == 'presBg')) thisTerm <- thisTerm[-which(thisTerm == 'presBg')]
			if (any(thisTerm == '~')) thisTerm <- thisTerm[-which(thisTerm == '~')]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			if (length(thisTerm) == 0) thisTerm <- 1

			term <- c(term, thisTerm)

		}

		term <- unique(term)
		out <- data.frame(term)
		out$frequencyInBestModels <- 0

		# tally usage of best term(s)
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1]
			thisTerm <- strsplit(thisTerm, ' ')[[1]]
			thisTerm <- thisTerm[3:length(thisTerm)]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			for (countTerm in seq_along(thisTerm)) {
				whichTerm <- which(out$term == thisTerm[countTerm])
				out$frequencyInBestModels[whichTerm] <- out$frequencyInBestModels[whichTerm] + 1
			}

		}

		out <- out[order(out$frequencyInBestModels, decreasing = TRUE), ]
		out$proportionOfModels <- out$frequencyInBestModels / length(x$tuning)

	### MAXENT
	##########

	} else if (trainFxName == tolower('trainMaxEnt')) {

		## tally frequency of feature class combinations across best models and mean regularization associated with each set of features

		# make vector of all possible feature combinations
		simpleFeats <- c('l', 'p', 'q', 'h')
		featCombos <- vector('list', length(simpleFeats))
		for (i in seq_along(simpleFeats)) {
			featCombos[[i]] <- utils::combn(simpleFeats, i, paste, collapse = '')
		}

		featCombos <- unique(unlist(featCombos))

		featFreqs <- featRegMultSums <- rep(0, length(featCombos))
		names(featFreqs) <- names(featRegMultSums) <- featCombos

		# tally number of times each feature combination used in top models
		targetFeatLengths <- nchar(featCombos)
		for (k in seq_along(tuning)) {

			thisModelFeats <- tuning[[k]]$classes[1]
			thisModelFeats <- strsplit(thisModelFeats, '')[[1]]
			modelFeatLength <- length(thisModelFeats)

			modelFeatInTarget <- rep(TRUE, length(featCombos)) # flags if potential features are *all* used in model
			names(modelFeatInTarget) <- featCombos
			for (countModelFeat in seq_along(thisModelFeats)) {

				thisModelFeatInTarget <- grepl(pattern=thisModelFeats[countModelFeat], featCombos)
				modelFeatInTarget <- (modelFeatInTarget & thisModelFeatInTarget)

			}

			modelFeatInTarget <- modelFeatInTarget & (modelFeatLength == targetFeatLengths)

			featFreqs[modelFeatInTarget] <- featFreqs[modelFeatInTarget] + 1
			featRegMultSums[modelFeatInTarget] <- featRegMultSums[modelFeatInTarget] + tuning[[k]]$regMult[1]

		}

		featRegMultMeans <- featRegMultSums / featFreqs
		# featFreqs <- featFreqs / sum(featFreqs)

		out <- data.frame(
			featureSet = featCombos,
			frequencyInBestModels = featFreqs,
			meanRegMult = featRegMultMeans
		)

	### NS
	######

	} else if (trainFxName == tolower('trainNs')) {

		# get predictor names
		strings <- colnames(tuning[[1]])
		strings <- strings[grepl(strings, pattern='splines::ns')]
		strings <- strsplit(strings, 'splines::ns')
		preds <- character()
		for (i in seq_along(strings)) {
			string <- strings[[i]][2]
			pred <- substr(string, 2, regexpr(string, pattern=' df') - 2)
			preds <- c(preds, pred)
		}

		preds <- sort(unique(preds))

		# create table to hold information on degrees of freedom for each top model
		subOut <- data.frame(DUMMY=NA)
		colnames(subOut) <- preds[1]
		if (length(preds) > 1) {
			for (i in 2L:length(preds)) {
				subOut$DUMMY <- NA
				colnames(subOut)[ncol(subOut)] <- preds[i]
			}
		}
		subOut <- subOut[rep(1, length(tuning)), , drop=FALSE]

		# find df for each predictor
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			for (pred in preds) {

				colPattern <- paste0('splines::ns\\(', pred, ', df')
				cols <- colnames(thisTuning)[grepl(colnames(thisTuning), pattern=colPattern)]

				# get degrees of freedom used to evaluate this variable
				dfs <- numeric()
				for (countCol in seq_along(cols)) {
					before <- paste0('splines::ns\\(', pred, ', df')
					thisDf <- strsplit(cols[countCol], before)[[1]][2]
					thisDf <- strsplit(thisDf, ' ')[[1]][3]
					thisDf <- substr(thisDf, 1, nchar(thisDf) - 1)
					thisDf <- as.numeric(thisDf)
					dfs <- c(dfs, thisDf)
				}

				vals <- thisTuning[1, cols, drop=FALSE]
				bestDf <- c(!apply(vals, 1, FUN=is.na))
				bestDf <- if (all(!bestDf)) {
					NA
				} else {
					max(dfs[bestDf])
				}

				subOut[k, pred] <- bestDf

			} # next predictor

		} # next model

		out <- data.frame()
		for (pred in preds) {
			thisPredOut <- data.frame(
				term=pred,
				frequencyInBestModels=sum(!is.na(subOut[ , pred])),
				minDf=min(subOut[ , pred], na.rm=TRUE),
				quantile25thDf=stats::quantile(subOut[ , pred], 0.25, na.rm=TRUE),
				meanDf=mean(subOut[ , pred], na.rm=TRUE),
				medianDf=stats::median(subOut[ , pred], na.rm=TRUE),
				quantile75thDf=stats::quantile(subOut[ , pred], 0.75, na.rm=TRUE),
				maxDf=max(subOut[ , pred], na.rm=TRUE)
			)

			out <- rbind(out, thisPredOut)
		
		}

	} else if (trainFxName == tolower('trainRf')) {
	
		# get predictor names
		strings <- colnames(tuning[[1]])
		out <- data.frame()
		for (k in seq_along(tuning)) {
		
			out <- rbind(
				out,
				data.frame(
					mtry = tuning[[k]]$mtry[1L],
					metric = tuning[[k]][1L, metric]
				)
			)
		
		}
		
		names(out)[2L] <- metric
		
		
		mtries[k] <- tuning[[k]]$mtry
		
	}

	rownames(out) <- NULL
	out

}
