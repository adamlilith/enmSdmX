#' Summarize distribution/niche model cross-validation object
#'
#' This function summarizes models calibrated using the \code{\link[enmSdmX]{trainByCrossValid}} function. It returns aspects of the best models across k-folds (the particular aspects depends on the kind of models used).
#' @param x The output from the \code{\link{trainByCrossValid}} function (which is a list). Note that the object \emph{must} include a sublist named \code{tuning}.
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
#' @param interceptOnly Logical. If \code{TRUE} (default) and the top models in each case were intercept-only models, return an emppty data frame (with a warning). If \code{FALSE}, return results using the first model in each fold that was not an intercept-only model. This is only used if the training function was a generalized linear model (GLM), natural splines model (NS), or generalized additive model (GAM).
#'
#' @return Data frame with statistics on the best set of models across k-folds. Depending on the model algorithm, this could be:
#' \itemize{
#' 	\item BRTs (boosted regression trees): Learning rate, tree complexity, and bag fraction.
#' 	\item GLMs (generalized linear models): Frequency of use of each term in the best models.
#' 	\item Maxent: Frequency of times each specific combination of feature classes was used in the best models plus mean master regularization multiplier for each feature set.
#' 	\item NSs (natural splines): Data frame, one row per fold and one column per predictor, with values representing the maximum degrees of freedom used for each variable in the best model of each fold.
#' 	\item RFs (random forests): Data frame, one row per fold, with values representing the optimal value of \code{numTrees} and \code{mtry} (see \code{\link[ranger]{ranger}}).
#' }
#' @seealso \code{\link[enmSdmX]{trainByCrossValid}}, \code{\link[enmSdmX]{trainBRT}}, \code{\link[enmSdmX]{trainGAM}}, \code{\link[enmSdmX]{trainGLM}}, \code{\link[enmSdmX]{trainMaxEnt}}, \code{\link[enmSdmX]{trainNS}}, \code{\link[enmSdmX]{trainRF}}
#'
#' @example man/examples/trainByCrossValid_examples.r
#'
#' @export
summaryByCrossValid <- function(
	x,
	metric = 'cbiTest',
	decreasing = TRUE,
	interceptOnly = TRUE
) {

	trainFxName <- tolower(x$meta$trainFxName)
	tuning <- x$tuning

	### order each tuning attempt by performance
	for (k in seq_along(tuning)) {

		tuningOrder <- order(tuning[[k]][ , metric], decreasing=decreasing)
		tuning[[k]] <- tuning[[k]][tuningOrder, ]

	}

	### BRT
	if (trainFxName == tolower('trainBRT')) {

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
			out$treeComplexity <- pmax(1, round(out$treeComplexity))
			out$nTrees <- round(out$nTrees)

			value <- data.frame(value=c('min', '25th percent quantile', 'mean', 'median', '75th percent quantile', 'max'))
			out <- omnibus::insertCol(value, out, 1)

		} else {
			out <- data.frame()
		}


	} else if (trainFxName == tolower('trainGAM')) {

		# get list of terms in best models
		term <- character()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			foundViableModel <- FALSE
			while (foundViableModel) {

				thisTerm <- thisTuning$model[1L]
				thisTerm <- strsplit(thisTerm, ' + ', fixed=TRUE)[[1L]]
				thisTerm <- thisTerm[2L:length(thisTerm)]
				if (any(thisTerm == 'presBg')) thisTerm <- thisTerm[-which(thisTerm == 'presBg')]
				if (any(thisTerm == '~')) thisTerm <- thisTerm[-which(thisTerm == '~')]
				if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
				if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
				if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

				if (length(thisTerm) == 0 & interceptOnly) {
					thisTerm <- '1'
					foundViableModel <- TRUE
				} else if (length(thisTerm) == 0 & !interceptOnly) {
					thisTuning <- thisTuning[-1, drop = FALSE]
				} else {
					foundViableModel <- TRUE
				}

			}


			term <- c(term, thisTerm)

		} # next fold

		term <- unique(term)

		if (all(term == '1')) {

			warning('All top models were intercept-only.')
			out <- data.frame()
		} else {

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
			out <- out[order(out$frequencyInBestModels, decreasing = TRUE), ]

		} # if at least one top model wasn't intercept-only

	} else if (trainFxName == tolower('trainGLM')) {

		# get list of terms in models
		term <- character()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			foundViableModel <- FALSE
			while (!foundViableModel) {

				thisTerm <- thisTuning$model[1L]
				thisTerm <- strsplit(thisTerm, '~')[[1L]]
				thisTerm <- strsplit(thisTerm, ' \\+ ')[[1L]]
				if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]

				if (length(thisTerm) == 0 & interceptOnly) {
					thisTerm <- '1'
					foundViableModel <- TRUE
				} else if (length(thisTerm) == 0 & !interceptOnly) {
					thisTuning <- thisTuning[-1, drop = FALSE]
				} else {
					foundViableModel <- TRUE
				}

			}

			term <- c(term, thisTerm)

		} # next fold

		term <- unique(term)

		if (all(term == '1')) {
			warning('All top models were intercept-only.')
			out <- data.frame()
		} else {

			out <- data.frame(term)
			out$frequencyInBestModels <- 0

			# tally usage of best term(s)
			for (k in seq_along(tuning)) {

				thisTuning <- tuning[[k]]

				thisTerm <- thisTuning$model[1L]

				# natural splines model
				if (grepl(thisTerm, pattern='splines')) {

					thisTerm <- strsplit(thisTerm, ' \\+\\ ')[[1L]]
					if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]

				# "normal" GLM model
				} else  {

					thisTerm <- strsplit(thisTerm, ' ')[[1]]
					thisTerm <- thisTerm[3:length(thisTerm)]
					if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
					if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
					if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

				}

				for (countTerm in seq_along(thisTerm)) {
					if (!is.na(thisTerm[countTerm])) {
						whichTerm <- which(out$term == thisTerm[countTerm])
						out$frequencyInBestModels[whichTerm] <- out$frequencyInBestModels[whichTerm] + 1
					}
				}

			}

			out <- out[order(out$frequencyInBestModels, decreasing = TRUE), ]
			out$proportionOfModels <- out$frequencyInBestModels / length(x$tuning)

		} # at least one top model was not intercept-only

	### MAXENT
	##########

	} else if (trainFxName %in% tolower(c('trainMaxEnt', 'trainMaxNet'))) {

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

		out <- data.frame(
			featureSet = featCombos,
			frequencyInBestModels = featFreqs,
			meanRegMult = featRegMultMeans
		)

		out <- out[order(out$frequencyInBestModels, decreasing = TRUE), ]


	# ### NS
	# ######

	# } else if (trainFxName == tolower('trainNS')) {

		# # get predictor names
		# strings <- colnames(tuning[[1]])
		# strings <- strings[grepl(strings, pattern='splines::ns')]
		# strings <- strsplit(strings, 'splines::ns')
		# preds <- character()
		# for (i in seq_along(strings)) {
			# string <- strings[[i]][2]
			# pred <- substr(string, 2, regexpr(string, pattern=' df') - 2)
			# preds <- c(preds, pred)
		# }

		# preds <- sort(unique(preds))

		# # create table to hold information on degrees of freedom for each top model
		# subOut <- data.frame(DUMMY=NA)
		# colnames(subOut) <- preds[1]
		# if (length(preds) > 1) {
			# for (i in 2L:length(preds)) {
				# subOut$DUMMY <- NA
				# colnames(subOut)[ncol(subOut)] <- preds[i]
			# }
		# }
		# subOut <- subOut[rep(1, length(tuning)), , drop=FALSE]

		# # find df for each predictor
		# for (k in seq_along(tuning)) {

			# thisTuning <- tuning[[k]]

			# for (pred in preds) {

				# colPattern <- paste0('splines::ns\\(', pred, ', df')
				# cols <- colnames(thisTuning)[grepl(colnames(thisTuning), pattern=colPattern)]

				# # get degrees of freedom used to evaluate this variable
				# dfs <- numeric()
				# for (countCol in seq_along(cols)) {
					# before <- paste0('splines::ns\\(', pred, ', df')
					# thisDf <- strsplit(cols[countCol], before)[[1]][2]
					# thisDf <- strsplit(thisDf, ' ')[[1]][3]
					# thisDf <- substr(thisDf, 1, nchar(thisDf) - 1)
					# thisDf <- as.numeric(thisDf)
					# dfs <- c(dfs, thisDf)
				# }

				# vals <- thisTuning[1, cols, drop=FALSE]
				# bestDf <- c(!apply(vals, 1, FUN=is.na))
				# bestDf <- if (all(!bestDf)) {
					# NA
				# } else {
					# max(dfs[bestDf])
				# }

				# subOut[k, pred] <- bestDf

			# } # next predictor

		# } # next model

		# out <- data.frame()
		# for (pred in preds) {
			# thisPredOut <- data.frame(
				# term=pred,
				# frequencyInBestModels=sum(!is.na(subOut[ , pred])),
				# minDf=min(subOut[ , pred], na.rm=TRUE),
				# quantile25thDf=stats::quantile(subOut[ , pred], 0.25, na.rm=TRUE),
				# meanDf=mean(subOut[ , pred], na.rm=TRUE),
				# medianDf=stats::median(subOut[ , pred], na.rm=TRUE),
				# quantile75thDf=stats::quantile(subOut[ , pred], 0.75, na.rm=TRUE),
				# maxDf=max(subOut[ , pred], na.rm=TRUE)
			# )

			# out <- rbind(out, thisPredOut)

		# }

	} else if (trainFxName == tolower('trainRF')) {

		# get list of terms in best models
		params <- data.frame()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			params <- rbind(
				params,
				data.frame(
					numTrees = thisTuning$numTrees[1L],
					mtry = thisTuning$mtry[1L],
					oobError = thisTuning$oobError[1L]
				)
			)

		}

		# summarize best models
		if (nrow(params) > 0L) {

			out <- rbind(
				apply(params, 2, min, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.25, na.rm=TRUE),
				apply(params, 2, mean, na.rm=TRUE),
				apply(params, 2, stats::median, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.75, na.rm=TRUE),
				apply(params, 2, max, na.rm=TRUE)
			)

			out <- as.data.frame(out)
			out$mtry <- pmax(1, round(out$mtry))
			out$numTrees <- round(out$numTrees)

			value <- data.frame(value=c('min', '25th percent quantile', 'mean', 'median', '75th percent quantile', 'max'))
			out <- omnibus::insertCol(value, out, 1)

		} else {
			out <- data.frame()
		}


	}

	rownames(out) <- NULL
	out

}
