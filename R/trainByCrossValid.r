#' Calibrate a distribution/niche model using cross-validation
#'
#' This function is an extension of any of the \code{trainXYZ} functions for calibrating species distribution and ecological niche models. This function uses the \code{trainXYZ} function to calibrate and evaluate a suite of models using cross-validation. The models are evaluated against withheld data to determine the optimal settings for a "final" model using all available data. The function returns a set of models and/or a table with statistics on each model. The statistics represent various measures of model accuracy, and are calculated against training and test sites (separately).
#'
#' @param data Data frame or matrix. Response variable and environmental predictors (and no other fields) for presences and non-presence sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character vector or integer vector. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data} as predictors.
#' @param folds Either a numeric vector, or matrix or data frame which specify which rows in \code{data} belong to which folds:
#' \itemize{
#' 	\item If a vector, there must be one value per row in \code{data}. If there are \emph{K} unique values in the vector, then \emph{K} unique models will be trained. Each model will use all of the data except for rows that match a particular value in the \code{folds} vector. For example, if \code{folds = c(1, 1, 1, 2, 2, 2, 3, 3, 3)}, then three models will be trained, one with all rows that match the 2s and 3s, one with all rows matching 1s and 2s, and one will all rows matching 1s and 3s. The models will be evaluated against the training data and against the withheld data. Use \code{NA} to exclude rows from all testing/training. The default is to construct 5 folds of roughly equal size.
#' \item If a matrix or data frame, there must be one row per row in \code{data}. Each column corresponds to a different model to be trained. For a given column there should be only two unique values, plus possibly \code{NA}s. Of the two values, the lesser value will be used to identify the calibration data and the greater value the evaluation data. Rows with \code{NA}s will be ignored and not used in training or testing.  For example, a particular column could contain 1s, 2, and \code{NA}s. Data rows corresponding to 1s will be used as training data, data rows corresponding to 2s as test data, and rows with \code{NA} are dropped. The \code{NA} flag is useful for creating spatially-structured cross-validation folds where training and test sites are separated (spatially) by censored (ignored) data.
#' }
#' @param trainFx Function, name of the \code{trainXYZ} function to use. Currently the functions/algorithms supported are \code{\link[enmSdmX]{trainBrt}}, \code{\link[enmSdmX]{trainGam}}, \code{\link[enmSdmX]{trainGlm}}, \code{\link[enmSdmX]{trainMaxEnt}}, and \code{\link[enmSdmX]{trainNs}}.
#' @param ... Arguments to pass to the "trainXYZ" function.
#' @param metrics Character vector, names of evaluation metrics to calculate. If \code{weightEvalTrain} and/or \code{weightEvalTest} is \code{TRUE}, then the "train" and "test" version of each metric will be weighted versions of each. The default is to calculate all of:
#' \itemize{
#' 	\item \code{'logLoss'}: Log loss. Higher (less negative) values imply better fit.
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdmX]{evalContBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdmX]{evalAUC}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdmX]{evalTSS}}.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainSe95'} and/or \code{'trainSe90'}: Sensitivity at the threshold that ensures either 95 or 90 percent of all training presences are classified as "present" (training sensitivity = 0.95 or 0.9).
#' }
#' @param weightEvalTrain Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "train" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of sites. Note that this applies \emph{only} to the calculation of evaluation statistics, not to model calibration.  If \code{w} is supplied, they will be used for model calibration.
#' @param weightEvalTest Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "test" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied, they will be used for model calibration.
#' @param na.rm Logical, if \code{TRUE} then remove \code{NA} predictions before calculating evaluation statistics. If \code{FALSE} (default), propagate \code{NA}s (meaning if predictions contain \code{NA}s, then the evaluation statistic will most likely also be \code{NA}.)
#' @param outputModels If \code{TRUE}, then return all models (in addition to tables reporting tuning paramaeters and evaluation metrics). \emph{WARNING}: Depending on the type of model and amount of data, using \code{'models'} may produce objects that are very large in memory.
#' @param verbose Numeric. If 0 show no progress updates. If > 0 then show minimal progress updates for this function only. If > 1 show detailed progress for this function. If > 2 show detailed progress plus detailed progress for the \code{trainXYZ} function.
#'
#' @return A list object with several named elements:
#' \itemize{
#' 		\item \code{meta}: Meta-data on the model call.
#' 		\item \code{folds}: The \code{folds} object.
#' 		\item \code{models} (if \code{outputModels} is \code{TRUE}): A list of model objects, one per  data fold.
#'		\item \code{tuning}: One data frame per k-fold, each containing evaluation statistics for all candidate models in the fold.
#' }
#'
#' @references
#' Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation} 24:38-49. \doi{10.1017/S0376892997000088}
#' La Rest, K., Pinaud, D., Monestiez, P., Chadoeuf, J., and Bretagnolle, V.  2014.  Spatial leave-one-out cross-validation for variable selection in the presence of spatial autocorrelation. \emph{Global Ecology and Biogeography} 23:811-820. \doi{https://doi.org/10.1111/geb.12161}
#' Radosavljevic, A. and Anderson, R.P.  2014.  Making better Maxent models of species distributions: complexity, overfitting and evaluation.  \emph{Journal of Biogeography} 41:629-643. \doi{10.1111/jbi.12227}
#'
#' @details In some cases models do not converge (e.g., boosted regression trees and generalized additive models sometimes suffer from this issue). In this case the model will be skipped, but a data frame with the k-fold and model number in the fold will be returned in the $meta element in the output. If np models converged, then this data frame will be empty.
#'
#' @seealso \code{\link[enmSdmX]{trainBrt}}, \code{\link[enmSdmX]{trainGam}}, \code{\link[enmSdmX]{trainGlm}}, \code{\link[enmSdmX]{trainMaxEnt}}, \code{\link[enmSdmX]{trainMaxNet}}, \code{\link[enmSdmX]{trainNs}}, \code{\link[enmSdmX]{trainRf}}
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
trainByCrossValid <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	folds = dismo::kfold(data),
	trainFx = enmSdmX::trainGlm,
	...,
	metrics = c('logLoss', 'cbi', 'auc', 'tss', 'msss', 'mdss', 'minTrainPres', 'trainSe95', 'trainSe90'),
	weightEvalTrain = TRUE,
	weightEvalTest = TRUE,
	na.rm = FALSE,
	outputModels = TRUE,
	verbose = 1
) {

	ellipses <- list(...)
	hasWeights <- ('w' %in% names(ellipses))

	# response and predictors
	if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
	if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

	# get info about folds
	if (inherits(folds, c('matrix', 'data.frame'))) {

		foldsType <- 'custom'
		foldCodes <- stats::na.omit(sort(unique(c(as.matrix(folds)))))
		if (length(foldCodes) != 2) stop('"folds" matrix must have at least two unique values (aside from NAs).')
		numFolds <- ncol(folds)

		# codes for train/test sets
		trainCode <- foldCodes[1]
		testCode <- foldCodes[2]

	} else {
		foldsType <- 'simple'
		foldCodes <- stats::na.omit(sort(unique(folds)))
		numFolds <- length(foldCodes)
	}

	# create list of threshold types to be analyzed
	# have to do this because evalThresholdStats uses slightly different names
	# since it does not know if a set of presences is test/train
	threshTypes <- character()
	if ('msss' %in% metrics) threshTypes <- c(threshTypes, 'msss')
	if ('mdss' %in% metrics) threshTypes <- c(threshTypes, 'mdss')
	if ('minTrainPres' %in% metrics) threshTypes <- c(threshTypes, 'minTrainPres')
	if ('trainSe95' %in% metrics) threshTypes <- c(threshTypes, 'trainSe95')
	if ('trainSe90' %in% metrics) threshTypes <- c(threshTypes, 'trainSe90')

	### by fold
	###########

	# models and their evaluation statistics
	tuning <- models <- list()

	nonConvergedModels <- data.frame()

	for (k in 1L:numFolds) {

		if (verbose > 0) omnibus::say('Modeling k = ', k, ' on ', date(), '...', post=ifelse(verbose > 1, 2, 1), pre=ifelse(verbose > 1, 2, 0))

		### get training/testing data according to type of folds being used
		###################################################################

			# SIMPLE folds
			if (foldsType == 'simple') {

				# copy data
				thisData <- data

				# drop NAs
				nas <- which(is.na(folds))
				if (length(nas) > 0) {
					thisData <- thisData[-nas, , drop=FALSE]
					if (hasWeights) w <- w[-nas]
				}

				# train/test data
				testCode <- foldCodes[k]
				trainData <- data[folds != testCode, c(resp, preds), drop=FALSE]
				testData <- data[folds == testCode, c(resp, preds), drop=FALSE]
				if (hasWeights) {
					trainWeights <- ellipses$w[folds != testCode]
					testWeights <- ellipses$w[folds == testCode]
				} else {
					trainWeights <- rep(1, nrow(trainData))
					testWeights <- rep(1, nrow(testData))
				}

			# CUSTOM folds
			} else {

				# copy data and get folds codes
				thisData <- data[ , c(resp, preds)]
				thisFolds <- folds[ , k, drop=TRUE]

				if (hasWeights) thisWeights <- ellipses$w

				# drop NAs
				nas <- which(is.na(thisFolds))
				if (length(nas) > 0) {
					thisData <- thisData[-nas, , drop=FALSE]
					if (hasWeights) thisWeights <- thisWeights[-nas]
					thisFolds <- thisFolds[-nas]
				}

				# train/test data
				trainData <- thisData[which(thisFolds == trainCode), , drop=FALSE]
				testData <- thisData[which(thisFolds == testCode), , drop=FALSE]
				if (hasWeights) {
					trainWeights <- thisWeights[which(thisFolds == trainCode)]
					testWeights <- thisWeights[which(thisFolds == testCode)]
				} else {
					trainWeights <- rep(1, nrow(trainData))
					testWeights <- rep(1, nrow(testData))
				}

			}

			args <- utils::modifyList(ellipses, list(data=trainData, resp=resp, preds=preds, w=trainWeights, verbose=verbose > 2, out=c('models', 'tuning')))

		### train model
		###############

			thisOut <- do.call(trainFx, args=args)
			kModels <- thisOut$models
			kTuning <- thisOut$tuning

			rm(thisOut)
			gc()

			kTuning <- omnibus::insertCol(data.frame(k = rep(k, nrow(kTuning))), kTuning, at=1)

		### indices and weights for this fold
		#####################################

			# training
			whichTrainPres <- which(trainData[ , resp] == 1)
			whichTrainContrast <- which(trainData[ , resp] == 0)

			# "training" evaluation weights
			if (hasWeights & weightEvalTrain) {
				trainPresWeights <- trainWeights[whichTrainPres]
				trainContrastWeights <- trainWeights[whichTrainContrast]
			} else {
				trainPresWeights <- rep(1, length(whichTrainPres))
				trainContrastWeights <- rep(1, length(whichTrainContrast))
			}

			# testing
			whichTestPres <- which(testData[ , resp] == 1)
			whichTestContrast <- which(testData[ , resp] == 0)

			# "testing" evaluation weights
			if (hasWeights & weightEvalTest) {
				testPresWeights <- testWeights[whichTestPres]
				testContrastWeights <- testWeights[whichTestContrast]
			} else {
				testPresWeights <- rep(1, length(whichTestPres))
				testContrastWeights <- rep(1, length(whichTestContrast))
			}

			insert <- data.frame(
				numTrainPres = if (na.rm) { length(stats::na.omit(whichTrainPres)) } else { length(whichTrainPres) },
				numTrainContrast = if (na.rm) { length(stats::na.omit(whichTrainContrast)) } else { length(whichTrainContrast) },
				numTestPres = if (na.rm) { length(stats::na.omit(whichTestPres)) } else { length(whichTestPres) },
				numTestContrast = if (na.rm) { length(stats::na.omit(whichTestContrast)) } else { length(whichTestContrast) },
				trainPresWeights = sum(trainPresWeights, na.rm=na.rm),
				trainContrastWeights = sum(trainContrastWeights, na.rm=na.rm),
				testPresWeights = sum(testPresWeights, na.rm=na.rm),
				testContrastWeights = sum(testContrastWeights, na.rm=na.rm)
			)

			insert <- insert[rep(1, nrow(kTuning)), ]

			kTuning <- omnibus::insertCol(insert, kTuning, at='k', before=FALSE)

		### evaluate model vs training data
		###################################

			for (countModel in seq_along(kModels)) {

				kModel <- kModels[[countModel]]

				# if model converged
				if (is.null(kModel)) {

					nonConvergedModels <- rbind(
						nonConvergedModels,
						omnibus::insertCol(
							data.frame(modelNumber = countModel),
							into=kTuning[countModel, ],
							at='k',
							before=FALSE
						)
					)

				} else {

					# predict to training data
					predToTrain <- predictEnmSdm(model=kModel, newdata=trainData, ...)
					# predToTrain <- predictEnmSdm(model=kModel, newdata=trainData)
					predToTrainPres <- predToTrain[whichTrainPres]
					predToTrainContrast <- predToTrain[whichTrainContrast]

					# predict to testing data
					predToTest <- predictEnmSdm(model=kModel, newdata=testData, ...)
					# predToTest <- predictEnmSdm(model=kModel, newdata=testData)
					predToTestPres <- predToTest[whichTestPres]
					predToTestContrast <- predToTest[whichTestContrast]

					# log loss
					if ('logLoss' %in% metrics) {

						metricTrain <- -1 * mean(c(trainPresWeights * log(predToTrainPres), trainContrastWeights * log(1 - predToTrainContrast)), na.rm=na.rm)
						metricTest <- -1 * mean(c(testPresWeights * log(predToTestPres), testContrastWeights * log(1 - predToTestContrast)), na.rm=na.rm)

						if (countModel == 1) kTuning$logLossDelta <- kTuning$logLossTest <- kTuning$logLossTrain <- NA
						kTuning$logLossTrain[countModel] <- metricTrain
						kTuning$logLossTest[countModel] <- metricTest
						kTuning$logLossDelta[countModel] <- metricTrain - metricTest

					}

					# CBI
					if ('cbi' %in% metrics) {

						# metricTrain <- evalContBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
						# metricTest <- evalContBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)

						metricTrain <- evalContBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
						metricTest <- evalContBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)

						if (countModel == 1) kTuning$cbiDelta <- kTuning$cbiTest <- kTuning$cbiTrain <- NA
						kTuning$cbiTrain[countModel] <- metricTrain
						kTuning$cbiTest[countModel] <- metricTest
						kTuning$cbiDelta[countModel] <- metricTrain - metricTest

					}

					# AUC
					if ('auc' %in% metrics) {

						metricTrain <- evalAUC(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
						metricTest <- evalAUC(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)

						if (countModel == 1) kTuning$aucDelta <- kTuning$aucTest <- kTuning$aucTrain <- NA
						kTuning$aucTrain[countModel] <- metricTrain
						kTuning$aucTest[countModel] <- metricTest
						kTuning$aucDelta[countModel] <- metricTrain - metricTest

					}

					# TSS
					if ('tss' %in% metrics) {

						# metricTrain <- evalTSS(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
						# metricTest <- evalTSS(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)

						metricTrain <- evalTSS(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
						metricTest <- evalTSS(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)

						if (countModel == 1) kTuning$tssDelta <- kTuning$tssTest <- kTuning$tssTrain <- NA
						kTuning$tssTrain[countModel] <- metricTrain
						kTuning$tssTest[countModel] <- metricTest
						kTuning$tssDelta[countModel] <- metricTrain - metricTest

					}

					# threshold-dependent
					if (length(threshTypes) > 0) {

						for (thisThreshType in threshTypes) {

							# threshold
							if (thisThreshType == 'minTrainPres') {
								threshCode <- 'minPres'
								sensitivity <- 0
							} else if (thisThreshType == 'trainSe95') {
								threshCode <- 'sensitivity'
								sensitivity <- 0.95
							} else if (thisThreshType == 'trainSe90') {
								threshCode <- 'sensitivity'
								sensitivity <- 0.9
							} else {
								threshCode <- thisThreshType
								sensitivity <- 0
							}

							# get thresholds
							# thresholds <- evalThreshold(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
							thresholds <- evalThreshold(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, sensitivity=sensitivity, na.rm=na.rm)

							thold <- thresholds[[threshCode]]

							# evaluate
							trainEval <- evalThresholdStats(thold, pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
							testEval <- evalThresholdStats(thold, pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)

							# remember
							if (countModel == 1) {

								kTuning$TEMP1 <- kTuning$TEMP2 <- kTuning$TEMP3 <- kTuning$TEMP4 <- kTuning$TEMP5 <- kTuning$TEMP6 <- kTuning$TEMP7 <- NA
								names(kTuning)[(ncol(kTuning) - 6):ncol(kTuning)] <- paste0(thisThreshType, c('Thold', 'SeTrain', 'SeTest', 'SeDelta', 'SpTrain', 'SpTest', 'SpDelta'))

							}

							kTuning[countModel, paste0(thisThreshType, 'Thold')] <- thold
							kTuning[countModel, paste0(thisThreshType, 'SeTrain')] <- trainEval[['sensitivity']]
							kTuning[countModel, paste0(thisThreshType, 'SeTest')] <- testEval[['sensitivity']]
							kTuning[countModel, paste0(thisThreshType, 'SeDelta')] <- trainEval[['sensitivity']] - testEval[['sensitivity']]
							kTuning[countModel, paste0(thisThreshType, 'SpTrain')] <- trainEval[['specificity']]
							kTuning[countModel, paste0(thisThreshType, 'SpTest')] <- testEval[['specificity']]
							kTuning[countModel, paste0(thisThreshType, 'SpDelta')] <- trainEval[['specificity']] - testEval[['specificity']]

						} # next threshold-dependent evaluation metric!

					} else { # threshold-specific metrics
						sensitivity <- NA
					}

				} # if model converged

			} # next model in this k-fold

		tuning[[k]] <- kTuning
		if (outputModels) models[[k]] <- kModels
		gc()

	} # next fold

	### return
	##########
	output <- list()

	# meta data
	meta <- list(
		resp = resp,
		preds = preds,
		metrics = metrics,
		sensitivity = sensitivity,
		weightEvalTrain = weightEvalTrain,
		weightEvalTest = weightEvalTest,
		na.rm = na.rm,
		outputModels = outputModels,
		nonConvergedModels = nonConvergedModels
	)

	if (length(ellipses) > 0) {
		for (i in seq_along(ellipses)) {
			meta <- c(meta, ellipses[i])
		}
	}

	output$meta <- meta
	output$folds <- folds

	# models and tuning
	if (outputModels) output$models <- models
	output$tuning <- tuning

	# class(output) <- c(class(output), 'crossValid')
	output

}
