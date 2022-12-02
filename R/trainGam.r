#' Calibrate a generalized additive model (GAM)
#'
#' This function constructs a GAM piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}).
#' @param gamma Initial penalty to degrees of freedom to use (larger ==> smoother fits).
#' @param construct Logical. If \code{TRUE} then construct model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{initialTerms}.
#' @param select Logical. If \code{TRUE} then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only if \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is \code{TRUE}.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{dredge()} is 31, so if this number is >31 and \code{select} is \code{TRUE} then it is forced to 31 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param interaction Character or \code{NULL}. Type of interaction term to use (\code{te}, \code{ts}, \code{s}, etc.). See \code{?te} (for example) for help on any one of these. If \code{NULL} then interactions are not used.
#' @param w Either logical in which case TRUE causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by AICc.
#' }
#' @param verbose Logical. If TRUE then display intermediate results on the display device.
#' @param ... Extra arguments (not used).
#' @return If \code{out = 'model'} this function returns an object of class \code{gam}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gam} object and the data frame.
#' @seealso \code{\link[mgcv]{gam}}
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


trainGam <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	gamma = 1,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 10,
	initialTerms = 8,
	interaction = 'te',
	w = TRUE,
	out = 'model',
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########

		# ellipses <- list(...)

		# force number of starting terms to 31 or less
		if (select & initialTerms > 31) {
			initialTerms <- 31
			warning('initialTerms must be 31 or less. Forcing to 31.')
		}

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

	#############
	## weights ##
	#############

		# model weights
		if (is.logical(w)) {
			if (w && (family %in% c('binomial', 'quasibinomial'))) {
				posCases <- sum(data[ , resp, drop=TRUE] == 1)
				negCases <- sum(data[ , resp, drop=TRUE] == 0)
				w <- c(rep(1, posCases), rep(posCases / (posCases + negCases), negCases))
			} else {
				w <- rep(1, nrow(data))
			}
		} else if (inherits(w, 'character')) {
			w <- data[ , w, drop=TRUE]
		}
		w <- w / max(w)

	################################
	## initial model construction ##
	################################

		# create starting formula
		form <- paste0(resp, ' ~ 1')

		if (construct) {

			### SINGLE-variable terms
			for (thisPred in preds) { # for each predictor test single-variable terms

				term <- if (class(data[ , thisPred]) != 'factor') {
					paste0('s(', thisPred, ', bs=\'cs\')')
				} else {
					thisPred
				}

				thisThisForm <- paste0(form, ' + ', term)

				thisAic <- MuMIn::AICc(
					mgcv::gam(
						formula=stats::as.formula(thisThisForm),
						family=family,
						data=data,
						method='ML',
						optimizer=c('outer', 'newton'),
						scale=-1,
						select=TRUE,
						gamma=gamma,
						weights=w,
						na.action='na.fail',
						...
					)
				)

				# remember
				gamFrame <- if (exists('gamFrame', inherits=FALSE)) {
					rbind(gamFrame, data.frame(term=term, AICc=thisAic))
				} else {
					data.frame(term=term, AICc=thisAic)
				}

			} # next single-variable term

			### TWO-variable terms
			if (length(preds) > 1 & !is.null(interaction)) {
				
				for (thisPred in preds[1:(length(preds)-1)]) { # for each predictor test two-variable terms

					for (thatPred in preds[ (which(preds==thisPred) + 1):length(preds) ]) { # for each second predictor test two-variable terms

						# create term
						term <- if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) != 'factor') {

							term <- paste0(interaction, '(', thisPred, ', ', thatPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) != 'factor') {

							paste0(interaction, '(', thatPred, ', by=', thisPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) == 'factor') {

							paste0(interaction, '(', thisPred, ', by=', thatPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) == 'factor') {

							paste0(thisPred, ' * ', thatPred)

						}

						thisAic <- MuMIn::AICc(
							mgcv::gam(
								formula=stats::as.formula(paste0(form, ' + ', term)),
								family=family,
								data=data,
								method='ML',
								optimizer=c('outer', 'newton'),
								scale=-1,
								select=TRUE,
								gamma=gamma,
								weights=w,
								na.action='na.fail',
								...
							)
						)

						# remember
						gamFrame <- if (exists('gamFrame', inherits=FALSE)) {
							rbind(gamFrame, data.frame(term=term, AICc=thisAic))
						} else {
							data.frame(term=term, AICc=thisAic)
						}

					}  # for each second predictor test two-variable terms

				} # for each predictor test two-variable terms
				
			} # if interactions

			# sort by AIC
			gamFrame <- gamFrame[order(gamFrame$AICc), ]

			# print AICc frame
			if (verbose) {

				omnibus::say('GAM construction results for each term tested:', level=2)
				print(gamFrame)
				omnibus::say('')

			}

			### no model selection
			######################
			if (!select) {

				## construct final model
				form <- paste0(form, ' + ', gamFrame$term[1]) # add first term

				# for each set of presences > min num required, add a term
				if (floor(sum(data[ , resp]) / presPerTermInitial ) - 1 > 1 & initialTerms > 1) {

					termsToAdd <- 2:min(initialTerms, c(floor(sum(data[ , resp]) / presPerTermInitial ) - 1, nrow(gamFrame) - 1))

					form <- paste0(form, ' + ', paste0(gamFrame$term[termsToAdd], collapse=' + '))

				} # if there are sufficient presences for additional terms beyond first
				
				# train FULL model
				model <- mgcv::gam(
					formula=stats::as.formula(form),
					family=family,
					data=data,
					method='ML',
					optimizer=c('outer', 'newton'),
					select=TRUE,
					gamma=gamma,
					weights=w,
					...
				)
				

			### model selection
			###################
			
			} else if (select) {
			
				if (verbose) omnibus::say('Selecting best model:', level=2)

				# create grid of model terms
				n <- nrow(gamFrame)
				gamFrame$minPres <- seq(presPerTermFinal, presPerTermFinal * n, by=presPerTermFinal)
				
				numPres <- sum(data[ , resp, drop=TRUE])
				totalTerms <- max(which(gamFrame$minPres <= numPres))
				gamFrame <- gamFrame[1L:totalTerms, , drop=FALSE]
				
				grid <- list()
				for (i in 1L:nrow(gamFrame)) grid <- c(grid, list(yesNo = c(FALSE, TRUE)))
				grid <- expand.grid(grid)
				colnames(grid) <- gamFrame$term
			
				### do all models
				bestAicc <- Inf
				tuning <- data.frame()
				if ('models' %in% out) models <- list()
				for (countModel in 1L:nrow(grid)) {
				
					form <- paste0(resp, ' ~ 1')
					thisGridRow <- unlist(grid[countModel, ])
					if (sum(thisGridRow) > 0) form <- paste(form, '+', paste(colnames(grid)[thisGridRow], collapse = ' + '))
					
					if (verbose) say('   Evaluating: ', form)
					
					thisModel <- mgcv::gam(
						formula=stats::as.formula(form),
						family=family,
						data=data,
						method='ML',
						optimizer=c('outer', 'newton'),
						select=TRUE,
						gamma=gamma,
						weights=w#,
						# ...
					)

					# tuning table
					thisAicc <- MuMIn::AICc(thisModel)

					dev <- thisModel$deviance
					nullDev <- thisModel$null.deviance
					devExplained <- (nullDev - dev) / nullDev
				
					tuning <- rbind(
						tuning,
						data.frame(
							model = form,
							AICc = thisAicc,
							deviance = dev,
							devianceExplained = devExplained
						)
					)
					
					# remember best model
					if ('model' %in% out) {
						if (thisAicc < bestAicc) {
							model <- thisModel
							bestAicc <- thisAicc
						}
					}
					
					# remember model
					if ('models' %in% out) models[[length(models) + 1L]] <- thisModel
				
				} # next model
			
				modelOrder <- order(tuning$AICc)
				tuning <- tuning[modelOrder, , drop=FALSE]
				if ('models' %in% out) models <- models[modelOrder]
				rownames(tuning) <- NULL
			
			} # if selecing best model

		# NO AUTOMATED MODEL CONSTRUCTION
		# use all single-variable terms and two-variable terms
		} else if (!construct) {

			# single terms
			for (thisPred in preds) {

				if (class(data[ , thisPred]) != 'factor') {
					# form <- paste0(form, ' + s(', thisPred, ', bs=\'cs\', k=basisK)')
					form <- paste0(form, ' + s(', thisPred, ')')
				} else {
					form <- paste0(form, ' + ', thisPred)
				}

			}

			# interaction terms
			if (length(preds) > 1 & !is.null(interaction)) {

				numPreds <- length(preds)
			
				for (countPred1 in 1:(numPreds - 1)) { # for each initial predictor

					thisPred <- preds[countPred1]
				
					for (countPred2 in (countPred1 + 1):numPreds) {
					
						thatPred <- preds[countPred2]

						form <- if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) != 'factor') {
							paste0( form, ' + ', interaction, '(', thisPred, ', ', thatPred, ')')
						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) != 'factor') {
							paste0( form, ' + ', interaction, '(', thatPred, ', by=', thisPred, ')')
						} else if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) == 'factor') {
							paste0( form, ' + ', interaction, '(', thisPred, ', by=', thatPred, ')')
						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) == 'factor') {
							paste0(form, ' + ', thisPred, ' * ', thatPred)
						}

					}

				}
				
			}

			# train FULL model
			model <- mgcv::gam(
				formula=stats::as.formula(form),
				family=family,
				data=data,
				method='ML',
				optimizer=c('outer', 'newton'),
				select=TRUE,
				gamma=gamma,
				weights=w,
				na.action='na.fail',
				...
			)

		} # if NOT doing automated model construction

	### summary
	if (verbose & 'model' %in% out) {

		omnibus::say('Best model:', level=2);
		print(summary(model))
		utils::flush.console()

	}

	### return
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
