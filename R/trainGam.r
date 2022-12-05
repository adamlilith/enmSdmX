#' Calibrate a generalized additive model (GAM)
#'
#' This function constructs a GAM piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#'
#' @param data Data frame.
#' @param resp Response variable. This is either the name of the column in \code{data} or an integer indicating the column in \code{data} that has the response varoable. The default is to use the first column in \code{data} as the response.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. The default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}).
#' @param gamma Initial penalty to degrees of freedom to use (larger ==> smoother fits).
#' @param scale A numeric value indicating the "scale" parameter (see argument \code{scale} in \code{\link[mgcv]{gam}}). The default is 0 (which allows a single smoother for Poisson and binomial error families and unknown scale for all others.)
#' @param construct If \code{TRUE} (default), then construct the model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{maxTerms}.
#' @param select If \code{TRUE} (default), then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if \code{construct} is \code{TRUE}).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only if \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is \code{TRUE}.
#' @param maxTerms Maximum number of terms to be used in any model, not including the intercept (default is 8). Used only if \code{construct} is \code{TRUE}.
#' @param interaction Character or \code{NULL}. Type of interaction term to use (\code{te}, \code{ts}, \code{s}, etc.). See \code{?te} (for example) for help on any one of these. If \code{NULL} then interactions are not used.
#' @param interceptOnly If \code{TRUE} (default) and model selection is enabled, then include an intercept-only model.
#' @param smoothingBasis Character. Indicates the type of smoothing basis. The default is \code{'cs'} (cubic splines), but see \code{\link[mgcv]{smooth.terms}} for other options. This is the value of argument \code{bs} in a \code{\link[mgcv]{s}} function.
#' @param w Weights. Any of:
#' \itemize{
#'	\item \code{TRUE}: Causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'})
#' 	\item \code{FALSE}: Each datum is assigned a weight of 1.
#'  \item A numeric vector of weights, one per row in \code{data}.
#' 	\item The name of the column in \code{data} that contains site weights.
#' }
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by AICc.
#' }
#' @param cores Integer >= 1. Number of cores to use when calculating multiple models. Default is 1.
#' @param parallelType Either \code{'doParallel'} (default) or \code{'doSNOW'}. Issues with parallelization might be solved by trying the non-default option.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device.
#' @param ... Extra arguments (not used).
#'
#' @return The object that is returned depends on the value of the \code{out} argument. It can be a model object, a data frame, a list of models, or a list of all two or more of these.
#'
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
	scale = 0,
	smoothingBasis = 'cs',
	interaction = 'te',
	interceptOnly = TRUE,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 10,
	maxTerms = 8,
	w = TRUE,
	out = 'model',
	cores = 1,
	parallelType = 'doParallel',
	verbose = FALSE,
	...
) {

	### setup
	#########

		# ellipses <- list(...)

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		w <- .calcWeights(w, data = data, resp = resp, family = family)

	### parallelization
	###################
			
		cores <- if (!construct) {
			1L
		} else {
			min(cores, parallel::detectCores(logical = FALSE))
		}

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

	### make list of candidate model terms
	######################################

		n <- if (family %in% c('binomial', 'quasibinomial')) {
			sum(data[ , resp, drop=TRUE])
		} else {
			nrow(data)
		}

		### create vector of terms
		terms <- character()
		
		# single-predictor terms
		for (thisPred in preds) {
		
			if (!inherits(data[ , thisPred], 'factor')) {
				terms <- c(terms, paste0('s(', thisPred, ', bs=\'', smoothingBasis, '\')'))
			} else {
				thisPred
			}
		}
		
		# interaction terms
		if (length(preds) > 1L & n >= 2 * presPerTermInitial) {
		
			for (countPred1 in 1L:(length(preds)-1L)) { # for each predictor test two-variable terms

				pred1 <- preds[countPred1]

				for (countPred2 in 2:length(preds)) { # for each second predictor test two-variable terms

					pred2 <- preds[countPred2]

					# create term
					if (!inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

						terms <- c(terms, paste0(interaction, '(', pred1, ', ', pred2, ', bs=\'', smoothingBasis, '\')'))

					} else if (inherits(data[ , pred1], 'factor') & !inherits(data[ , pred2], 'factor')) {

						terms <- c(terms, paste0(interaction, '(', pred2, ', by=', pred1, ', bs=\'', smoothingBasis,'\')'))

					} else if (!inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

						terms <- c(terms, paste0(interaction, '(', pred1, ', by=', pred2, ', bs=\'', smoothingBasis, '\')'))

					} else if (inherits(data[ , pred1], 'factor') & inherits(data[ , pred2], 'factor')) {

						terms <- c(terms, paste0(pred1, ' * ', pred2))

					}
					
				} # next second term
				
			} # next first term
			
		} # if more than one term
			
	## term-by-term model construction
	##################################
	if (construct) {

		assess <- foreach::foreach(
			i = seq_along(terms),
			.options.multicore = mcOptions,
			.combine = 'rbind',
			.inorder = FALSE,
			.export = c('.trainGamWorker')
		) %makeWork% {
			.trainGamWorker(
				i = i,
				forms = terms,
				data = data,
				resp = resp,
				family = family,
				gamma = gamma,
				scale = scale,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = FALSE,
				...
			)
				
		}
		
		assess <- assess[order(assess$AICc), , drop=FALSE]
		rownames(assess) <- NULL

		if (verbose) {
			omnibus::say('Term-by-term evaluation:', level=2)
			print(assess)
			utils::flush.console()
		}

		### no selection
		################
		
		# just return model with "best" terms without further selection
		
		if (!select) {
			
			requiredPresPerTerm <- presPerTermFinal * 1L:nrow(assess)
			
			these <- which(requiredPresPerTerm <= n)
			theseTerms <- assess$formula[these]
			form <- paste(theseTerms, collapse=' + ')
			thisForm <- paste0(resp, ' ~ 1 + ', form)
			
			model <- mgcv::gam(
				formula = stats::as.formula(thisForm),
				family = family,
				data = data,
				method = 'ML',
				optimizer = c('outer', 'newton'),
				scale = scale,
				select = TRUE,
				gamma = gamma,
				weights = w,
				...
			)
			
			AICc <- MuMIn::AICc(model)
			
			tuning <- data.frame(
				model = form,
				AICc = AICc
			)

			models <- NULL

			if (verbose) {

				omnibus::say('Final model (construction from best terms, but no selection):', level=2)
				print(summary(model))
				utils::flush.console()
				
			}

		### model selection
		###################
		
		# select best model from all possible subsets of "full" model with best terms
		
		} else {

			# make formulae for all possible models
			candidates <- list(x = c(TRUE, FALSE))
			if (length(terms) > 1L) for (i in 2L:length(terms)) candidates <- c(candidates, list(x = c(TRUE, FALSE)))
			candidates <- expand.grid(candidates, stringsAsFactors = FALSE)
			names(candidates) <- terms
			candidates <- candidates[rowSums(candidates) != 0, , drop=FALSE]

			numTerms <- rowSums(candidates)
			requiredPresPerTerm <- numTerms * presPerTermFinal
			candidates <- candidates[requiredPresPerTerm <= n, , drop = FALSE]
			
			numTerms <- rowSums(candidates)
			candidates <- candidates[numTerms <= maxTerms, , drop=FALSE]

			forms <- character()
			for (i in 1L:nrow(candidates)) {
				thisForm <- paste0('1 + ', paste(terms[unlist(candidates[i, , drop = TRUE])], collapse = ' + '))
				forms <- c(forms, thisForm)
			}
			if (interceptOnly) forms <- c(forms, '1')

			work <- foreach::foreach(
				i = seq_along(forms),
				.options.multicore = mcOptions,
				.combine = 'c',
				.inorder = FALSE,
				.export = c('.trainGamWorker')
			) %makeWork% {
				.trainGamWorker(
					i = i,
					forms = forms,
					data = data,
					resp = resp,
					family = family,
					gamma = gamma,
					scale = scale,
					w = w,
					insertIntercept = FALSE,
					paths = paths,
					modelOut = ('models' %in% out | 'model' %in% out),
					...
				)
			}
			
			# tuning table
			tuning <- data.frame(
				model = work[[1L]]$formula,
				AICc = work[[1L]]$AICc
			)
			
			if (length(work) > 1L) {
				for (i in 2L:length(work)) {
					
					tuning <- rbind(
						tuning,
						data.frame(
							model = work[[i]]$formula,
							AICc = work[[i]]$AICc
						)
					)
					
				}
			}
		
			bestOrder <- order(tuning$AICc)
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
		
		} # if selecting best model from subsets of "full" model

	### if not constructing model term-by-term (selection not possible)
	###################################################################
	} else {
	
		form <- paste(terms, collapse = ' + ')
		thisForm <- paste0(resp, ' ~ 1 + ', form)
	
		model <- mgcv::gam(
			formula = stats::as.formula(thisForm),
			family = family,
			data = data,
			method = 'ML',
			optimizer = c('outer', 'newton'),
			scale = scale,
			select = TRUE,
			gamma = gamma,
			weights = w,
			...
		)
		
		AICc <- MuMIn::AICc(model)
		
		tuning <- data.frame(
			model = form,
			AICc = AICc
		)
		
		models <- NULL
		
		if (select) warning('Model selection is not performed when argument "construct" is FALSE.')
		
		if (verbose) {
		
			omnibus::say('Model (no construction or selection):', level=2)
			print(summary(model))
			utils::flush.console()
		
		}
		
	} # if not constructing model term-by-term

	if (cores > 1L) parallel::stopCluster(cl)

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

#################
### train GAM ###
#################

.trainGamWorker <- function(
	i,
	forms, # formulae (without LHS and maybe without intercept)
	data,
	resp,
	family,
	gamma,
	scale,
	w,
	paths,
	insertIntercept, # if TRUE, add "1 +" to the RHS side
	modelOut, # if TRUE, return model *and* data frame with model
	...
) {

	 # need to call this to avoid "object '.doSnowGlobals' not found" error!!!
	.libPaths(paths)

	form <- forms[i]
	if (insertIntercept) {
		form <- if (form == '') {
			'1'
		} else {
			paste('1', form, sep=' + ')
		}
	}
	thisForm <- paste0(resp, ' ~ ', form)

	model <- mgcv::gam(
		formula = stats::as.formula(thisForm),
		family = family,
		data = data,
		method = 'ML',
		optimizer = c('outer', 'newton'),
		scale = scale,
		select = TRUE,
		gamma = gamma,
		weights = w,
		...
	)
	
	AICc <- MuMIn::AICc(model)
	
	# out
	out <- if (modelOut) {
		
		list(
			list(
				model = model,
				formula = form,
				AICc = AICc
			)
		)
		
	} else {
	
		data.frame(
			formula = form,
			AICc = AICc
		)
	
	}
		
	out

}
