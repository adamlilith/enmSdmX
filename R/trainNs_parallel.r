#' Calibrate a natural splines model
#'
#' This function constructs a natural-spline model piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param df Integer > 0 \emph{or} vector of integers > 0. Sets flexibility of model fit. See documentation for \code{\link[splines]{ns}}.  If \code{construct} is \code{TRUE}, then univariate models for each term will be evaluated using each value in \code{df}. Note that \code{NULL} is also valid, but it can create problems when used with other functions in this package (and usually defaults to \code{df=3} anyway).
#' @param interaction If \code{TRUE} (default), include two-way interaction terms.
#' @param presPerTermFinal Minimum number of presence sites per term in initial starting model.
#' @param maxTerms Maximum number of terms to be used in any model, not including the intercept (default is 8). Used only if \code{construct} is \code{TRUE}.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by AICc.
#' }
#' @param cores Number of cores to use. Default is 1.
#' @param parallelType Either \code{'doParallel'} (default) or \code{'doSNOW'}. Issues with parallelization might be solved by trying the non-default option.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device. Default is \code{FALSE}.
#' @param ... Arguments to send to \code{gam()} or \code{dredge()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{gam}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gam} object and the data frame.
#' @seealso \code{\link[splines]{ns}}, \code{\link[mgcv]{gam}}, \code{\link{trainGam}}
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
trainNs_parallel <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	df = 1:4,
	interaction = TRUE,
	method = 'glm.fit',
	presPerTermFinal = 10,
	maxTerms = 8,
	w = TRUE,
	cores = 1,
	parallelType = 'doParallel',
	out = 'model',
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########

		# degrees of freedom
		if (is.null(df)) df <- 'NULL'

		# response and predictors
		if (inherits(resp, c('integer', 'numeric'))) resp <- names(data)[resp]
		if (inherits(preds, c('integer', 'numeric'))) preds <- names(data)[preds]

		w <- .calcWeights(w, data = data, resp = resp)
		
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
		for (thisPred in preds) {
		
			for (thisDf in df) {
			
				term <- if (!inherits(data[ , thisPred], 'factor')) {
					paste0('splines::ns(', thisPred, ', df=', thisDf, ')')
				} else {
					thisPred
				}
				
				terms <- c(terms, term)
			}
		
		}

		# interaction terms
		if (interaction & length(preds) > 1L & n >= 2 * presPerTermFinal) {
		
			for (countPred1 in 1L:(length(preds) - 1L)) { # for each predictor test two-variable terms

				pred1 <- preds[countPred1]

				for (countPred2 in 2L:length(preds)) { # for each second predictor test two-variable terms

					pred2 <- preds[countPred2]
					terms <- c(terms, paste0('splines::ns(', preds[countPred1], ' * ', preds[countPred2], ', df=', df, ')'))
					
				} # next second term
				
			} # next first term
			
		} # if more than one term
		
		terms <- sort(terms)

	### model selection
	###################

		### make table of all possible terms
		
		# single-predictor terms
		candidates <- list(x = paste0('splines::ns(', preds[1L], ', df=', df, ')'))
		candidates$x <- c(candidates$x, NA)
		
		if (length(preds) > 1L) {
			for (i in 2L:length(preds)) {
				candidates <- c(
					candidates,
					list(x = c(paste0('splines::ns(', preds[i], ', df=', df, ')'), NA))
				)
			}
		
		}
		names(candidates) <- paste0('pred', seq_along(preds))
		
		# add interaction terms
		if (interaction & length(preds) > 1L & n >= 2 * presPerTermFinal) {
		
			for (countPred1 in 1L:(length(preds) - 1L)) { # for each predictor test two-variable terms

				pred1 <- preds[countPred1]

				for (countPred2 in 2L:length(preds)) { # for each second predictor test two-variable terms

					pred2 <- preds[countPred2]
					candidates$x <- c(paste0('splines::ns(', preds[countPred1], ' * ', preds[countPred2], ', df=', df, ')'), NA)
					names(candidates)[length(candidates)] <- paste0('pred', length(candidates))
					
				} # next second term
				
			} # next first term
			
		} # if more than one term


		candidates <- expand.grid(candidates, stringsAsFactors = FALSE)
		
		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[numTerms > 0, , drop=FALSE]

		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[numTerms <= maxTerms, , drop=FALSE]

		numTerms <- rowSums(!is.na(candidates))
		candidates <- candidates[n >= numTerms * presPerTermFinal, , drop=FALSE]
		
		forms <- rep(NA, nrow(candidates))
		for (i in 1L:nrow(candidates)) {
		
			terms <- which(!is.na(unlist(candidates[i , ])))
			terms <- unlist(candidates[i, terms])
			if (length(terms) > 1L) terms <- paste(terms, collapse = ' + ')
			forms[i] <- terms
		
		}
		
		forms <- paste0('1 + ', forms)
		if (interceptOnly) forms <- c(forms, '1')
		
		work <- foreach::foreach(
			i = seq_along(forms),
			.options.multicore = mcOptions,
			.combine = 'c',
			.inorder = FALSE,
			.export = c('.trainNsWorker')
		) %makeWork% {
			.trainNsWorker(
				i = i,
				forms = forms,
				data = data,
				resp = resp,
				family = family,
				method = method,
				w = w,
				insertIntercept = FALSE,
				paths = paths,
				modelOut = ('models' %in% out | 'model' %in% out),
				...
			)
		}
	
		if (cores > 1L) parallel::stopCluster(cl)

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
	
		aiccOrder <- order(tuning$AICc)
		if ('model' %in% out) model <- work[[aiccOrder[1L]]]$model
		if ('models' %in% out) {
			models <- list()
			models[[1]] <- work[[1L]]$model
			for (i in 2L:length(work)) models[[i]] <- work[[i]]$model
			models <- models[aiccOrder]
		}
		tuning <- tuning[aiccOrder, , drop = FALSE]
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

#################
### train GAM ###
#################

.trainNsWorker <- function(
	i,
	forms, # formulae (without LHS and maybe without intercept)
	data,
	resp,
	family,
	method,
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

	model <- stats::glm(
		formula = stats::as.formula(thisForm),
		data = data,
		family = family,
		method = method,
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
