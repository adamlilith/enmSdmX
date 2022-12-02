#' Calibrate a natural splines model
#'
#' This function constructs a natural-spline model piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param df Integer > 0 \emph{or} vector of integers > 0. Sets flexibility of model fit. See documentation for \code{\link[splines]{ns}}.  If \code{construct} is \code{TRUE}, then univariate models for each term will be evaluated using each value in \code{df}. Note that \code{NULL} is also valid, but it can create problems when used with other functions in this package (and usually defaults to \code{df=3} anyway).
#' @param construct Logical. If TRUE then construct model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{initialTerms}.
#' @param select Logical. If TRUE then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is \code{TRUE}.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{\link[MuMIn]{dredge}} is 31, so if this number is >31 and \code{select} is \code{TRUE} then it is forced to 31 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character vector. One or more values:
#' \itemize{
#' 	\item	\code{'model'}: Model with the lowest AICc.
#' 	\item	\code{'models'}: All models evaluated, sorted from lowest to highest AICc (lowest is best).
#' 	\item	\code{'tuning'}: Data frame with tuning patrameters, one row per model, sorted by AICc.
#' }
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
trainNs <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	df = 1:4,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 10,
	initialTerms = 8,
	w = TRUE,
	out = 'model',
	verbose = FALSE,
	...
) {

	ellipses <- list(...)

	# force number of starting terms to 31 or less
	if (select & initialTerms > 31) {
		initialTerms <- 31
		warning('initialTerms must be 31 or less. Forcing to 31.')
	}

	# degrees of freedom
	if (is.null(df)) df <- 'NULL'

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	#############
	## weights ##
	#############

	# model weights
	if (is.logical(w)) {
		if (w[1L] & (family %in% c('binomial', 'quasibinomial'))) {
			posCases <- sum(data[ , resp, drop=TRUE] == 1)
			negCases <- sum(data[ , resp, drop=TRUE] == 0)
			w <- c(rep(1, posCases), rep(posCases / (posCases + negCases), negCases))
		} else {
			w <- rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w, drop=TRUE]
	}
	w <<- w / max(w) # declare to global because dredge() and pdredge() have problems if it is not

	################################
	## initial model construction ##
	################################

	# create starting formula
	form <- paste0(resp, ' ~ 1')

	if (construct) {

		### SINGLE-variable terms
		for (thisPred in preds) { # for each predictor test single-variable terms

			for (thisDf in df) {
			
				term <- if (class(data[ , thisPred]) != 'factor') {
					paste0('splines::ns(', thisPred, ', df=', thisDf, ')')
				} else {
					thisPred
				}

				thisThisForm <- paste0(form, ' + ', term)

				thisModel <- stats::glm(stats::as.formula(thisThisForm), family=family, data=data, weights=w, ...)
				thisAicc <- MuMIn::AICc(thisModel)

				# remember
				tuning <- if (exists('tuning', inherits=FALSE)) {
					rbind(tuning, data.frame(term=term, AICc=thisAicc, df=thisDf))
				} else {
					data.frame(term=term, AICc=thisAicc, df=thisDf)
				}
				
			}

		} # next single-variable term

		# sort by AIC
		tuning <- tuning[order(tuning$AICc), ]

		# print AICc frame
		if (verbose) {

			omnibus::say('Model construction results for each term tested:', level=2)
			print(tuning)
			omnibus::say('')

		}

		## construct final model
		form <- paste0(form, ' + ', tuning$term[1]) # add first term

		# for each set of presences > min num required, add a term
		if (floor(sum(data[ , resp]) / presPerTermInitial ) - 1 > 1 & initialTerms > 1) {

			termsToAdd <- 2:max(2, min(initialTerms, c(floor(sum(data[ , resp]) / presPerTermInitial ) - 1, nrow(tuning) - 1)))

			form <- paste0(form, ' + ', paste0(tuning$term[termsToAdd], collapse=' + '))

		} # if there are sufficient presences for additional terms beyond first

	# NO AUTOMATED MODEL CONSTRUCTION
	# use all single-variable terms and two-variable terms
	} else {

		if (length(df) > 1) warning('Multiple values of "df" assigned. Using the first one.')
	
		# single terms
		for (thisPred in preds) {

			if (class(data[ , thisPred]) != 'factor') {
				form <- paste0(form, ' + splines::ns(', thisPred, ', df=', df[1], ')')
			} else {
				form <- paste0(form, ' + ', thisPred)
			}

		}

	} # if not doing automated model construction

	###########################################################################
	## train model ############################################################
	## while model hasn't converged and while gamma is <= tolerance value... ##
	###########################################################################

	# get GLM model... using automated scale selection with weights so influence of absences equals influence of presences... using tryCatch because sometimes for variables with too little variation the default df of the basis is too high, in which case it is reduced and attempted again (for univariate and bivariate models only)
	model <- stats::glm(stats::as.formula(form), family=family, data=data, weights=w, na.action='na.fail', ...)

	if (verbose) {

		omnibus::say('Starting full model:', level=2)
		print(summary(model))
		omnibus::say('')

	}

	#######################################################################################
	## if doing model construction, evaluate all possible models using AIC then get best ##
	#######################################################################################

	if (select) {

		if (verbose) omnibus::say('Calculating AICc across all possible models...')

		# calculate all possible models and rank by AIC
		lims <- c(0, max(1, min(c(floor(sum(data[ , resp]) / presPerTermFinal), initialTerms, nrow(tuning)))))

		tuningModels <- MuMIn::dredge(
			global.model=model,
			rank='AICc',
			m.lim=lims,
			trace=FALSE,
			...
		)

		# get model with best AIC
		model <- MuMIn::get.models(tuningModels, subset = 1)[[1]]
		if ('models' %in% out) models <- MuMIn::get.models(tuningModels, subset=TRUE)
		
		if (verbose) {
		
			omnibus::say('Final model:', level=2)
			print(summary(model))
			omnibus::say('')
		
		}

		rownames(tuningModels) <- NULL
	
	} # if model selection

	# return
	if (length(out) > 1) {
		output <- list()
		if ('models' %in% out) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuningModels
		output
	} else if (out == 'models') {
		models
	} else if (out == 'model') {
		model
	} else if (out == 'tuning') {
		tuningModels
	}
	
}
