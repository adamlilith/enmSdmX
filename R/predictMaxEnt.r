#' Predict a MaxEnt model object (with optional feature-level permutation)
#'
#' Takes a MaxEnt \code{lambda} object or a MaxEnt object and returns raw or logistic predictions.  Its output is the same as the \code{\link[terra]{predict}} function from the \pkg{terra} package, and in fact, is slower than the function from \pkg{terra}. However, this function does allow custom manipulations that those functions do not allow (e.g., permuting product features while leaving other features with the same variables intact).  This function does \emph{not} clamp predictions--beyond the range of the training data, it extends the prediction in the direction it was going (up/down/no change). The function is based on Peter D. Wilson's document "Guidelines for computing MaxEnt model output values from a lambdas file". The function has a special feature in that it allows you to permute single variables or combinations of variables in specific features before making predictions. This is potentially useful, for example, if you wanted to determine the relative importance of a quadratic feature for a particular variable in a Maxent model relative to the other features in the model.  You can also permute values of a variable regardless of which features they appear in. For product features, you can implement the permutation before or after the values are multiplied together (before often makes for bigger differences in predictions).
#' @param x Either a Maxent lambda object or a Maxent model object
#' @param data Data frame. Data to which to make predictions
#' @param type Character.  One of:
#' \itemize{
#' 		\item \code{'raw'}: Maxent "raw" values
#' 		\item \code{'logistic'}: Maxent logistic values
#' 		\item \code{'cloglog'} complementary log-log output (as per version 3.4.0+ of maxent--called "\code{maxnet()}" in the package of the same name)
#' }
#' @param perm Character vector. Name(s) of variable to permute before calculating predictions. This permutes the variables for \emph{all} features in which they occur.  If a variable is named here, it overrides permutation settings for each feature featType.  Note that for product features the variable is permuted before the product is taken. This permutation is performed before any subsequent permutations (i.e., so if both variables in a product feature are included in \code{perms}, then this is equivalent to using the \code{'before'} rule for \code{permProdRule}). Ignored if \code{NULL}.
#' @param permLinear Character list. Names(s) of variables to permute in linear features before calculating predictions.  Ignored if \code{NULL}.
#' @param permQuad Names(s) of variables to permute in quadratic features before calculating predictions.  Ignored if \code{NULL}.
#' @param permHinge Character vector. Names(s) of variables to permute in forward/reverse hinge features before calculating predictions.  Ignored if \code{NULL}.
#' @param permThresh Character list. Names(s) of variables to permute in threshold features before calculating predictions.  Ignored if \code{NULL}.
#' @param permProd Character list. A list object of \code{n} elements, each of which has two character elements naming the variables to permute if they occur in a product feature.  Depending on the value of \code{permProdRule}, the function will either permute the individual variables then calculate their product or calculate their product, then permute the product across observations.  Any other features containing the variables will produce values as normal.  Example: \code{permProd=list(c('precipWinter', 'tempWinter'), c('tempSummer', 'precipFall'))}.  The order of the variables in each element of \code{permProd} doesn't matter, so \code{permProd=list(c('temp', 'precip'))} is the same as \code{permProd=list(c('precip', 'temp'))}.  Ignored if \code{NULL}.
#' @param permProdRule Character. Rule for how permutation of product features is applied: \code{'before'} ==> Permute individual variable values then calculate product; \code{'after'} ==> calculate product then permute across these values. Ignored if \code{permProd} is \code{NULL}.
#' @param ... Extra arguments (not used).
#' @return Numeric.
#' @seealso \code{\link[dismo]{maxent}}
#' @examples
#'
#' # The example below shows a very basic modeling workflow. It has been 
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
#' mxMap <- predictEnmSdm(mx, madEnv)
#'
#'
#' @export
predictMaxEnt <- function(
	x,
	data,
	type = 'cloglog',
	perm = NULL,
	permLinear = NULL,
	permQuad = NULL,
	permHinge = NULL,
	permThresh = NULL,
	permProd = NULL,
	permProdRule = NULL,
	...
) {

	# get lambdas is object is a model object
	if (inherits(x, 'MaxEnt')) x <- x@lambdas

	### create meta value data frame
	################################
	metaLines <- x[which(grepl(x, pattern='linearPredictorNormalizer')):length(x)]
	meta <- data.frame(
		item=substr(metaLines, 1, unlist(gregexpr(metaLines, pattern=',')) - 1),
		value=as.numeric(substr(metaLines, unlist(gregexpr(metaLines, pattern=',')) + 2, nchar(metaLines)))
	)

	### create data frame of feature functions
	##########################################
	feats <- data.frame(
		feature=x[1L:(which(grepl(x, pattern='linearPredictorNormalizer')) - 1L)],
		featType='linear',
		var1=NA,
		var2=NA,
		crit=NA,
		lambda=NA,
		beta1=NA,
		beta2=NA
	)

	feats$feature <- as.character(feats$feature)
	feats$featType <- as.character(feats$featType)
	feats$crit <- as.numeric(feats$crit)
	feats$lambda <- as.numeric(feats$lambda)
	feats$beta1 <- as.numeric(feats$beta1)
	feats$beta2 <- as.numeric(feats$beta2)

	## identify feature types
	feats$featType[which(grepl(feats$feature, pattern='[(]') & grepl(feats$feature, pattern='[)]'))] <- 'threshold'
	feats$featType[which(grepl(feats$feature, pattern='\\^'))] <- 'quadratic'
	feats$featType[which(grepl(feats$feature, pattern='[*]'))] <- 'product'
	feats$featType[which(grepl(feats$feature, pattern='\''))] <- 'forward hinge'
	feats$featType[which(grepl(feats$feature, pattern='`'))] <- 'reverse hinge'

	## identify variables used
	for (i in 1L:nrow(feats)) {
	
		vars <- strsplit(feats$feature[i], split=', ')[[1]][[1]]
		
		vars <- gsub(vars, pattern='\\^2', replacement='')
		vars <- gsub(vars, pattern='[(]', replacement='')
		vars <- gsub(vars, pattern='[)]', replacement='')
		vars <- gsub(vars, pattern='`', replacement='')
		vars <- gsub(vars, pattern='\'', replacement='')
		vars <- strsplit(vars, split='<')[[1]]
		if (length(vars) == 2) vars <- vars[2]
		vars <- strsplit(vars, split='\\*')[[1]]
		
		feats$var1[i] <- vars[1]
		feats$var2[i] <- if (length(vars) == 2) { vars[2] } else { NA }
	
	}

	## extract "betas"
	commaPos <- unlist(gregexpr(pattern=',', feats$feature))
	feats$beta1 <- as.numeric(substr(feats$feature, start=commaPos[seq(2, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(3, 3 * nrow(feats), by=3)] - 1))
	feats$beta2 <- as.numeric(substr(feats$feature, start=commaPos[seq(3, 3 * nrow(feats), by=3)] + 2, stop=nchar(feats$feature)))

	## extract critical values (hinge/threshold)
	feats$crit[which(feats$featType=='threshold')] <- as.numeric(substr(feats$feature[which(feats$featType=='threshold')], start=2, stop=unlist(gregexpr(feats$feature[which(feats$featType=='threshold')], pattern='<')) - 1))
		
	feats$crit[which(feats$featType=='forward hinge')] <- feats$beta1[which(feats$featType=='forward hinge')]
	feats$crit[which(feats$featType=='reverse hinge')] <- feats$beta2[which(feats$featType=='reverse hinge')]

		
	## extract lambdas
	feats$lambda <- as.numeric(substr(feats$feature, start=commaPos[seq(1, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(2, 3 * nrow(feats), by=3)] - 1))

	## get values
	# value1 <- if (nrow(data) > 1L) {
		# t(data[ , feats$var1])
	# } else {
		# cbind(data[ , feats$var1])
	# }
	value1 <- t(data[ , feats$var1])
	rownames(value1) <- paste0('feat', 1:nrow(feats))
	value2 <- value1 * NA
		
	if ('product' %in% feats$featType) {
		value2Compact <- t(data[ , feats$var2[which(feats$featType=='product')]])
		value2[which(feats$featType=='product'), ] <- value2Compact
	}

	### permutations
	################

	## permute values REGARDLESS OF FEATURE TYPE
	if (!is.null(perm)) {

		for (i in 1:nrow(feats)) {
			if (any(perm %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
			if (any(perm %in% feats$var2[i])) value2[i, ] <- sample(value2[i, ], ncol(value2))
		}
			
	}

	## permute values in LINEAR features
	if (!is.null(permLinear) & 'linear' %in% feats$featType) {
		for (i in which(feats$featType %in% 'linear')) {
			if (any(permLinear %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in QUADRATIC features
	if (!is.null(permQuad) & 'quadratic' %in% feats$featType) {
		for (i in which(feats$featType %in% 'quadratic')) {
			if (any(permQuad %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in HINGE features
	if (!is.null(permHinge) & ('forward hinge' %in% feats$featType | 'reverse hinge' %in% feats$featType)) {
		for (i in sort(unique(c(which(feats$featType %in% 'forward hinge'), which(feats$featType %in% 'reverse hinge'))))) {
			if (any(permHinge %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in THRESHOLD features
	if (!is.null(permThresh) & 'threshold' %in% feats$featType) {
		for (i in which(feats$featType %in% 'threshold')) {
			if (any(permThresh %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in PRODUCT features
	if (!is.null(permProd) & 'product' %in% feats$featType) {

		# permute variables then calculate product
		if (permProdRule=='before') {
		
			value1Perm <- value1
			value2Perm <- value2
		
			for (i in which(feats$featType %in% 'product')) {

				# see if permutation is wanted for this particular var1 and var2 of this product feature
				wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
				
				if (sum(wantPairPerm)==2) {
					value1Perm[i, ] <- sample(value1Perm[i, ], ncol(value1Perm))
					value2Perm[i, ] <- sample(value2Perm[i, ], ncol(value2Perm))
				}

			}

			# calculate product
			prodMat <- value1Perm * value2Perm
		
		# calculate product then permute product
		} else if (permProdRule=='after') {
		
			# calculate product
			prodMat <- value1 * value2

			for (i in which(feats$featType %in% 'product')) {

				# see if permutation is wanted for this particular var1 and var2 of this product feature
				wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
				if (sum(wantPairPerm)==2) prodMat[i, ] <- sample(prodMat[i, ], ncol(prodMat))
				
			}
			
		}
			
	} else {

		# calculate unpermuted product
		prodMat <- value1 * value2
		
	}

	### calculate predictions
	#########################

	# initialize exponent
	S <- numeric(nrow(data))

	# for each datum
	for (i in 1:nrow(data)) {

		# for each featType of feature
		for (featType in c('linear', 'quadratic', 'product', 'forward hinge', 'reverse hinge', 'threshold')) {

			# if feature featType occurs in lambdas, evaluate
			if (featType %in% feats$featType) {

				thisValue1 <- value1[which(feats$featType==featType), i]
				thisValue2 <- value2[which(feats$featType==featType), i]
				if ('product' %in% feats$featType) thisProd <- prodMat[which(feats$featType==featType), i]
			
				thisLambda <- feats$lambda[which(feats$featType==featType)]
				thisCrit <- feats$crit[which(feats$featType==featType)]
				
				thisBeta1 <- feats$beta1[which(feats$featType==featType)]
				thisBeta2 <- feats$beta2[which(feats$featType==featType)]
				
				thisRange <- thisBeta2 - thisBeta1
				
				if (featType=='linear') S[i] <- S[i] + sum(thisLambda * (thisValue1 - thisBeta1) / thisRange, na.rm=TRUE)
				if (featType=='quadratic') S[i] <- S[i] + sum(thisLambda * (thisValue1^2 - thisBeta1) / thisRange, na.rm=TRUE)
				if (featType=='product') S[i] <- S[i] + sum(thisLambda * (thisProd - thisBeta1) / thisRange, na.rm=TRUE)
				# if (featType=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=TRUE)
				
				
				if (featType=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=TRUE)
				if (featType=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=TRUE)
				# if (featType=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=TRUE)
				
				
				if (featType=='threshold') S[i] <- S[i] + sum(thisLambda[which(thisValue1 > thisCrit)], na.rm=TRUE)

			}
			
		}

	}

	raw <- exp(S - meta$value[match('linearPredictorNormalizer', meta$item)]) / meta$value[match('densityNormalizer', meta$item)]

	if (type=='raw') {
		out <- raw
	} else if (type=='logistic') {
		out <- (raw * exp(meta$value[match('entropy', meta$item)])) / (1 + raw * exp(meta$value[match('entropy', meta$item)]))
	} else if (type=='cloglog') {
		entropy <- x[length(x)]
		entropy <- as.numeric(substr(entropy, regexpr(',', entropy) + 2, nchar(entropy)))
		out <- 1 - exp(-exp(entropy) * raw)
	}

	out

}
