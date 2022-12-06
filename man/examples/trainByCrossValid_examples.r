# The example below show a very basic modeling workflow. They have been 
# designed to work fast, not produce accurate, defensible models.
# The general idea is to calibrate a series of models and evaluate them
# against a withheld set of data. One can then use the series of models
# of the top models to better select a "final" model.

# Runing the entire set of commands can take quite a bit of time.

if (FALSE) {

library(sf)
library(terra)
set.seed(123)

### setup data
##############

# environmental rasters
rastFile <- system.file('extdata/madEnv.tif', package='enmSdmX')
madEnv <- rast(rastFile)
madEnv <- madEnv

crs <- sf::st_crs(madEnv)

# lemur occurrence data
data(lemurs)
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
occs <- sf::st_as_sf(occs, coords=c('longitude', 'latitude'), crs=crs)
occEnv <- extract(madEnv, occs, ID=FALSE)
occEnv <- occEnv[complete.cases(occEnv), ]
	
# create background sites (using just 1000 to speed things up!)
bgEnv <- terra::spatSample(madEnv, 3000)
bgEnv <- bgEnv[complete.cases(bgEnv), ]
bgEnv <- bgEnv[sample(nrow(bgEnv), 1000), ]

# collate occurrences and background sites
presBg <- data.frame(
   presBg = c(
      rep(1, nrow(occEnv)),
      rep(0, nrow(bgEnv))
   )
)

env <- rbind(occEnv, bgEnv)
env <- cbind(presBg, env)

predictors <- c('bio1', 'bio12')

# using "vector" form of "folds" argument
folds <- dismo::kfold(env, 3) # just 3 folds (for speed)

### calibrate models
####################

cores <- 1 # increase this to go faster, if your computer can handle it
parallelType <- 'doParallel' # if this doesn't work, try 'doSNOW'

## MaxEnt
mxx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainMaxEnt,
	regMult = 1:2, # too few values for valid model, but fast!
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize MaxEnt feature sets and regularization across folds
summaryByCrossValid(mxx)

## MaxNet
mnx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainMaxNet,
	regMult = 1:2, # too few values for valid model, but fast!
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize MaxEnt feature sets and regularization across folds
summaryByCrossValid(mnx)

## generalized linear models
glx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainGLM,
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize GLM terms in best models
summaryByCrossValid(glx)

## generalized additive models
gax <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainGAM,
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize GAM terms in best models
summaryByCrossValid(gax)

## natural splines
nsx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainNS,
	df = 1:2,
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize NS terms in best models
summaryByCrossValid(nsx)

## boosted regression trees
brtx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainBRT,
	learningRate = c(0.001, 0.0001), # too few values for reliable model(?)
	treeComplexity = c(2, 4), # too few values for reliable model, but fast
	minTrees = 1000,
	maxTrees = 1500, # too small for reliable model(?), but fast
	tryBy = 'treeComplexity',
	anyway = TRUE, # return models that did not converge
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize BRT parameters in best models
summaryByCrossValid(brtx)

## random forests
rfx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainRF,
	verbose = 1,
	cores = cores,
	parallelType = parallelType
)

# summarize RF parameters in best models
summaryByCrossValid(rfx)

}
