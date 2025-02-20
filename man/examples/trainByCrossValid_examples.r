# The example below show a very basic modeling workflow. It has been 
# designed to work fast, not produce accurate, defensible models.
# The general idea is to calibrate a series of models and evaluate them
# against a withheld set of data. One can then use the series of models
# of the top models to better select a "final" model.

\dontrun{
# Running the entire set of commands can take a few minutes. This can
# be sped up by increasing the number of cores used. The examples below use
# one core, but you can change that argument according to your machine's
# capabilities.

library(sf)
library(terra)
set.seed(123)

### setup data
##############

# environmental rasters
rastFile <- system.file('extdata/madClim.tif', package='enmSdmX')
madClim <- rast(rastFile)

# coordinate reference system
wgs84 <- getCRS('WGS84')

# lemur occurrence data
data(lemurs)
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
occs <- vect(occs, geom=c('longitude', 'latitude'), crs=wgs84)

occs <- elimCellDuplicates(occs, madClim)

occEnv <- extract(madClim, occs, ID = FALSE)
occEnv <- occEnv[complete.cases(occEnv), ]
	
# create background sites (using just 1000 to speed things up!)
bgEnv <- terra::spatSample(madClim, 3000)
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
folds <- predicts::kfold(env, 3) # just 3 folds (for speed)

### calibrate models
####################

cores <- 1 # increase this to go faster, if your computer handles it

## MaxEnt
mxx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainMaxEnt,
	regMult = 1:2, # too few values for valid model, but fast!
	verbose = 1,
	cores = cores
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
	cores = cores
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
	cores = cores
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
	cores = cores
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
	cores = cores
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
	cores = cores
)

# summarize BRT parameters across best models
summaryByCrossValid(brtx)

## random forests
rfx <- trainByCrossValid(
	data = env,
	resp = 'presBg',
	preds = c('bio1', 'bio12'),
	folds = folds,
	trainFx = trainRF,
	verbose = 1,
	cores = cores
)

# summarize RF parameters in best models
summaryByCrossValid(rfx)

}
