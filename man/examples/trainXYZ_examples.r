# The examples below show a very basic modeling workflow. They have been 
# designed to work fast, not produce accurate, defensible models.

\donttest{

library(sf)
library(terra)
set.seed(123)

### setup data
##############

# environmental rasters
rastFile <- system.file('extdata/madEnv.tif', package='enmSdmX')
madEnv <- rast(rastFile)
madEnv <- madEnv / 100 # values were rounded to nearest 100th then * by 100

crs <- sf::st_crs(madEnv)

# lemur occurrence data
data(lemurs)
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
occs <- sf::st_as_sf(occs, coords=c('longitude', 'latitude'), crs=crs)
occEnv <- extract(madEnv, occs, ID=FALSE)
occEnv <- occEnv[complete.cases(occEnv), ]
	
# create 10000 background sites (or as many as raster can support)
bgEnv <- terra::spatSample(madEnv, 20000)
bgEnv <- bgEnv[complete.cases(bgEnv), ]
bgEnv <- bgEnv[1:min(10000, nrow(bgEnv)), ]

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

### calibrate models
####################

# MaxEnt
mx <- trainMaxEnt(
	data = env,
	resp = 'presBg',
	preds = predictors,
	regMult = 1, # too few values for reliable model, but fast
	verbose = TRUE
)

# MaxNet
mn <- trainMaxNet(
	data = env,
	resp = 'presBg',
	preds = predictors,
	regMult = 1, # too few values for reliable model, but fast
	verbose = TRUE
)

# generalized linear model (GLM)
# Normally, we'd center and standardize variables before modeling.
gl <- trainGLM(
	data = env,
	resp = 'presBg',
	preds = predictors,
	verbose = TRUE
)

# generalized additive model (GAM)
ga <- trainGAM(
	data = env,
	resp = 'presBg',
	preds = predictors,
	verbose = TRUE
)

# natural splines
ns <- trainNS(
	data = env,
	resp = 'presBg',
	preds = predictors,
	df = 1:2, # too few values for reliable model(?)
	verbose = TRUE
)

# boosted regression trees
envSub <- env[1:1049, ] # subsetting data to run faster
brt <- trainBRT(
	data = envSub,
	resp = 'presBg',
	preds = predictors,
	learningRate = 0.001, # too few values for reliable model(?)
	treeComplexity = 2, # too few values for reliable model, but fast
	minTrees = 1200, # minimum trees for reliable model(?), but fast
	maxTrees = 1200, # too small for reliable model(?), but fast
	tryBy = 'treeComplexity',
	anyway = TRUE, # return models that did not converge
	verbose = TRUE
)

# random forests
rf <- trainRF(
	data = env,
	resp = 'presBg',
	preds = predictors,
	numTrees = c(100, 500), # using at least 500 recommended, but fast!
	verbose = TRUE
)

### make maps of models
#######################

mxMap <- predictEnmSdm(mx, madEnv)
mnMap <- predictEnmSdm(mn, madEnv)
glMap <- predictEnmSdm(gl, madEnv)
gaMap <- predictEnmSdm(ga, madEnv)
nsMap <- predictEnmSdm(ns, madEnv)
brtMap <- predictEnmSdm(brt, madEnv)
rfMap <- predictEnmSdm(rf, madEnv)

maps <- c(
	mxMap,
	mnMap,
	glMap,
	gaMap,
	nsMap,
	brtMap,
	rfMap
)

names(maps) <- c('MaxEnt', 'MaxNet', 'GLM', 'GAM', 'Natural Splines', 'BRTs', 'RFs')
fun <- function() plot(occs[1], col='black', add=TRUE)
plot(maps, fun = fun, nc = 4)

### compare model responses to BIO12 (mean annual precipitation)
################################################################

# make a data frame holding all other variables at mean across occurrences,
# varying only BIO12
occEnvMeans <- colMeans(occEnv, na.rm=TRUE)
occEnvMeans <- rbind(occEnvMeans)
occEnvMeans <- as.data.frame(occEnvMeans)
climFrame <- occEnvMeans[rep(1, 100), ]
rownames(climFrame) <- NULL

minBio12 <- min(env$bio12)
maxBio12 <- max(env$bio12)
climFrame$bio12 <- seq(minBio12, maxBio12, length.out=100)

predMx <- predictEnmSdm(mx, climFrame)
predMn <- predictEnmSdm(mn, climFrame)
predGl <- predictEnmSdm(gl, climFrame)
predGa <- predictEnmSdm(ga, climFrame)
predNat <- predictEnmSdm(ns, climFrame)
predBrt <- predictEnmSdm(brt, climFrame)
predRf <- predictEnmSdm(rf, climFrame)


plot(climFrame$bio12, predMx,
xlab='BIO12', ylab='Prediction', type='l', ylim=c(0, 1))

lines(climFrame$bio12, predMn, lty='solid', col='red')
lines(climFrame$bio12, predGl, lty='dotted', col='blue')
lines(climFrame$bio12, predGa, lty='dashed', col='green')
lines(climFrame$bio12, predNat, lty=4, col='purple')
lines(climFrame$bio12, predBrt, lty=5, col='orange')
lines(climFrame$bio12, predRf, lty=6, col='cyan')

legend(
   'topleft',
   inset = 0.01,
   legend = c(
	'MaxEnt',
	'MaxNet',
	'GLM',
	'GAM',
	'NS',
	'BRT',
	'RF'
   ),
   lty = c(1, 1:6),
   col = c(
	'black',
	'red',
	'blue',
	'green',
	'purple',
	'orange',
	'cyan'
   ),
   bg = 'white'
)

}