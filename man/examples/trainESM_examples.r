# NB: The examples below show a very basic modeling workflow. They have been 
# designed to work fast, not produce accurate, defensible models. They can
# take a few minutes to run.

library(mgcv)
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
	
# create 10000 background sites (or as many as raster can support)
bgEnv <- terra::spatSample(madClim, 20000)
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

# "traditional" ESMs with just 2 linear predictors
# just one model in this case because we have just 2 predictors
esm1 <- trainESM(
   data = env,
   resp = 'presBg',
   preds = predictors,
   family = stats::binomial(),
   scale = TRUE,
   w = TRUE
)

str(esm1, 1)
esm1$tuning

# extended ESM with other kinds of terms
esm2 <- trainESM(
   data = env,
   resp = 'presBg',
   preds = predictors,
   univariate = TRUE,
   quadratic = TRUE,
   interaction = TRUE,
   interceptOnly = TRUE,
   family = stats::binomial(),
   scale = TRUE,
   w = TRUE,
   verbose = TRUE
)

str(esm2, 1)
esm2$tuning
