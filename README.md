# enmSdmX
Tools for modeling niches and distributions of species

<img align="right" src="enmSdmX.png" height="223"/>

This package for R contains a suite of efficiency functions for species distribution modeling and ecological niche modeling. You can install this package from CRAN using:

`install.packages('enmSdmX', dependencies = TRUE)`

You can install the development version of this package using:

`remotes::install_github('adamlilith/enmSdmX', dependencies=TRUE)`  

You may need to install the `remotes` package first.

# Functions #

Most functions are prefixed with the main data type they create or manipulate:
* `crs*`: Create or lookup coordinate reference systems
* `eval*`: Evaluate a SDM/ENM or niche overlap
* `rast*`: Create or manipulate rasters
* `points*`: Create or manipulate spatial points
* `predict*`: Predict from a model object
* `train*`: Train a SDM/ENM
}


### Data preparation ###
* `pointGeoFold`: Assign geographically-distinct k-folds
* `pointElimCellDups`: Eliminate duplicate points in each cell of a raster

### Using spatially imprecise records
* `nearestGeogPoints`: Minimum convex polygon from a set of spatial polygons and/or points ("nearest geographic point" method)
* `nearestEnvPoints`:  Extract "most conservative" environments from points and/or polygons ("nearest environmental point" method)
* `pointImprecision`: Calculate maximum possible coordinate precision

### Bias correction
* `pointDistWeights`: Proximity-based weighting for occurrences to correct for spatial bias
* `pointGeoThin`: Thin geographic points deterministically or randomly

### Model training ###
* `trainByCrossValid` and `summaryByCrossValid`: Calibrate a distribution/niche model using cross-validation
* `trainBRT`: Boosted regression trees (BRTs)
* `trainGAM`: Generalized additive models (GAMs)
* `trainGLM`: Generalized linear models (GLMs)
* `trainMaxEnt`: MaxEnt models
* `trainMaxNet`: MaxNet models
* `trainNS`: Natural splines (NSs)
* `trainRF`: Random forests (RFs)  

### Model prediction ###
* `predictEnmSdm` Predict most model types using default settings; parallelized
* `predictMaxEnt` Predict MaxEnt model
* `predictMaxNet` Predict MaxNet model

### Model evaluation ###
* `evalAUC`: AUC (with/out site weights)
* `evalMultiAUC`: Multivariate version of AUC (with/out site weight)
* `evalContBoyce`: Continuous Boyce Index (with/out site weights)
* `evalThreshold`: Thresholds to convert continuous predictions to binary predictions (with/out site weights)
* `evalThresholdStats`: Model performance statistics based on thresholded predictions (with/out site weights)
* `evalTjursR2`: Tjur's R2 (with/out site weights)
* `evalTSS`: True Skill Statistic (TSS) (with/out site weights)
* `modelSize`: Number of response values in a model object

### Niche overlap ###
* `evalNicheOverlap`: Niche overlap metrics
* `compareResponse`: Compare niche model responses to a single variable

### Functions for rasters ###
* `rastGetValueByCell`: Retrieve raster values(s) by cell number
* `rastInterpolate`: Interpolate a stack of rasters
* `rastLongLat`: Generate rasters with values of longitude/latitude for cell values
* `rastSquareCells`: Create a raster with square cells from an object with an extent
* `rastVelocity`: Velocity of "movement" of mass across a series of rasters
* `rastSample` : Sample raster with/out replacement
* `rastSetValueByCell`: Set raster values(s) by cell number
* `squareRastCells`: Resample a raster so cells are square

### Coordinate reference systems ###
* `crss`: Coordinate reference systems
* `crsGet`: Return a WKT2 (well-known text) string using a nickname
* `crsLambert`: Create a custom Lambert azimuthal equal-area projection
* `crsVertical`: Create a custom "vertical near-side" projection

### Geographic utility functions ###
* `crsVertical`: Generate a "Vertical Near-Side" projection WKT2 string
* `plotExtent`: Create a spatial polygon the same size as a plot region
* `decimalToDms`: Convert decimal coordinate to degrees-minutes-seconds
* `dmsToDecimal`: Convert degrees-minutes-seconds coordinate to decimal
* `extentToVect`: Convert extent to a spatial polygon
* `spatVectorToSpatial`: Convert SpatVector object to a Spatial* object

### Data
* `lemurs`: Lemur occurrences
* `mad0`: Madagascar spatial object
* `mad1`: Madagascar spatial object
* `madEnv`: Madagascar climate rasters

# Citation #
As of December 2022 there is no package-specific publication for `enmSdmX`, but the package was first used and cited in:

Smith, A.B., Murphy, S.J., Henderson, D., and Erickson, K.D. 2023. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth.  <italic>Global Ecology and Biogeography</italic> In press. <a href='http://dx.doi.org/10.1101/2021.06.10.447988'>Open access pre-print</a>

<i>Abstract</i>
<i>Aim</i> Museum and herbarium specimen records are frequently used to assess species’ conservation status and responses to climate change. Typically, occurrences with imprecise geolocality information are discarded because they cannot be matched confidently to environmental conditions, and are thus expected to increase uncertainty in downstream analyses. However, using only precisely georeferenced records risks undersampling of species’ environmental and geographic distributions. We present two related methods to allow the use of imprecisely georeferenced occurrences in biogeographic analysis.

<i>Innovation</i> Our two procedures assign imprecise records to the 1) locations or 2) climates that are closest to the geographic or environmental centroid of the precise records of a species. For virtual species, including imprecise records alongside precise records improved the accuracy of ecological niche models projected to the present and the future, especially for species with ~20 or fewer precise occurrences. Using only precise records underestimates loss in suitable habitat and overestimates the amount of suitable habitat in both the present and future. Including imprecise records also improves estimates of niche breadth and extent of occurrence. An analysis of 44 species of North American <i>Asclepias</i> (Apocynaceae) yielded similar results.

<i>Main conclusions</i> Existing studies examining the effects of spatial imprecision compare outcomes based on precise records to the same records with spatial error added to them. However, in real-world cases, analysts possess a mix of precise and imprecise records and must decide whether to retain or discard the latter. Discarding imprecise records can undersample species’ geographic and environmental distributions and lead to mis-estimation of responses to past and future climate change. Our method, for which we provide a software implementation in the enmSdmX package for R, is simple to employ and can help leverage the large number of specimen records that are typically deemed “unusable” because of spatial imprecision in their geolocation.
