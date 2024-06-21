#' enmSdmX: Species distribution modeling and ecological niche modeling
#'
#' Tools for implementing species distribution models and ecological niche models, including: bias correction, spatial cross-validation, model evaluation, raster interpolation, biotic "velocity" (speed and direction of movement of a "mass" represented by a raster), and tools for using spatially imprecise records. The heart of the package is a set of "training" functions which automatically optimize model complexity based number of available occurrences. These algorithms include MaxEnt, MaxNet, boosted regression trees/gradient boosting machines, generalized additive models, generalized linear models,	natural splines, and random forests. To enhance interoperability to and from other packages, the package does not create any new classes. The package works with PROJ6 geodetic objects and coordinate reference systems.\cr \cr
#'
#' Create an issue on \href{https://github.com/adamlilith/enmSdmX/issues}{GitHub}.
#'
#' @details
#' @section Using imprecisely-georeferenced occurrences:
#' 		\code{\link{coordImprecision}}: Coordinate imprecision \cr
#' 		\code{\link{nearestEnvPoints}}: Extract "most conservative" environments from points and/or polygons \cr
#' 		\code{\link{nearestGeogPoints}}: Create a minimum convex polygon from a set of spatial polygons and/or points \cr
#'
#' @section Data preparation:
#' 		\code{\link{geoFold}}: Assign geographically-distinct k-folds \cr
#' 		\code{\link{geoFoldContrast}}: Assign geographically-distinct k-folds to background or absence sites\cr
#' 		\code{\link{elimCellDuplicates}}: Eliminate duplicate points in each cell of a raster \cr
#'
#' @section Bias correction:
#'		\code{\link{geoThin}}: Thin geographic points deterministically or randomly \cr
#'		\code{\link{weightByDist}}: Proximity-based weighting for occurrences (points) to correct for spatial bias \cr
#'
#' @section Model calibration:
#' 		\code{\link{trainByCrossValid}}: and \code{\link{summaryByCrossValid}}: Implement a \code{trainXYZ} function across calibration folds (which are distinct from evaluation folds). \cr
#' 		\code{\link{trainBRT}}: Boosted regression trees (BRTs) \cr
#' 		\code{\link{trainESM}}: Ensembles of small models (ESMs) \cr
#' 		\code{\link{trainGLM}}: Generalized linear models (GLMs) \cr
#' 		\code{\link{trainGLM}}: Generalized linear models (GLMs) \cr
#' 		\code{\link{trainMaxEnt}}: MaxEnt models \cr
#'		\code{\link{trainMaxNet}}: MaxNet models
#' 		\code{\link{trainNS}}: Natural spline (NS) models \cr
#' 		\code{\link{trainRF}}: Random forest (RF) models \cr
#'
#' @section Model prediction:
#' 		\code{\link{predictEnmSdm}}: Predict most model types using default settings \cr
#' 		\code{\link{predictMaxEnt}}: Predict MaxEnt model \cr
#' 		\code{\link{predictMaxNet}}: Predict MaxNet model \cr
#'
#' @section Model evaluation:
#' 		\code{\link{evalAUC}}: AUC (with/out site weights) \cr
#' 		\code{\link{evalMultiAUC}}: Multivariate version of AUC (with/out site weight) \cr
#' 		\code{\link{evalContBoyce}}: Continuous Boyce Index (with/out site weights) \cr
#' 		\code{\link{evalThreshold}}: Thresholds to convert continuous predictions to binary predictions (with/out site weights) \cr
#' 		\code{\link{evalThresholdStats}}: Model performance statistics based on thresholded predictions (with/out site weights) \cr
#' 		\code{\link{evalTjursR2}}: Tjur's R2 (with/out site weights) \cr
#' 		\code{\link{evalTSS}}: True Skill Statistic (TSS) (with/out site weights) \cr
#' 		\code{\link{modelSize}}: Number of response values in a model object \cr
#'
#' @section Functions for rasters:
#' 		\code{\link{bioticVelocity}}: Velocity of movement across a series of rasters \cr
#'		\code{\link{getValueByCell}}: Get value(s) in raster cell(s) by cell number \cr
#'		\code{\link{globalx}}: "Friendly" wrapper for terra::global() for calculating raster statistics \cr
#' 		\code{\link{interpolateRasts}}: Interpolate a stack of rasters \cr
#' 		\code{\link{longLatRasts}}: Generate rasters with values of longitude/latitude for cell values \cr
#' 		\code{\link{sampleRast}}: Sample raster with/out replacement \cr
#'		\code{\link{setValueByCell}}: Set value(s) in raster cell(s) by cell number \cr
#' 		\code{\link{squareCellRast}}: Create a raster with square cells \cr
#'
#' @section Niche overlap and similarity:
#' 		\code{\link{compareResponse}}: Compare niche model responses to a single variable \cr
#' 		\code{\link{nicheOverlapMetrics}}: Niche overlap metrics \cr
#'
#' @section Coordinate reference systems:
#' 		\code{\link{getCRS}}: Return a WKT2 string (coordinate reference system string) using a nickname \cr
#' 		\code{\link{crss}}: Table of coordinate reference systems and nicknames \cr
#' 		\code{\link{customAlbers}}: Create a custom Albers conic equal-area projection \cr
#' 		\code{\link{customLambert}}: Create a custom Lambert azimuthal equal-area projection \cr
#' 		\code{\link{customVNS}}: Create a custom "vertical near-side" projection \cr
#'
#' @section Geographic utility functions:
#'		\code{\link{countPoints}}: Number of points in a "spatial points" object \cr
#' 		\code{\link{decimalToDms}}: Convert decimal coordinate to degrees-minutes-seconds \cr
#' 		\code{\link{dmsToDecimal}}: Convert degrees-minutes-seconds coordinate to decimal \cr
#' 		\code{\link{extentToVect}}: Convert extent to polygon \cr
#'		\code{\link{plotExtent}}: Create a `SpatialPolygon` the same size as a plot region \cr
#'		\code{\link{spatVectorToSpatial}}: Convert \code{SpatVector} object to a \code{Spatial}* object. \cr
#'
#' @section Data:
#' 		\code{\link{canada}}: Outline of Canada \cr
#' 		\code{\link{lemurs}}: Lemur occurrences \cr
#' 		\code{\link{mad0}}: Madagascar spatial object \cr
#' 		\code{\link{mad1}}: Madagascar spatial object \cr
#' 		\code{\link{madClim}}: Madagascar climate rasters for the present \cr
#' 		\code{\link{madClim2030}}: Madagascar climate rasters for the 2030s \cr
#' 		\code{\link{madClim2050}}: Madagascar climate rasters for the 2050s \cr
#' 		\code{\link{madClim2070}}: Madagascar climate rasters for the 2070s \cr
#' 		\code{\link{madClim2090}}: Madagascar climate rasters for the 2090s \cr
#'
#' @docType package
#' @author Adam B. Smith
#' @name enmSdmX
#' @keywords internal 
"_PACKAGE"
NULL
