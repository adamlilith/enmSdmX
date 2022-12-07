#' enmSdmX: Species distribution modeling and ecological niche modeling
#'
#' This package contains tools for modeling the distributions and niches of species or species-like entities. Its main features are a set of functions for training and evaluating SDMs/ENMs, and calculating "velocity" of entities represented by a time series of rasters.
#'
#' Create an issue on \href{https://github.com/adamlilith/enmSdmX/issues}{GitHub}.
#'
#' @details
#' @section Data preparation:
#' 		\code{\link{pointElimCellDups}}: Eliminate duplicate points in each cell of a raster \cr
#'
#' @section Using imprecisely-georeferenced occurrences:
#' 		\code{\link{nearestEnvPoints}}: Extract "most conservative" environments from points and/or polygons \cr
#' 		\code{\link{nearestGeogPoints}}: Minimum convex polygon from a set of spatial polygons and/or points \cr
#'
#' @section Bias correction:
#'		\code{\link{pointDistWeights}}: Proximity-based weighting for occurrences (points) to correct for spatial bias \cr
#'		\code{\link{pointGeoThin}}: Deterministic geographic thinning of points \cr
#'
#' @section Model calibration:
#' 		\code{\link{trainByCrossValid}}: and \code{\link{summaryByCrossValid}}: Implement a \code{trainXYZ} function across calibration folds (which are distinct from evaluation folds). \cr
#' 		\code{\link{trainBRT}}: Boosted regression trees (BRTs) \cr
#' 		\code{\link{trainGAM}}: Generalized additive models (GAMs) \cr
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
#'		\code{\link{rastGetValueByCell}}: Get value(s) in raster cell(s) by cell number \cr
#' 		\code{\link{rastInterpolate}}: Interpolate a stack of rasters \cr
#' 		\code{\link{rastLongLat}}: Generate rasters with values of longitude/latitude for cell values \cr
#' 		\code{\link{rastSquareCells}}: Create a raster with square cells \cr
#' 		\code{\link{rastVelocity}}: Velocity of movement across a series of rasters \cr
#' 		\code{\link{rastSample}}: Sample raster with/out replacement \cr
#'		\code{\link{rastSetValueByCell}}: Set value(s) in raster cell(s) by cell number \cr
#'
#' @section Niche overlap and similarity:
#' 		\code{\link{evalNicheOverlap}}: Niche overlap metrics \cr
#' 		\code{\link{compareResponse}}: Compare niche model responses to a single variable \cr
#'
#' @section Coordinate reference systems:
#' 		\code{\link{crss}}: Table of coordinate reference systems \cr
#' 		\code{\link{crsGet}}: Return a WKT2 string (coordinate reference system string) using a nickname \cr
#' 		\code{\link{crsLambert}}: Create a custom Lambert azimuthal equal-area projection \cr
#' 		\code{\link{crsVertical}}: Create a custom "vertical near-side" projection \cr
#'
#' @section Geographic utility functions:
#' 		\code{\link{coordImprecision}}: Coordinate imprecision \cr
#' 		\code{\link{decimalToDms}}: Convert decimal coordinate to degrees-minutes-seconds \cr
#' 		\code{\link{dmsToDecimal}}: Convert degrees-minutes-seconds coordinate to decimal \cr
#' 		\code{\link{extentToVect}}: Convert extent to polygon \cr
#'		\code{\link{plotExtent}}: Create a `SpatialPolygon` the same size as a plot region \cr
#'		\code{\link{spatVectorToSpatial}}: Convert \code{SpatVector} object to a \code{Spatial}* object. \cr
#'
#' @section Data:
#' 		\code{\link{lemurs}}: Lemur occurrences \cr
#' 		\code{\link{mad0}}: Madagascar spatial object \cr
#' 		\code{\link{mad1}}: Madagascar spatial object \cr
#' 		\code{\link{madEnv}}: Madagascar climate rasters \cr
#'
#' @docType package
#' @author Adam B. Smith
#' @name enmSdmX
#' @keywords internal 
"_PACKAGE"
NULL
