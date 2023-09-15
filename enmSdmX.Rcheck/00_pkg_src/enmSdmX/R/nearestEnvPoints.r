#' Extract "most conservative" environments from points and/or polygons
#'
#' This function implements the "nearest environmental point" method (Smith et al. 2023) to enable the use of occurrence records geolocated only to a general place (e.g., a country or province), along with occurrences georeferenced with little error.  The function returns environments from a set of precisely-geolocated points plus the environment associated with each imprecise record.
#'
#' @param rasts	A \code{SpatRaster} or "stack" of \code{SpatRaster}s. Please also see argument \code{pca}.
#' @param pts A set of spatial points of class \code{SpatVector} or \code{sf}.
#' @param polys A set of spatial polygons of class \code{SpatVector} or \code{sf}.
#' @param centerFrom Indicates how to locate the "reference" centroid used to identify single points on each polygon. This is only relevant if both \code{pts} and \code{polys} are specified.
#' \itemize{
#' 	\item \code{'pts'}: The default is to use the environmental centroid of \code{pts}, which finds the centroid of \code{pts}, then finds the location on the border of each polygon closest to this centroid.
#' 	\item \code{'polys'}: This option will first calculate the environmental centroid of each polygon, then the centroid of these points, and then find the location on the border of each polygon closest to this point.
#' 	\item \code{'both'}: This option first calculates the environmental centroid of each polygon, then finds the joint centroid of these points plus of \code{pts}, and lastly locates on the border of each polygon the point closest to this grand centroid.
#' }
#' @param pca If \code{TRUE} (default) and there is more than one raster specified in \code{rasts}, then a principal components analysis (PCA) is applied to the values of the rasters before finding the closest points. The returned values are those of the original rasters and the PC scores.
#' @param numPcs The number of PC axes used to find environmental centroids. This is only used if \code{pca} is \code{TRUE}. By default, all axes are used.
#' @param center,scale Settings for \code{\link[stats]{prcomp}}. These indicate if, when calculating the PCA, variables should first be centered and scaled (both \code{TRUE} by default). If the values in \code{rasts} are not of the same units, this should almost always be \code{TRUE}. They are ignored if \code{pca} is \code{FALSE}.
#' @param rule Determines how to identify the single environmental point to associate with each polygon. Options include:
#' \itemize{
#'	\item	\code{'nearest'} (default): Returns the environmental point \emph{closest} to the centroid (i.e., the "nearest environmental point").
#'	\item	\code{'farthest'}: Returns the environmental point \emph{farthest} from the centroid (i.e., opposite of the "nearest" point)
#' }
#' @param na.rm If \code{TRUE} (default), ignore \code{NA}s when extracting from rasters (e.g., if a point or polygon falls onto an \code{NA} cell). If \code{FALSE}, then any \code{NA}s that overlap a point or polygon will result in an error.
#' @param out Determines what is returned. Only used if both \code{pts} and \code{polys} are provided.
#' \itemize{
#' 		\item \code{'both'} (default): Returns all environmental points. If \emph{n} is the number of points in \code{pts} and \emph{m} the number of polygons in \code{polys}, then the first \code{n} rows in the returned data frame refer to the environments of the \code{pts} and the subsequent \emph{m} to each \code{poly}.
#'		\item \code{'pts'}: Returns the environmental values associated with each point.
#'		\item \code{'polys'}: Returns the environmental values on each \code{poly} polygon closest to the given center.
#'	}
#'
#' @details This function locates a set of points from the environments covered by each polygon using the following procedure, the details of which depend on what arguments are specified:
#' \itemize{
#' \item Only \code{pts} is specified: Environments are taken directly from the locations of \code{pts} in environmental space.
#' \item Only \code{polys} is specified: Environments are taken from the closest environment of all the environments associated with each each polygon that is closest to the environmental centroid of the environmental centroids of the polygons (that may be confusing, but it is not a typo).
#' \item \code{pts} and \code{polys} are specified: Environments are taken from the locations of \code{pts} plus the environment from each polygon closest to the environmental centroid of \code{pts}. By default, the function uses the environmental centroid of the precise occurrences in step (1), but this can be changed to the environmental centroid of the centroids of the polygons or the environmental centroid of the points defined by the union of precise occurrence points plus the environmental centroids of the polygons.
#' }
#'
#' The function can alternatively return the points on the vertices of the MCP, or points on the input polygons closest to the reference centroid.
#'
#' @return A data frame.
#'
#' @references
#' Smith, A.B., Murphy, S.J., Henderson, D., and Erickson, K.D. 2023. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth.  \emph{Global Ecology and Biogeography} In press. Open access pre-print: \doi{10.1101/2021.06.10.447988}
#'
#' @seealso \code{\link{nearestGeogPoints}} for the "nearest geographic point" method, a related approach for geographic space.
#'
#' @example man/examples/nearestEnvPoints_example.r
#'
#' @export
nearestEnvPoints <- function(
	rasts,
	pts = NULL,
	polys = NULL,
	centerFrom = 'pts',
	pca = TRUE,
	numPcs = 3,
	center = TRUE,
	scale = TRUE,
	rule = 'nearest',
	na.rm = TRUE,
	out = 'both'
) {

	if (inherits(rasts, c('Raster', 'RasterStack', 'RasterBrick'))) {
		warning('Please note that support for the "raster" package will be deprecated in future releases.')
		rasts <- terra::rast(rasts)
	}

	# PCA
	if (pca) {
		numPcs <- min(numPcs, terra::nlyr(rasts))

		env <- as.data.frame(rasts)
		pcModel <- stats::prcomp(env, center=center, scale.=scale)
		rastsPcs <- terra::predict(rasts, pcModel)
		rasts <- c(rasts, rastsPcs)

		vars <- paste0('PC', 1L:numPcs)
	} else {
		vars <- names(rasts)
	}

	if (!is.null(pts)) {
		if (!inherits(pts, 'SpatVector')) pts <- terra::vect(pts)
		envPts <- terra::extract(rasts, pts)
		envPts$ID <- paste0('point ', envPts$ID)
	}
	if (!is.null(polys)) {
		if (!inherits(polys, 'SpatVector')) polys <- terra::vect(polys)
		envPolys <- terra::extract(rasts, polys)
		envPolys$ID <- paste0('poly ', envPolys$ID)
	}

	### calculate centroid
	if (is.null(pts) & is.null(polys)) {
		stop('Either "polys", "pts", or both must be specified.')
	} else if (!is.null(pts) & is.null(polys)) {
		output <- envPts
	} else if (is.null(pts) & !is.null(polys)) {

		# centroid and distances to centroid
		centroid <- colMeans(envPolys[ , vars], na.rm = na.rm)
		dists <- sweep(envPolys[ , vars, drop=FALSE], 2, centroid)
		dists <- dists^2
		dists <- rowSums(dists)
		# dists <- sqrt(dists) # skipping bc unnecessary to find closest point

		# for each polygon environment, get nearest to centroid
		output <- data.frame()
		ids <- unique(envPolys$ID)
		for (id in ids) {

			theseDists <- dists[envPolys$ID == id]
			closest <- if (rule == 'nearest') {
				which.min(theseDists)
			} else if (rule == 'farthest') {
				which.max(theseDists)
			} else {
				stop('Argument "rule" must be either "nearest" or "farthest".')
			}
			thisEnvPolys <- envPolys[envPolys$ID == id, , drop=FALSE]
			nearestPolyEnv <- thisEnvPolys[closest, , drop=FALSE]
			output <- rbind(output, nearestPolyEnv)

		}

	} else {

		env <- if (centerFrom == 'pts') {
			envPts
		} else if (centerFrom == 'polys') {
			envPolys
		} else {
			rbind(envPts, envPolys)
		}

		# centroid and distances to centroid
		centroid <- colMeans(env[ , vars], na.rm = na.rm)

		dists <- sweep(envPolys[ , vars, drop=FALSE], 2L, centroid)
		dists <- dists^2
		dists <- rowSums(dists)
		# dists <- sqrt(dists) # skipping bc unnecessary to find closest point

		# nearest polygon environments to centroid
		nearestPolyEnvs <- data.frame()
		ids <- unique(envPolys$ID)
		for (id in ids) {

			theseDists <- dists[envPolys$ID == id]
			closest <- if (rule == 'nearest') {
				which.min(theseDists)
			} else if (rule == 'farthest') {
				which.max(theseDists)
			} else {
				stop('Argument "rule" must be either "nearest" or "farthest".')
			}
			thisEnvPolys <- envPolys[envPolys$ID == id, , drop=FALSE]
			thisNearestPolyEnv <- thisEnvPolys[closest, , drop=FALSE]
			nearestPolyEnvs <- rbind(nearestPolyEnvs, thisNearestPolyEnv)

		}

		output <- if (out == 'both') {
			rbind(envPts, nearestPolyEnvs)
		} else if (out == 'polys') {
			nearestPolyEnvs
		} else if (out == 'pts') {
			envPts
		} else {
			stop('Argument "out" must be one of "both", "polys", or "pts".')
		}

	} # if both points and polys exist

	if (na.rm) output <- output[stats::complete.cases(output), , drop=FALSE]
	rownames(output) <- NULL
	output

}
