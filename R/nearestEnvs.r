#' Extract "most conservative" environments from points and/or polygons
#'
#' This function implements the "nearest environmental point" method (Smith et al. in review) to enable the use of occurrence records geolocated only to a general place (e.g., a country or province), along with occurrences georeferenced with little error.  The function returns environments from a set of precisely-geolocated points plus the environment associated with each imprecise record.
#'
#' @param rasts	A \code{SpatRaster} or "stack" of \code{SpatRaster}s. Please also see argument \code{pca}!
#' @param pts A set of spatial points of class \code{SpatVector} or \code{sf}.
#' @param polys A set of spatial polygons of class \code{SpatVector} or \code{sf}.
#' @param centerFrom Indicates how to locate the "reference" centroid used to identify single points on each polygon. This is only relevant if both \code{pts} and \code{polys} are specified.
#' \itemize{
#' 	\item \code{'pts'}: The default is to use the environmental centroid of \code{pts}, which finds the centroid of \code{pts}, then finds the location on the border of each polygon closest to this centroid.
#' 	\item \code{'polys'}: This option will first calculate the environmental centroid of each polygon, then the centroid of these points, and then find the location on the border of each polygon closest to this point.
#' 	\item \code{'both'}: This option first calculates the environmental centroid of each polygon, then finds the joint centroid of these points plus of \code{pts}, and lastly locates on the border of each polygon the point closest to this grand centroid.
#' }
#' @param pca If \code{TRUE} (default) and there is more than one raster specified in \code{rasts}, then a principal components analysis (PCA) is applied to the values of the rasters before finding the closest points. The returned values are those of the original rasters and the PC scores.
#' @param numPcs The number of PC aces used to find environmental centroids. This is only used if \code{pca} is \code{TRUE}. By default, all axes are used.
#' @param center,scale Settings for \code{\link[stats]{prcomp}}. These indicate if, when calculating the PCA, variables should first be centered and scaled (both \code{TRUE} by default). If the values in \code{rasts} are not of the same units, these should almost always be \code{TRUE}. They are ignored if \code{pca} is \code{FALSE}.
#' @param na.rm If \code{TRUE} (default), ignore \code{NA}s when extracting from rasters (e.g., if a point or polygon falls onto an \code{NA} cell). If \code{FALSE}, then any \code{NA}s that overlap a point or polygon will result in an error.
#' @param return Determines what is returned:
#' \itemize{
#' 		\item \code{'allPoints'} (default): Returns all points. If \emph{n} is the number of points in \code{pts} and \emph{m} the number of polygons in \code{polys}, then the first \code{n} rows in the returned data frame refer to the environments of the \code{pts} and the subsequent \emph{m} to each \code{poly}.
#'		\item \code{'polyPoints'}: Returns the environmental values on each \code{poly} polygon closest to the given center.
#'	}
#'
#' @details This function locates a set of points from the environments covered by each polygon using the following procedure, the details of which depend on what arguments are specified:
#' \itemize{
#' \item Only \code{pts} is specified: Environments are taken directly from the locations of \code{pts}.
#' \item Only \code{polys} is specified: Environments are taken from the closest environment of all the environments associated with each each polygon that is closest to the environmental centroid of the environmental centroids of the polygons.
#' \item \code{pts} and \code{polys} are specified: Environments are taken from the locations of \code{pts} plus the environment from each polygon closest to the environmental centroid of \code{pts}*.
#' }
#'
#' * By default, the function uses the environmental centroid of the precise occurrences in step (1), but this can be changed to the environmental centroid of the centroids of the polygons or the environmental centroid of the points defined by the union of precise occurrence points plus the environmental centroids of the polygons.
#'
#' The function can alternatively return the points on the vertices of the MCP, or points on the input polygons closest to the reference centroid.
#'
#' @return A data frame.
#'
#' @references Smith, A.B., Murphy, S., Henderson, D., and Erickson, K.D. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth. \doi{10.1101/2021.06.10.447988}
#'
#' @seealso \code{\link{mcpFromPointsPolys}} for a related application in geographic space.
#'
#' @examples
#'
#' # This is a contrived example based on red-bellied lemurs in Madagascar
#' # represented by points data and (pretend) Faritras-level occurrences.
#'
#'
#' data(lemurs)
#' redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' ll <- c('longitude', 'latitude')
#' pts <- sf::st_as_sf(redBelly[ , ll], crs=4326, coords=ll)
#'
#' faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
#' 'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
#' 'Analamanga', 'Itasy')
#' data(mad1)
#' polys <- mad1[mad1$NAME_2 %in% faritras, ]
#'
#' rasts <- syste.file('ex/madEnv.tif', package='enmSdmX')
#'
#' # plot environment of points and environments of each polygon closest to
#' # centroid of environments of points
#'
#' envAll <- nearestEnvs(rasts, pts = redBelly, polys = polys, numPcs = 2)
#' envPolys <- nearestEnvs(rasts, pts = redBelly, polys = polys, numPcs = 2,
#' 	return = 'polyPoints')
#' allPolyEnvs <- extract(rasts, polys)
#' plot(envAll$bio1, envAll$bio12, pch=16, col='orange',
#' 	xlab='bio1', ylab='bio12')
#' points(envPolys$bio1, envPolys$bio12, pch=21, bg='orange')
#' legend(
#' 	'topright',
#' 	inset = 0.01,
#' 	legend = c('point', 'polygon'),
#' 	pch = c(16, 21),
#' 	pt.bg = c(NA, 'orange')
#' )
#'
#' # compare to all environments across all polygons
#' allPolyEnvs <- extract(rasts, polys)
#' plot(allPolyEnvs$bio1, allPolyEnvs$bio12, pch=16, col='orange',
#' 	xlab='bio1', ylab='bio12')
#' points(envAll$bio1, envAll$bio12, pch=16)
#' points(envPolys$bio1, envPolys$bio12, pch=21, bg='orange')
#' legend(
#' 	'bottomleft',
#' 	inset = 0.01,
#' 	legend = c('point', 'polygon (closest)', 'polygon (all)'),
#' 	pch = c(16, 21, 16),
#' 	col = c('black', 'black', 'orange'),
#' 	pt.bg = c(NA, 'orange')
#' )
#'
#'
#' @export

nearestEnvs <- function(
	rasts,
	pts = NULL,
	polys = NULL,
	centerFrom = 'pts',
	pca = TRUE,
	numPcs = terra::nlyr(rasts),
	center = TRUE,
	scale = TRUE,
	na.rm = TRUE,
	return = 'allPoints'
) {

	if (inherits(rasts, c('Raster', 'RasterStack', 'RasterBrick'))) {
		warning('Please note that support for the "raster" package will be deprecated in future releases.')
		rasts <- terra::rast(rasts)
	}

	# PCA
	if (pca) {
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
		out <- envPts
	} else if (is.null(pts) & !is.null(polys)) {

		# centroid and distances to centroid
		centroid <- colMeans(envPolys[ , vars], na.rm=TRUE)
		dists <- sweep(envPolys[ , vars, drop=FALSE], 2, centroid)
		dists <- dists^2
		dists <- rowSums(dists)
		# dists <- sqrt(dists) # skipping bc unnecessary to find closest point

		# for each polygon environment, get nearest to centroid
		out <- data.frame()
		ids <- unique(envPolys$ID)
		for (id in ids) {

			theseDists <- dists[envPolys$ID == id]
			closest <- which.min(theseDists)
			thisEnvPolys <- envPolys[envPolys$ID == id, , drop=FALSE]
			nearestPolyEnv <- thisEnvPolys[closest, , drop=FALSE]
			out <- rbind(out, nearestPolyEnv)

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
		centroid <- colMeans(env[ , vars], na.rm=TRUE)

		dists <- sweep(envPolys[ , vars, drop=FALSE], 2L, centroid)
		dists <- dists^2
		dists <- rowSums(dists)
		# dists <- sqrt(dists) # skipping bc unnecessary to find closest point

		# nearest polygon environments to centroid
		nearestPolyEnvs <- data.frame()
		ids <- unique(envPolys$ID)
		for (id in ids) {

			theseDists <- dists[envPolys$ID == id]
			closest <- which.min(theseDists)
			thisEnvPolys <- envPolys[envPolys$ID == id, , drop=FALSE]
			thisNearestPolyEnv <- thisEnvPolys[closest, , drop=FALSE]
			nearestPolyEnvs <- rbind(nearestPolyEnvs, thisNearestPolyEnv)

		}

		out <- if (return == 'allPoints') {
			rbind(envPts, nearestPolyEnvs)
		} else if (return == 'polyPoints') {
			nearestPolyEnvs
		}

	} # if both points and polys exist

	if (na.rm) out <- out[complete.cases(out), , drop=FALSE]
	out

}
