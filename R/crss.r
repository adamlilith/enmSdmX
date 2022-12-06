#' Coordinate reference systems (CRSs)
#'
#' Outline of Madagascar from GADM
#'
#' @docType data
#'
#' @usage data(crss)
#'
#' @format An object of class \code{data.table}. This is a table with "named" coordinate referenbce systems and their well-known-text (WKT2) representation. It can be used as-is, or with \code{\link{crsGet}}. The fields are as:
#' \itemize{
#'	\item \code{long}: "Long" name of the CRS
#'	\item \code{short1} and \code{short2}: "Short" names of the CRS
#'	\item \code{region}: Region for which CRS is fit
#'	\item \code{projected}: Is the CRS projected or not?
#'	\item \code{projectionGeometry}: Type of projection (\code{NA}, 'cylindrical', 'conic', or 'planar')
#'	\item \code{datum}: Datum
#'	\item \code{type}: Either 'CRS' or 'data'. The former are proper CRSs, and the latter are those used by popular datasets.
#'	\item \code{wkt2}: WKT2 string.
#'	\item \code{notes}: Notes.
#'}
#'
#' @keywords coordinate reference system, CRS
#'
#' @examples
#'
#' data(crss)
#'
'crss'
