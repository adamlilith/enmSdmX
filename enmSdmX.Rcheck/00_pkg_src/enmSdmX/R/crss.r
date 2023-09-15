#' Coordinate reference systems (CRSs) and nicknames
#'
#' A table of commonly-used coordinate reference systems, their nicknames, and WKT2 (well-known text) strings
#'
#' @docType data
#'
#' @usage data(crss)
#'
#' @format An object of class \code{data.frame}. This is a table with "named" coordinate referenbce systems and their well-known-text (WKT2) representation. It can be used as-is, or with \code{\link{getCRS}} to quickly get a WKT for a particular CRS. The fields are as:
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
#' @keywords coordinate projection CRS
#'
#' @examples
#'
#' data(crss)
#' getCRS('North America Albers', nice = TRUE)
#'
'crss'
