#' @name mad1
#'
#' @title Madagascar spatial object
#'
#' @description Outlines of regions ("Faritra") of Madagascar from GADM. The geometry has been simplified from the version available in GADM, so pleased do not use this for "official" analyses.
#'
#' @docType data
#'
#' @usage data(mad1, package='enmSdmX')
#'
#' @format An object of class \code{sf}.
#'
#' @keywords Madagascar
#'
#' @source \href{https://gadm.org}{GADM}
#' 
#' @examples
#'
#' library(sf)
#' data(mad1)
#' mad1
#' plot(st_geometry(mad1), main='Madagascar')
#'
NULL
