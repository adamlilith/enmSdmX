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
#' @source \href{www.gadm.org}{GADM}
#' 
#' @examples
#' data(mad1)
#' mad1
#' plot(mad1['NAME_2'], main='Malagasy Faritra')
NULL
