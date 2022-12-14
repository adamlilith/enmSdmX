% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatVectorToSpatial.r
\name{spatVectorToSpatial}
\alias{spatVectorToSpatial}
\title{Convert SpatVector to Spatial*}
\usage{
spatVectorToSpatial(x)
}
\arguments{
\item{x}{\code{SpatVector} object.}
}
\value{
Object of class \code{Spatial}.
}
\description{
This function converts a \code{SpatVector} object from the \pkg{terra} package to a \code{Spatial} object of the appropriate class (\code{SpatialPoints}, \code{SpatialPointsDataFrame}, \code{SpatialPolygons}, or \code{SpatialPolygonsDataFrame}) from the \pkg{sp} package. Note that \pkg{sp} is to be retired in 2023, so this function is to be come useful only for legacy applications.
}
\examples{

library(terra)
f <- system.file('ex/lux.shp', package='terra')
v <- vect(f)
spat <- spatVectorToSpatial(v)
class(spat)

}
