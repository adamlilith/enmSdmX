#' Return a "proj6" string for a datum and projection
#'
#' This function returns the "proj6" string for a particular datum and possibly projection.
#' @param x Character. Name of proj4 string to return. Spaces, case, and dashes are ignored, but to make the codes more memorable, some of the codes are shown as having them.
#'
#' @section Unprojected:
#' \itemize{
#' 	\item \strong{NAD27}: \code{NAD27}
#' 	\item \strong{NAD83}: \code{NAD83}
#'	\item \strong{WGS72}: \code{wgs72}
#' 	\item \strong{WGS84}: \code{WGS84}
#' }
#'
#' @section Hybrid global/local:
#' \itemize{
#'	\item \strong{Web Mercator}, \strong{Google Web Mercator}, \strong{Spherical Mercator}, \strong{WGS 84 Web Mercator}, or \strong{WGS84/Pseudo-Mercator}: \code{Web Mercator} or \code{WM}
#' }
#'
#' @section Global:
#' \itemize{
#'	\item \strong{Equal Earth}: \code{ee} or \code{Equal Earth}
#'	\item \strong{Eckert IV}: \code{Eckert 4} or \code{Eckert IV}
#'	\item \strong{Good-Homosline (Land)}: \code{GH Land} or \code{Good-Homosline Land}
#' 	\item \strong{Hammer} or \strong{Hammer-Aitoff}: \code{Hammer} or \code{Hammer-Aitoff}
#' 	\item \strong{Hammer-Wagner} (same as Wagner VII): \code{Wagner VII} or {Hammer-Wagner}
#' 	\item \strong{Mollweide}: \code{Mollweide}
#' 	\item \strong{Plate Carree}: \code{Plate Carree} or \code{PC}
#'  \item \strong{Putnin\'s P2} (same as Wagner IV): \code{Putnins P2} or \code{Wagner IV}
#'	\item \strong{Robinson}: \code{Robinson}
#'  \item \strong{Wagner IV} (same as Putnin\'s P2): \code{Wagner IV} or \code{Putnins P2}
#'  \item \strong{Wagner VII} (same as Hammer-Wagner): \code{Wagner VII} or \code{Hammer-Wagner}
#' }
#'
#' @section Continents and Oceans:
#' \itemize{
#'	\item \strong{Africa Albers}: \code{Albers Africa}
#'	\item \strong{Africa Lambert Conformal Conic}: \code{Africa LCC} or \code{Africa Lambert Conformal Conic}

#'	\item \strong{Asia Lambert Conformal Conic}: \code{Asia LCC} or \code{Asia Lambert Conformal Conic}

#'	\item \strong{Australia Albers}: \code{Australia Albers}
#'	\item \strong{Australia Lambert Conformal Conic}: \code{Australia LCC} or \code{Australia Lambert Conformal Conic}

#'	\item \strong{Europe Albers}: \code{Europe Albers}
#'	\item \strong{Europe Lambert Conformal Conic}: \code{Europe LCC} or \code{Europe Lambert Conformal Conic}

#'	\item \strong{North America Albers}: \code{NA Albers} or \code{North America Albers}
#'	\item \strong{North America Lambert Conformal Conic}: \code{NA LCC} or \code{North America Lambert Conformal Conic}

#'	\item \strong{South America Albers}: \code{SA Albers} or \code{South America Albers}
#'	\item \strong{South America Lambert Conformal Conic}: \code{SA LCC} or \code{South America Lambert Conformal Conic}

#' }
#'
#' @section Regional/country:
#' \itemize{
#'	\item \strong{Asia (North) Lambert Conformal Conic}: \code{N Asia LCC} or \code{North Asia Lambert Conformal Conic}
#'	\item \strong{Asia (South) Lambert Conformal Conic}: \code{S Asia LCC} or \code{South Asia Lambert Conformal Conic}
#'	\item \strong{Australia Albers}: \code{Australia Albers}

#'	\item \strong{Central America Albers}: \code{Central America Albers}
#'	\item \strong{Central America Lambert Conformal Conic}: \code{Central America LCC} or \code{Central America Lambert Conformal Conic}

#'	\item \strong{Madagascar Albers}: \code{Madagascar Albers}
#'	\item \strong{Madagascar Lambert Conformal Conic}: \code{Madagascar LCC} or \code{Madagascar Lambert Conformal Conic}

#'	\item \strong{United States (coterminous, or CONUS) Lambert Conformal Conic}: \code{CONUS LCC} or \code{CONUS Lambert Conformal Conic}

#'	\item \strong{North Pole Lambert Azimuthal}: \code{N Pole LA} or \code{North Pole Lambert Azimuthal}
#'	\item \strong{South Pole Lambert Azimuthal}: \code{S Pole LA} or \code{South Pole Lambert Azimuthal}
#' }
#'
#' @section Data-set specific (often same as one of the above, but easy to remember if you are using a specific data set):
#'	\itemize{
#'	\item \strong{ClimateNA}: \code{ClimateNA}
#'	\item \strong{CHELSA}: \code{CHELSA} (actually WGS84)
#'	\item \strong{DayMet}: \code{DayMet} (actually Lambert Conformal Conic for North America)
#'	\item \strong{GBIF}: \code{GBIF} (actually WGS84)
#'	\item \strong{PRISM}: \code{PRISM} (actually WGS84)
#' }

#' @section Universal Trans Mercator (UTM):
#'  UTM projections require the datum (e.g., NAD27, NAD83, WGS84, etc.), the hemisphere (north or south) and the zone. The "code" is as:
#' \code{utm <coordinate system> <hemisphere> <zone>}. For example, \code{UTM NAD83 North 9} returns the CRS for UTM zone 9 in the northern hemisphere using the NAD83 coordinate reference system. Likewise, \code{UTM NAD27 South 6} returns the CRS for UTM zone 6 in the southern hemisphere using the NAD27 coordinate system.
#'
#' @return A WKT (well-known text) object.
#' @examples
#' getCRS('wgs84')
#' getCRS('prism')
#' getCRS('mollweide')
#' @export

getCRS <- function(
	x
) {

	x <- tolower(x)
	x <- gsub(x, pattern=' ', replacement='')
	x <- gsub(x, pattern='-', replacement='')

	out <- if (is.na(x) || x == '') {
		NA

		### world / unprojected
		#######################

		# Eckert IV
		} else if (x %in% c('eckert4', 'eckertiv')) {
			'PROJCS["ProjWiz_Custom_Eckert_IV",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Eckert_IV"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

		# Equal Earth
		} else if (x == 'equalearth') {
			'PROJCS["ProjWiz_Custom_Equal_Earth",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Equal_Earth"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

		# Good Homosline
		} else if (x %in% c('ghland', 'goodhomoslineland')) {
			'PROJCS["World_Goode_Homolosine_Land",
			GEOGCS["WGS 84",
			DATUM["WGS_1984",
			SPHEROID["WGS 84",6378137,298.257223563,
			AUTHORITY["EPSG","7030"]],
			AUTHORITY["EPSG","6326"]],
			PRIMEM["Greenwich",0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Interrupted_Goode_Homolosine"],
			PARAMETER["central_meridian",0],
			PARAMETER["false_easting",0],
			PARAMETER["false_northing",0],
			UNIT["metre",1,
			AUTHORITY["EPSG","9001"]],
			AXIS["Easting",EAST],
			AXIS["Northing",NORTH],
			AUTHORITY["ESRI","54052"]]'

		# Hammer-Aitoff
		} else if (x %in% c('hammer', 'hammeraitoff')) {
			'PROJCS["ProjWiz_Custom_Hammer_Aitoff",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Hammer_Aitoff"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

		# Mollweide
		} else if (x == 'mollweide') {
			'PROJCS["ProjWiz_Custom_Mollweide",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Mollweide"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

		# Wagner IV / Putnins P2
		} else if (x %in% c('wagner4', 'putninsp2')) {
			'PROJCS["ProjWiz_Custom_Wagner_IV",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Wagner_IV"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

		# Wagner VII / Hammer-Wagner
		} else if (x %in% c('wagner7', 'hammerWagner')) {
			'PROJCS["ProjWiz_Custom_Wagner_VII",
			GEOGCS["GCS_WGS_1984",
			DATUM["D_WGS_1984",
			SPHEROID["WGS_1984",6378137.0,298.257223563]],
			PRIMEM["Greenwich",0.0],
			UNIT["Degree",0.0174532925199433]],
			PROJECTION["Wagner_VII"],
			PARAMETER["False_Easting",0.0],
			PARAMETER["False_Northing",0.0],
			PARAMETER["Central_Meridian",0],
			UNIT["Meter",1.0]]'

	### hybrid global/local
	#######################

		# Web Mercator
		} else if (x %in% c('webmercator', 'wm')) {
			'PROJCS["WGS 84 / Pseudo-Mercator",
			GEOGCS["WGS 84",
			DATUM["WGS_1984",
			SPHEROID["WGS 84",6378137,298.257223563,
			AUTHORITY["EPSG","7030"]],
			AUTHORITY["EPSG","6326"]],
			PRIMEM["Greenwich",0,
			AUTHORITY["EPSG","8901"]],
			UNIT["degree",0.0174532925199433,
			AUTHORITY["EPSG","9122"]],
			AUTHORITY["EPSG","4326"]],
			PROJECTION["Mercator_1SP"],
			PARAMETER["central_meridian",0],
			PARAMETER["scale_factor",1],
			PARAMETER["false_easting",0],
			PARAMETER["false_northing",0],
			UNIT["metre",1,
			AUTHORITY["EPSG","9001"]],
			AXIS["X",EAST],
			AXIS["Y",NORTH],
			EXTENSION["PROJ4","+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"],
			AUTHORITY["EPSG","3857"]]'

	### continents / oceans
	#######################

		# Asia Lambert Conformal Conic
		} else if (x %in% c('asialambertconformalconic', 'asialcc')) {
			'PROJCS["Asia_Lambert_Conformal_Conic",
			GEOGCS["GCS_WGS_1984",
			DATUM["WGS_1984",
			SPHEROID["WGS_1984",6378137,298.257223563]],
			PRIMEM["Greenwich",0],
			UNIT["Degree",0.017453292519943295]],
			PROJECTION["Lambert_Conformal_Conic_2SP"],
			PARAMETER["False_Easting",0],
			PARAMETER["False_Northing",0],
			PARAMETER["Central_Meridian",105],
			PARAMETER["Standard_Parallel_1",30],
			PARAMETER["Standard_Parallel_2",62],
			PARAMETER["Latitude_Of_Origin",0],
			UNIT["Meter",1],
			AUTHORITY["EPSG","102012"]]'

		# Europe Lambert Conformal Conic
		} else if (x %in% c('eulcc', 'europelambertconformalconic')) {
			'PROJCS["Europe_Lambert_Conformal_Conic",
			GEOGCS["ED50",
			DATUM["European_Datum_1950",
			SPHEROID["International 1924",6378388,297,
			AUTHORITY["EPSG","7022"]],
			AUTHORITY["EPSG","6230"]],
			PRIMEM["Greenwich",0,
			AUTHORITY["EPSG","8901"]],
			UNIT["degree",0.0174532925199433,
			AUTHORITY["EPSG","9122"]],
			AUTHORITY["EPSG","4230"]],
			PROJECTION["Lambert_Conformal_Conic_2SP"],
			PARAMETER["latitude_of_origin",30],
			PARAMETER["central_meridian",10],
			PARAMETER["standard_parallel_1",43],
			PARAMETER["standard_parallel_2",62],
			PARAMETER["false_easting",0],
			PARAMETER["false_northing",0],
			UNIT["metre",1,
			AUTHORITY["EPSG","9001"]],
			AXIS["Easting",EAST],
			AXIS["Northing",NORTH],
			AUTHORITY["ESRI","102014"]]'

	### regions/countries
	#####################

		# North Pole Lambert Azimuthal
		} else if (x %in% c('npolella', 'North Pole Lambert Azimuthal')) {
			'CONVERSION["North Pole Lambert Azimuthal Equal Area (Bering Sea)",
			METHOD["Lambert Azimuthal Equal Area",
			ID["EPSG",9820]],
			PARAMETER["Latitude of natural origin",90,
			ANGLEUNIT["degree",0.0174532925199433],
			ID["EPSG",8801]],
			PARAMETER["Longitude of natural origin",180,
			ANGLEUNIT["degree",0.0174532925199433],
			ID["EPSG",8802]],
			PARAMETER["False easting",0,
			LENGTHUNIT["metre",1],
			ID["EPSG",8806]],
			PARAMETER["False northing",0,
			LENGTHUNIT["metre",1],
			ID["EPSG",8807]],
			ID["EPSG",17295]]'

		# North Asia Lambert Conformal Conic
		} else if (x %in% c('northasialambertconformalconic', 'nasialcc')) {
			'PROJCS["Asia_North_Lambert_Conformal_Conic",
			GEOGCS["GCS_WGS_1984",
			DATUM["WGS_1984",
			SPHEROID["WGS_1984",6378137,298.257223563]],
			PRIMEM["Greenwich",0],
			UNIT["Degree",0.017453292519943295]],
			PROJECTION["Lambert_Conformal_Conic_2SP"],
			PARAMETER["False_Easting",0],
			PARAMETER["False_Northing",0],
			PARAMETER["Central_Meridian",95],
			PARAMETER["Standard_Parallel_1",15],
			PARAMETER["Standard_Parallel_2",65],
			PARAMETER["Latitude_Of_Origin",30],
			UNIT["Meter",1],
			AUTHORITY["EPSG","102027"]]'

		# South Asia Lambert Conformal Conic
		} else if (x %in% c('southasialambertconformalconic', 'sasialcc')) {
			'PROJCS["Asia_South_Lambert_Conformal_Conic",
			GEOGCS["WGS 84",
			DATUM["WGS_1984",
			SPHEROID["WGS 84",6378137,298.257223563,
			AUTHORITY["EPSG","7030"]],
			AUTHORITY["EPSG","6326"]],
			PRIMEM["Greenwich",0,
			AUTHORITY["EPSG","8901"]],
			UNIT["degree",0.0174532925199433,
			AUTHORITY["EPSG","9122"]],
			AUTHORITY["EPSG","4326"]],
			PROJECTION["Lambert_Conformal_Conic_2SP"],
			PARAMETER["latitude_of_origin",-15],
			PARAMETER["central_meridian",125],
			PARAMETER["standard_parallel_1",7],
			PARAMETER["standard_parallel_2",-32],
			PARAMETER["false_easting",0],
			PARAMETER["false_northing",0],
			UNIT["metre",1,
			AUTHORITY["EPSG","9001"]],
			AXIS["Easting",EAST],
			AXIS["Northing",NORTH],
			AUTHORITY["ESRI","102030"]]'

	### data sets
	#############


	}

	if (asCRS) out <- sp::CRS(out)
	out

}
