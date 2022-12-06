
library(sf)
library(terra)

### example using SpatVector inputs (terra package)
###################################################

# Tananarive (Paris) / Laborde Grid - EPSG:29701
wgs84 <- crsGet('WGS84')
madProj <- crsGet('Madagascar Albers')

data(mad1)
mad1 <- vect(mad1)
mad1 <- project(mad1, madProj)

data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- vect(redBelly[ , ll], geom=ll, crs=wgs84)
redBelly <- project(redBelly, madProj)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
polys <- mad1[mad1$NAME_2 %in% faritras, ]

mcpPolys <- nearestGeogPoints(polys = polys)
mcpPts <- nearestGeogPoints(pts = redBelly, polys = NULL)
mcpPolysPoints <- nearestGeogPoints(pts = redBelly, polys = polys)

# extent of occurrence in m2
expanse(mcpPolys)
expanse(mcpPts)
expanse(mcpPolysPoints)

plot(mad1, border='gray')
plot(polys, col='gray80', add=TRUE)
plot(mcpPolysPoints, col=scales::alpha('green', 0.4), add=TRUE)
plot(mcpPolys, col=scales::alpha('purple', 0.4), add=TRUE)
plot(mcpPts, col=scales::alpha('red', 0.4), add=TRUE)
plot(redBelly, pch=16, add=TRUE)

legend('topleft', 
legend=c('Presences', '"Occupied" Faritras',
'MCP w/ polygons', 'MCP w/ points', 'MCP w/ polygons & points'),
fill=c(NA, 'gray', scales::alpha('purple', 0.4),
scales::alpha('red', 0.4),
scales::alpha('green', 0.4)),
pch=c(16, NA, NA, NA, NA),
border=c(NA, 'black', 'black', 'black', 'black'))

### example using sf* inputs (sf package)
#########################################

# Tananarive (Paris) / Laborde Grid - EPSG:29701
madProj <- sf::st_crs(crsGet('Madagascar Albers'))
wgs84 <- crsGet('WGS84')

data(mad1)
mad1 <- sf::st_transform(mad1, madProj)

data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- sf::st_as_sf(redBelly[ , ll], crs=wgs84, coords=ll)
redBelly <- sf::st_transform(redBelly, madProj)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
polys <- mad1[mad1$NAME_2 %in% faritras, ]

mcpPolys <- nearestGeogPoints(polys = polys, terra = FALSE)
mcpPts <- nearestGeogPoints(pts = redBelly, polys = NULL, terra = FALSE)
mcpPolysPoints <- nearestGeogPoints(pts = redBelly, polys = polys,
terra = FALSE)

# extent of occurrence in m2
sf::st_area(mcpPolys)
sf::st_area(mcpPts)
sf::st_area(mcpPolysPoints)

plot(sf::st_geometry(mad1))
plot(sf::st_geometry(polys), col='gray80', add=TRUE)
plot(sf::st_geometry(mcpPolysPoints), col=scales::alpha('green', 0.4),
add=TRUE)
plot(mcpPts, col=scales::alpha('red', 0.4), add=TRUE)
plot(mcpPolys, col=scales::alpha('purple', 0.4), add=TRUE)
plot(redBelly, pch=16, add=TRUE)

legend('topleft', 
legend=c('Presences', '"Occupied" Faritras',
'MCP w/ polygons', 'MCP w/ points', 'MCP w/ polygons & points'),
fill=c(NA, 'gray', scales::alpha('purple', 0.4),
scales::alpha('red', 0.4),
scales::alpha('green', 0.4)),
pch=c(16, NA, NA, NA, NA),
border=c(NA, 'black', 'black', 'black', 'black'))

### NOTE
# Using SpatVector input (terra package) yields EOOs that are slightly
# larger than using Spatial* (sp) or sf (sf) objects (by about 0.03-0.07%
# in this example). The difference arises because terra::expanse() yields a
# different value than sf::st_area.
