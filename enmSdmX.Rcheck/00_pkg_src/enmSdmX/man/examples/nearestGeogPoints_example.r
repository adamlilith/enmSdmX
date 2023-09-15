
library(sf)
library(terra)

#######################################################
### example using SpatVector inputs (terra package) ###
#######################################################

### prepare data
################

# Get coordinate reference systems:
# * WGS84
# * Tananarive (Paris) / Laborde Grid - EPSG:29701
wgs84 <- getCRS('WGS84')
madProj <- getCRS('Madagascar Albers')

# outline of Madagascar faritras
data(mad1)
mad1 <- vect(mad1)
mad1 <- project(mad1, madProj)

# lemur point data
data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- vect(redBelly, geom=ll, crs=wgs84)
redBelly <- project(redBelly, madProj)

# *fake* lemur farita-level data
faritras <- c('Toamasina', 'Atsimo-Atsinana',
'Amoron\'i mania', 'Sava', 'Itasy')
polys <- mad1[mad1$NAME_2 %in% faritras, ]

### apply Nearest Geographic Point method
#########################################

# get three kinds of minimum convex polygons (MCPs):

# MCP using just polygons
mcpPolys <- nearestGeogPoints(polys = polys)

# MCP using just points
mcpPts <- nearestGeogPoints(pts = redBelly)

# MCP using points & polys
mcpPolysPoints <- nearestGeogPoints(pts = redBelly, polys = polys)

# compare extent of occurrence (EOO) in m2
expanse(mcpPolys)
expanse(mcpPts)
expanse(mcpPolysPoints)

### plot minimum convex polygons
################################

# MCP from precise occurrences only
plot(mad1, border='gray', main='MCP points only')
plot(polys, col='gray80', add=TRUE)
plot(mcpPts, col=scales::alpha('red', 0.4), add=TRUE)
plot(redBelly, pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('red', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

# MCP from imprecise occurrences only
plot(mad1, border='gray', main='MCP polys only')
plot(polys, col='gray80', add=TRUE)
plot(mcpPolys, col=scales::alpha('orange', 0.4), add=TRUE)
plot(redBelly, pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('orange', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

# MCP from precise and imprecise occurrences
plot(mad1, border='gray', main='MCP polys + points')
plot(polys, col='gray80', add=TRUE)
plot(mcpPolysPoints, col=scales::alpha('green', 0.4), add=TRUE)
plot(redBelly, pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('green', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

############################################
### example using sf inputs (sf package) ###
############################################

### prepare data
################

# Get coordinate reference systems:
# * WGS84
# * Tananarive (Paris) / Laborde Grid - EPSG:29701
madProj <- sf::st_crs(getCRS('Madagascar Albers'))
wgs84 <- getCRS('WGS84')

# outline of Madagascar faritras
data(mad1)
mad1 <- sf::st_transform(mad1, madProj)

# lemur point occurrence data
data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- sf::st_as_sf(redBelly[ , ll], crs=wgs84, coords=ll)
redBelly <- sf::st_transform(redBelly, madProj)

# *fake* farita-level occurrences
faritras <- c('Toamasina', 'Atsimo-Atsinana',
'Amoron\'i mania', 'Sava', 'Itasy')
polys <- mad1[mad1$NAME_2 %in% faritras, ]

### apply Nearest Geographic Point method
#########################################

# get three kinds of minimum convex polygons (MCPs):

# MCP using just polygons
mcpPolys <- nearestGeogPoints(polys = polys, terra = FALSE)

# MCP using just points
mcpPts <- nearestGeogPoints(pts = redBelly, terra = FALSE)

# MCP using points & polys
mcpPolysPoints <- nearestGeogPoints(pts = redBelly, polys = polys,
terra = FALSE)

# extent of occurrence (EOO) in m2
sf::st_area(mcpPolys)
sf::st_area(mcpPts)
sf::st_area(mcpPolysPoints)

### plot minimum convex polygons
################################

# MCP from precise occurrences only
plot(st_geometry(mad1), border='gray', main='MCP points only')
plot(st_geometry(polys), col='gray80', add=TRUE)
plot(st_geometry(mcpPts), col=scales::alpha('red', 0.4), add=TRUE)
plot(st_geometry(redBelly), pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('red', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

# MCP from imprecise occurrences only
plot(st_geometry(mad1), border='gray', main='MCP points only')
plot(st_geometry(polys), col='gray80', add=TRUE)
plot(st_geometry(mcpPolys), col=scales::alpha('orange', 0.4), add=TRUE)
plot(st_geometry(redBelly), pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('orange', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

# MCP from precise and imprecise occurrences
plot(st_geometry(mad1), border='gray', main='MCP points only')
plot(st_geometry(polys), col='gray80', add=TRUE)
plot(st_geometry(mcpPolysPoints), col=scales::alpha('green', 0.4), add=TRUE)
plot(st_geometry(redBelly), pch=21, bg='red', add=TRUE)

legend('topleft', 
legend=c('Precise occurrence', 'Imprecise occurrence', 'MCP'),
fill=c(NA, 'gray', scales::alpha('green', 0.4)),
pch=c(21, NA, NA),
pt.bg=c('red', NA, NA),
border=c(NA, 'black', 'black'))

### NOTE
# Using SpatVector input (terra package) yields EOOs that are slightly
# larger than using Spatial* (sp) or sf (sf) objects (by about 0.03-0.07%
# in this example). The difference arises because terra::expanse() yields a
# different value than sf::st_area.
