library(sf)

### Madagascar
data(mad0)

alb <- customAlbers(mad0)
lamb <- customLambert(mad0)
vert <- customVNS(mad0)

madAlb <- st_transform(mad0, alb)
madLamb <- st_transform(mad0, lamb)
madVert <- st_transform(mad0, vert)

oldPar <- par(mfrow=c(2, 2))

plot(st_geometry(mad0), main='Unprojected (WGS84)')
plot(st_geometry(madAlb), main='Albers')
plot(st_geometry(madLamb), main='Lambert')
plot(st_geometry(madVert), main='Vertical')

par(oldPar)

### Canada
# The effect is more noticeable when plotting large areas,
# especially if they lie near the poles.
# This example can take a few minutes to run and plot.

library(terra)

canFile <- system.file('extdata', 'canada_level0_gadm41.gpkg', package='enmSdmX')
can <- vect(canFile)

alb <- customAlbers(can)
lamb <- customLambert(can)
vert <- customVNS(can)

canAlb <- project(can, alb)
canLamb <- project(can, lamb)
canVert <- project(can, vert)

oldPar <- par(mfrow=c(2, 2))

plot(can, main = 'Unprojected (WGS84)')
plot(canAlb, main = 'Albers')
plot(canLamb, main = 'Lambert')
plot(canVert, main = 'Vertical')
	
par(oldPar)
