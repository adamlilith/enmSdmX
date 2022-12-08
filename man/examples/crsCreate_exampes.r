library(sf)

# Madagascar
data(mad0)

alb <- crsAlbers(mad0)
lamb <- crsLambert(mad0)
vert <- crsVertical(mad0)

madAlb <- st_transform(mad0, alb)
madLamb <- st_transform(mad0, lamb)
madVert <- st_transform(mad0, vert)

par(mfrow=c(2, 2))
plot(st_geometry(mad0), main='Unprojected (WGS84)')
plot(st_geometry(madAlb), main='Albers')
plot(st_geometry(madLamb), main='Lambert')
plot(st_geometry(madVert), main='Vertical')

# The effect is more noticable when plotting large areas,
# especially if they lie near the poles.
if (FALSE) {

library(geodata)
library(terra)

can <- gadm('CAN', level=0, path=tempdir()) # outline of Canada

alb <- crsAlbers(can)
lamb <- crsLambert(can)
vert <- crsVertical(can)

canAlb <- project(can, alb)
canLamb <- project(can, lamb)
canVert <- project(can, vert)

par(mfrow=c(2, 2))
plot(can, main='Unprojected (WGS84)')
plot(canAlb, main='Albers')
plot(canLamb, main='Lambert')
plot(canVert, main='Vertical')
	
}