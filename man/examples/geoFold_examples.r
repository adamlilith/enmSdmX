library(sf)

# lemur occurrence data
data(mad0)
data(lemurs)
crs <- crsGet('WGS84')
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
ll <- c('longitude', 'latitude')
occs <- st_as_sf(occs, coords = ll, crs = crsGet('WGS84'))

folds1 <- pointGeoFold(occs, k = 3)
folds2 <- pointGeoFold(occs, k = 3, minIn = 10)
folds3 <- pointGeoFold(occs, k = 3, minIn = 10, collapseFolds = TRUE)
folds4 <- pointGeoFold(occs, k = 3, minIn = 10, collapseFolds = TRUE, kMin = 2)

# map
par(mfrow=c(2, 2))
plot(st_geometry(occs), cex = 1.4, col = folds1, main = 'Swapping')
plot(st_geometry(mad0), add = TRUE)

plot(st_geometry(occs), cex = 1.4, col = folds2, main = 'Swapping, minIn = 10')
plot(st_geometry(mad0), add = TRUE)

plot(st_geometry(occs), cex = 1.4, col = folds3, main = 'Collapse, minIn = 10')
plot(st_geometry(mad0), add = TRUE)

plot(st_geometry(occs), cex = 1.4, col = folds4, main = 'Collapse --> swapping, minIn = 10')
plot(st_geometry(mad0), add = TRUE)

# number of points per fold
table(folds1) # default
table(folds2) # each fold has >= 10 (minIn)
table(folds3) # can yield just one fold when collapseFolds = TRUE
table(folds4) # collapsing then switch to swapping when k = 2
