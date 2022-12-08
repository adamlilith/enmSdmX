library(sf)

# lemur occurrence data
data(mad0)
data(lemurs)
crs <- crsGet('WGS84')
ll <- c('longitude', 'latitude')

# use occurrences of all species... easier to see on map
occs <- st_as_sf(lemurs, coords = ll, crs = crsGet('WGS84'))

folds2 <- pointGeoFold(occs, k = 2, minIn = 51)
folds4 <- pointGeoFold(occs, k = 4, minIn = 25)

# map folds
par(mfrow = c(1, 2))
plot(st_geometry(occs), pch=folds2, col=folds2, main = '2 g-folds')
plot(st_geometry(mad0), border = 'gray', add = TRUE)

plot(st_geometry(occs), pch=folds4, col=folds4, main = '4 g-folds')
plot(st_geometry(mad0), border = 'gray', add = TRUE)

# inspect number of sites per fold
table(folds2) # 2 folds
table(folds4) # 4 folds
