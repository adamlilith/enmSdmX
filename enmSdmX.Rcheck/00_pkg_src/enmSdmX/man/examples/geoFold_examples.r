library(sf)
library(terra)

# lemur occurrence data
data(mad0)
data(lemurs)
crs <- getCRS('WGS84')
ll <- c('longitude', 'latitude')

# use occurrences of all species... easier to see on map
occs <- st_as_sf(lemurs, coords = ll, crs = getCRS('WGS84'))

# create 100 background points
mad0 <- vect(mad0)
bg <- spatSample(mad0, 100)

### assign 3 folds to occurrences and to background sites
k <- 3
minIn <- floor(nrow(occs) / k) # maximally spread between folds

presFolds <- geoFold(occs, k = k, minIn = minIn)
bgFolds <- geoFoldContrast(bg, pres = occs, presFolds = presFolds)

# number of sites per fold
table(presFolds)
table(bgFolds)

# map
plot(mad0, border = 'gray', main = paste(k, 'geo-folds'))
plot(bg, pch = 3, col = bgFolds + 1, add = TRUE)
plot(st_geometry(occs), pch = 20 + presFolds, bg = presFolds + 1, add = TRUE)

legend(
	'bottomright',
	legend = c(
		'presence fold 1',
		'presence fold 2',
		'presence fold 3',
		'background fold 1',
		'background fold 2',
		'background fold 3'
	),
	pch = c(21, 22, 23, 3, 3),
	col = c(rep('black', 3), 2, 3),
	pt.bg = c(2, 3, 4, NA, NA)
)
