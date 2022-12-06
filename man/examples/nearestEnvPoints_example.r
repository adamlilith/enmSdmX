# This is a contrived example based on red-bellied lemurs in Madagascar.
# Point locations (which are real data) will be assumed to be "precise"
# records. We will designate a set of Faritas ("counties") to represent
# "imprecise" occurrences that can only be georeferenced to a geopolitical
# unit.

library(sf)
library(terra)

data(lemurs)
precise <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
wgs84 <- crsGet('WGS84')
precise <- sf::st_as_sf(precise[ , ll], coords=ll, crs=wgs84)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
data(mad1)
imprecise <- mad1[mad1$NAME_2 %in% faritras, ]

rastFile <- system.file('extdata/madEnv.tif', package='enmSdmX')
rasts <- rast(rastFile)

### Plot environment of points and environments of each polygon closest to
### centroid of environments of points. In this example, we use the first two
### principal component axes to characterize the niche.

envPtsPolys <- nearestEnvPoints(rasts, pts = precise, polys = imprecise,
	pca = TRUE,	numPcs = 2)
envPolys <- nearestEnvPoints(rasts, pts = precise, polys = imprecise, numPcs = 2,
	out = 'polys')
envPts <- nearestEnvPoints(rasts, pts = precise, polys = imprecise, numPcs = 2,
	out = 'pts')
allPolyEnvs <- extract(rasts, imprecise)
plot(envPtsPolys$PC1, envPtsPolys$PC2, pch=16, col='black',
	xlab='PC1', ylab='PC2')
points(envPolys$PC1, envPolys$PC2, pch=21, bg='orange')
legend(
	'bottomleft',
	inset = 0.01,
	legend = c('precise', 'imprecise (closest)'),
	pch = c(16, 21),
	col = c('black', 'black'),
	pt.bg = c('orange', 'orange')
)

### compare identified environments to all environments across all polygons
###########################################################################
env <- as.data.frame(rasts)
pca <- stats::prcomp(env, center=TRUE, scale.=TRUE)

allPolyEnvs <- extract(rasts, imprecise, ID = FALSE)
allPolyEnvsPcs <- predict(pca, allPolyEnvs)
allPolyEnvs <- cbind(allPolyEnvs, allPolyEnvsPcs)

plot(allPolyEnvs$PC1, allPolyEnvs$PC2, pch=16, col='orange',
	xlab='PC1', ylab='PC2')
points(envPts$PC1, envPts$PC2, pch=16)
points(envPolys$PC1, envPolys$PC2, pch=1)
legend(
	'bottomleft',
	inset = 0.01,
	legend = c('precise', 'imprecise (closest)', 'imprecise (all)'),
	pch = c(16, 21, 16),
	col = c('black', 'black', 'orange'),
	pt.bg = c(NA, 'orange')
)

### display niches (minimum convex hulls) estimated
### using just precise or precise + imprecise records
#####################################################
pcs <- c('PC1', 'PC2')
preciseIndices <- chull(envPts[ , pcs])
preciseImpreciseIndices <- chull(envPtsPolys[ , pcs])

preciseIndices <- c(preciseIndices, preciseIndices[1])
preciseImpreciseIndices <- c(preciseImpreciseIndices,
	preciseImpreciseIndices[1])

preciseOnlyNiche <- envPts[preciseIndices, pcs]
preciseImpreciseNiche <- envPtsPolys[preciseImpreciseIndices, pcs]

# plot
plot(allPolyEnvs$PC1, allPolyEnvs$PC2, pch=16, col='orange',
	xlab='PC1', ylab='PC2')
points(envPts$PC1, envPts$PC2, pch=16)
points(envPolys$PC1, envPolys$PC2, pch=1)
lines(preciseImpreciseNiche, col='coral4', lwd=2)
lines(preciseOnlyNiche, lty='dotted')

legend(
	'bottomleft',
	inset = 0.01,
	legend = c('precise', 'imprecise (closest)', 'imprecise (all)', 
		'MCP imprecise-only', 'MCP precise + imprecise'),
	pch = c(16, 21, 16, NA, NA),
	col = c('black', 'black', 'orange', 'black', 'coral4'),
	pt.bg = c(NA, 'orange', NA, NA, NA),
	lwd = c(NA, NA, NA, 1, 2),
	lty = c(NA, NA, NA, 'dotted', 'solid')
)
