\donttest{
# NB These examples can take a few minutes to run.
# To illustrate calculation and interpretation of biotic velocity,
# we will calibrate a SDM for the Red-Bellied Lemur and project
# the model to the present and successive future climates. The time series
# of rasters is then used to calculate biotic velocity.

library(sf)
library(terra)

### process environmental rasters
#################################

# get rasters
rastFile <- system.file('extdata/madClim.tif', package='enmSdmX')
madClim <- rast(rastFile)

rastFile <- system.file('extdata/madClim2030.tif', package='enmSdmX')
madClim2030 <- rast(rastFile)

rastFile <- system.file('extdata/madClim2050.tif', package='enmSdmX')
madClim2050 <- rast(rastFile)

rastFile <- system.file('extdata/madClim2070.tif', package='enmSdmX')
madClim2070 <- rast(rastFile)

rastFile <- system.file('extdata/madClim2090.tif', package='enmSdmX')
madClim2090 <- rast(rastFile)

# The bioticVelocity() function needs rasters to be in equal-area
# projection, so we will project them here.
madAlbers <- getCRS('madAlbers') # Albers projection for Madagascar
madClim <- project(madClim, madAlbers)
madClim2030 <- project(madClim2030, madAlbers)
madClim2050 <- project(madClim2050, madAlbers)
madClim2070 <- project(madClim2070, madAlbers)
madClim2090 <- project(madClim2090, madAlbers)

# Coordinate reference systems:
wgs84 <- getCRS('WGS84') # WGS84
madAlbers <- getCRS('madAlbers') # Madagascar Albers

# lemur occurrence data
data(lemurs)
occs <- lemurs[lemurs$species == 'Eulemur fulvus', ]
occs <- vect(occs, geom=c('longitude', 'latitude'), crs=wgs84)
occs <- project(occs, madAlbers)

# eliminate cell duplicates
occs <- elimCellDuplicates(occs, madClim)

# extract environment at occurrences
occEnv <- extract(madClim, occs, ID = FALSE)
occEnv <- occEnv[complete.cases(occEnv), ]
	
# create background sites (using just 1000 to speed things up!)
bgEnv <- terra::spatSample(madClim, 3000)
bgEnv <- bgEnv[complete.cases(bgEnv), ]
bgEnv <- bgEnv[sample(nrow(bgEnv), 1000), ]

# collate occurrences and background sites
presBg <- data.frame(
   presBg = c(
      rep(1, nrow(occEnv)),
      rep(0, nrow(bgEnv))
   )
)

env <- rbind(occEnv, bgEnv)
env <- cbind(presBg, env)

### calibrate model
###################

predictors <- c('bio1', 'bio12')

# MaxEnt
mx <- trainMaxEnt(
	data = env,
	resp = 'presBg',
	preds = predictors,
	regMult = 1, # too few values for reliable model, but fast
	cores = 2
)

### project to present and future climate
#########################################

predPresent <- predictEnmSdm(mx, madClim)
pred2030 <- predictEnmSdm(mx, madClim2030)
pred2050 <- predictEnmSdm(mx, madClim2050)
pred2070 <- predictEnmSdm(mx, madClim2070)
pred2090 <- predictEnmSdm(mx, madClim2090)

plot(predPresent, main = 'Present Suitability')

# plot change in suitability between present and 2090s
delta <- pred2090 - predPresent
plot(delta, main = 'Change in Suitability')

### calculate biotic velocity
#############################

series <- c(
	predPresent,
	pred2030,
	pred2050,
	pred2070,
	pred2090
)

names(series) <- c('present', 't2030', 't2050', 't2070', 't2090')
plot(series)

times <- c(1985, 2030, 2050, 2070, 2090)
quants <- c(0.10, 0.90)

bv <- bioticVelocity(
	x = series,
	times = times,
	quants = quants,
	cores = 2
)
 
bv

### centroid velocities

# centroid (will always be >= 0)
# fastest centroid movement around 2060
plot(bv$timeMid, bv$centroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Speed (m / y)', main = 'Centroid Speed')
  
# velocity northward/southward through time
# shows northward shift because always positive, fastest around 2060
plot(bv$timeMid, bv$nsCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)', main = 'Centroid N/S Velocity')
  
# velocity eastward (positive)/westward (negative) through time
# movement eastward (positive) first, then westward (negative)
plot(bv$timeMid, bv$ewCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)', main = 'Centroid E/W Velocity')

### map of centroid location through time
# shows centroid moves slightly northward through time
plot(delta, main = 'Centroid Location &\nChange in Suitability')
points(bv$centroidLong[1], bv$centroidLat[1], pch = 1)
points(bv$centroidLong[4], bv$centroidLat[4], pch = 16)
lines(bv$centroidLong, bv$centroidLat)
legend(
  'bottomright',
  legend = c(
    'start (~1985)',
	'stop (~2090)',
	'trajectory'
  ),
  pch = c(1, 16, NA),
  lwd = c(NA, NA, 1)
)

### velocities of portions of range north/south/east/west of centroid
# positive ==> northward shift
# negative ==> southward shift
  
# portion of range north of centroid
# shows northward expansion because always positive
plot(bv$timeMid, bv$nCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Northern Part of Range')
  
# portion of range south of centroid
# shows northward contraction because always positive
plot(bv$timeMid, bv$sCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Southern Part of Range')
  
# portion of range east of centroid
# shows eastern portion moves farther east
plot(bv$timeMid, bv$eCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Eastern Part of Range')
  
# portion of range west of centroid
# shows western portion moves east
plot(bv$timeMid, bv$wCentroidVelocity, type = 'l',
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Western Part of Range')

### velocities of range margins

# from south to north, 10th and 90th quantiles of density
# positive ==> northward shift
# negative ==> southward shift
# shows both northern and southern range margins shift northward
# because always positive... northern margin shift usually slower
ylim <- range(bv$nsQuantVelocity_quant0p1, bv$nsQuantVelocity_quant0p9)

plot(bv$timeMid, bv$nsQuantVelocity_quant0p1, type = 'l', ylim = ylim,
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Northern/Southern Range Margins')
lines(bv$timeMid, bv$nsQuantVelocity_quant0p9, lty = 'dashed')
legend(
  'bottomright',
  legend = c('Southern Margin', 'Northern Margin'),
  lty = c('solid', 'dashed')
)
  
# from east to west, 10th and 90th quantiles of density
# positive ==> eastward shift
# negative ==> westward shift
ylim <- range(bv$ewQuantVelocity_quant0p1, bv$ewQuantVelocity_quant0p9)

plot(bv$timeMid, bv$ewQuantVelocity_quant0p1, type = 'l', ylim = ylim,
  xlab = 'Year', ylab = 'Velocity (m / y)',
  main = 'Eastern/Western Range Margins')
lines(bv$timeMid, bv$ewQuantVelocity_quant0p9, lty = 'dashed')
legend(
  'bottomright',
  legend = c('Eastern Margin', 'Western Margin'),
  lty = c('solid', 'dashed')
)
  
  
### summary statistics

# mean density across cells through time
plot(bv$timeMid, bv$mean, type = 'l',
  xlab = 'Year', ylab = 'Mean Density',
  main = 'Mean Density')

# sum of density across cells through time
plot(bv$timeMid, bv$sum, type = 'l',
  xlab = 'Year', ylab = 'Sum of Density',
  main = 'Sum of Density')

### change metrics

# average change in suitability from one time period to next
# shows average conditions getting worse
plot(bv$timeMid, bv$simpleMeanDiff, type = 'l',
  xlab = 'Year', ylab = 'Mean Change in Suitability')
  
# average absolute change in suitability from one time period to next
# shows average absolute change declining
plot(bv$timeMid, bv$meanAbsDiff, type = 'l',
  xlab = 'Year', ylab = 'Mean Absolute Change in Suitability')
  
# root-mean square difference from one time period to the next
# shows difference between successive rasters declines through time
plot(bv$timeMid, bv$rmsd, type = 'l',
  xlab = 'Year', ylab = 'RMSD')
  
### raster similarity
# most indicate that successive rasters are similar through time
ylim <- range(bv$godsoeEsp, bv$schoenerD, bv$warrenI, bv$cor, bv$warrenI)
plot(bv$timeMid, bv$godsoeEsp, type = 'l', lty = 1, col = 1,
  xlab = 'Year', ylab = 'Raster similarity', ylim = ylim)
lines(bv$timeMid, bv$schoenerD, lty = 2, col = 2)
lines(bv$timeMid, bv$warrenI, lty = 3, col = 3)
lines(bv$timeMid, bv$cor, lty = 4, col = 4)
lines(bv$timeMid, bv$rankCor, lty = 5, col = 5)

legend(
  'right',
  legend = c(
    'Godsoe\'s ESP',
    'Schoener\'s D',
    'Warren\'s I',
    'Correlation',
    'Rank Correlation'
  ),
  col = 1:5,
  lty = 1:5
)
  

# values of 10th and 90th quantiles across cells through time
# shows most favorable cells becoming less favorable
# least favorable cells remain mainly unchanged
ylim <- range(bv$quantile_quant0p1, bv$quantile_quant0p9)

plot(bv$timeMid, bv$quantile_quant0p1, type = 'l', ylim = ylim,
  xlab = 'Year', ylab = 'Quantile Value',
  main = 'Quantiles across Cells')
lines(bv$timeMid, bv$quantile_quant0p9, lty = 'dashed')

legend(
  'topright',
  legend = c('10th quantile', '90th quantile'),
  lty = c('solid', 'dashed')
)

### map of northern/southern range margins through time

# range of longitude shown in plot
madExtent <- ext(madClim)
xExtent <- as.vector(madExtent)[1:2]

plot(predPresent, main = 'North/South Range Margin Location')
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p9[1], bv$nsQuantLat_quant0p9[1]))
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p9[2], bv$nsQuantLat_quant0p9[2]), lty = 'dashed')
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p9[3], bv$nsQuantLat_quant0p9[3]), lty = 'dotdash')
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p9[4], bv$nsQuantLat_quant0p9[4]), lty = 'dotted')

lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p1[1], bv$nsQuantLat_quant0p1[1]))
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p1[2], bv$nsQuantLat_quant0p1[2]), lty = 'dashed')
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p1[3], bv$nsQuantLat_quant0p1[3]), lty = 'dotdash')
lines(c(xExtent[1], xExtent[2]),
  c(bv$nsQuantLat_quant0p1[4], bv$nsQuantLat_quant0p1[4]), lty = 'dotted')

legend(
  'bottomright',
  legend = c(
    '1980s',
	'2030s',
	'2050s',
	'2070s',
	'2090s'
  ),
  lty = c('solid', 'dashed', 'dotdash', 'dotted')
)

### map of eastern/western range margins through time

# range of longitude shown in plot
madExtent <- ext(madClim)
yExtent <- as.vector(madExtent)[3:4]

plot(predPresent, main = 'North/South Range Margin Location')
lines(c(bv$ewQuantLong_quant0p9[1], bv$ewQuantLong_quant0p9[1]),
  c(yExtent[1], yExtent[2]))
lines(c(bv$ewQuantLong_quant0p9[2], bv$ewQuantLong_quant0p9[2]),
  c(yExtent[1], yExtent[2]), lty = 'dashed')
lines(c(bv$ewQuantLong_quant0p9[3], bv$ewQuantLong_quant0p9[3]),
  c(yExtent[1], yExtent[2]), lty = 'dotdash')
lines(c(bv$ewQuantLong_quant0p9[4], bv$ewQuantLong_quant0p9[4]),
  c(yExtent[1], yExtent[2]), lty = 'dotted')

lines(c(bv$ewQuantLong_quant0p1[1], bv$ewQuantLong_quant0p1[1]),
  c(yExtent[1], yExtent[2]))
lines(c(bv$ewQuantLong_quant0p1[2], bv$ewQuantLong_quant0p1[2]),
  c(yExtent[1], yExtent[2]), lty = 'dashed')
lines(c(bv$ewQuantLong_quant0p1[3], bv$ewQuantLong_quant0p1[3]),
  c(yExtent[1], yExtent[2]), lty = 'dotdash')
lines(c(bv$ewQuantLong_quant0p1[4], bv$ewQuantLong_quant0p1[4]),
  c(yExtent[1], yExtent[2]), lty = 'dotted')

legend(
  'bottomright',
  legend = c(
    '1980s',
	'2030s',
	'2050s',
	'2070s',
	'2090s'
  ),
  lty = c('solid', 'dashed', 'dotdash', 'dotted')
)

}

