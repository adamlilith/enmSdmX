# enmSdmX
Tools for modeling niches and distributions of species

<img align="right" src="enmSdmX.png"/>

This package contains a suite of efficiency functions for species distribution modeling and ecological niche modeling. You can install this package from CRAN using:

`install.packages('enmSdmX', dependencies = TRUE)`

You can install the development version of this package in R using:

`remotes::install_github('adamlilith/enmSdmX', dependencies=TRUE)`  

### Data preparation ###
* `elimCellDups`: Eliminate duplicate points in each cell of a raster

### Using spatially imprecise records
* `mcpFromPointsPolys`: Minimum convex polygon from a set of spatial polygons and/or points ("nearest geographic point" method)
* `nearestEnvs:  Extract "most conservative" environments from points and/or polygons ("nearest environmental point" method)

### Bias correction
* `pointDistWeights`: Proximity-based weighting for occurrences to correct for spatial bias

### Model training ###
* `trainByCrossValid` and `summaryByCrossValid`: Calibrate a distribution/niche model using cross-validation
* `trainBRT`: Boosted regression trees (BRTs)
* `trainGAM`: Generalized additive models (GAMs)
* `trainGLM`: Generalized linear models (GLMs)
* `trainMaxEnt: MaxEnt models
* `trainMaxNet`: MaxNet models
* `trainNS`: Natural splines (NSs)
* `trainRF`: Random forests (RFs)  

### Model prediction ###
* `predictEnmSdm` Predict most model types using default settings; parallelized
* `predictMaxEnt` Predict MaxEnt model
* `predictMaxNet` Predict MaxNet model

### Model evaluation ###
* `evalAUC`: AUC (with/out site weights)
* `evalMultiAUC`: Multivariate version of AUC (with/out site weight)
* `evalContBoyce`: Continuous Boyce Index (with/out site weights)
* `evalThreshold`: Thresholds to convert continuous predictions to binary predictions (with/out site weights)
* `evalThresholdStats`: Model performance statistics based on thresholded predictions (with/out site weights)
* `evalTjursR2`: Tjur's R2 (with/out site weights)
* `evalTSS`: True Skill Statistic (TSS) (with/out site weights)
* `modelSize`: Number of response values in a model object

### Niche overlap ###
* `compareNiches`: Niche overlap metrics
* `compareResponse`: Compare niche model responses to a single variable

### Functions for rasters ###
* `getValueByCell`: Retrieve raster values(s) by cell number
* `interpolateRasts`: Interpolate a stack of rasters
* `longLatRasts`: Generate rasters with values of longitude/latitude for cell values
* `makeSquareCells`: Create a raster with square cells from an object with an extent
* `rastVelocity`: Velocity of "movement" of mass across a series of rasters
* `sampleRast` : Sample raster with/out replacement
* `setValueByCell`: Set raster values(s) by cell number
* `squareRastCells`: Resample a raster so cells are square

### Range area based on minimum convex polygons ###
* `mcpFromPolys`: Minimum convex polygon from a set of polygons *and* points

### Geographic utility functions ###
* `coordImprecision`: Calculate maximum possible coordinate precision
* `makePlotPoly`: Create a `SpatialPolygon` the same size as a plot region
* `decimalToDms`: Convert decimal coordinate to degrees-minutes-seconds
* `dmsToDecimal`: Convert degrees-minutes-seconds coordinate to decimal
* `extentToVect`: Convert extent to a spatial polygon
* `getCRS`: Return a proj4string (coordinate reference system string) using a nickname
* `svToSpatial`: Convert SpatVector object to a Spatial* object

### Data
* `crss`: Coordinate reference systems
* `lemurs`: Lemur occurrences
* `mad0`: Madagascar spatial object
* `mad1`: Madagascar spatial object
* `madEnv`: Madagascar climate rasters

## Citation ##
As of October 2020 there is no package-specific publication for `enmSdmX`, but the package was first used and cited in:

Morelli*, T.L., Smith*, A.B., Mancini, A.N., Balko, E. A., Borgenson, C., Dolch, R., Farris, Z., Federman, S., Golden, C.D., Holmes, S., Irwin, M., Jacobs, R.L., Johnson, S., King, T., Lehman, S., Louis, E.E. Jr., Murphy, A., Randriahaingo, H.N.T., Lucien, Randriannarimanana, H.L.L., Ratsimbazafy, J., Razafindratsima, O.H., and Baden, A.L. 2020. The fate of Madagascar’s rainforest habitat.  **Nature Climate Change** 10:89-96. * Equal contribution https://doi.org/10.1038/s41558-019-0647-x

**Abstract**. Madagascar has experienced extensive deforestation and overharvesting, and anthropogenic climate change will compound these pressures. Anticipating these threats to endangered species and their ecosystems requires considering both climate change and habitat loss effects. The genus **Varecia** (ruffed lemurs), which is composed of two Critically Endangered forest-obligate species, can serve as a status indicator of the biodiverse eastern rainforest of Madagascar. Here, we combined decades of research to show that the suitable habitat for ruffed lemurs could be reduced by 29–59% from deforestation, 14–75% from climate change (representative concentration pathway 8.5) or 38–93% from both by 2070. If current protected areas avoid further deforestation, climate change will still reduce the suitable habitat by 62% (range: 38–83%). If ongoing deforestation continues, the suitable habitat will decline by 81% (range: 66–93%). Maintaining and enhancing the integrity of protected areas, where rates of forest loss are lower, will be essential for ensuring persistence of the diversity of the rapidly diminishing Malagasy rainforests.
