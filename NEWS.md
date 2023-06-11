# enmSdmX 1.0.8 2023-06-11
- Reworked `trainGLM()` to handle cases with large numbers of predictors.

# enmSdmX 1.0.7 2023-05-22
- Added `geoFoldContrast()` for assigning geo-folds to background or absence sites

# enmSdmX 1.0.6 2023-05-22
- Expanded capacity of .calcWeights() (a hidden function) to handle different values for `family`

# enmSdmX 1.0.5 2023-05-12
- Fixed bug in `predictEnmSdm()` when using a BRT model and predicting to a `SpatRaster` (thank you, Nikki C!)
- Added help page for troubleshooting running functions that support parallel operation

# enmSdmX 1.0.4 2023-04-10

- Fixed failing example in `customCRS()` when GADM server is down
- Fixed bug in `geoThin()` which returned input if it was a `data.frame`*
- Fixed issue in `geoThin()` which returned a `data.frame` lacking coordinates if input was a `data.frame`*
- Users can select clustering method in `geoThin()` and `geoFold()`*

`*` Thank you, Pascal T!

# enmSdmX 1.0.3 2023-03-07

- Removed Fedora issue with external link to `raster` package
- Added `rJava` dependency

# enmSdmX 1.0.2 2023-03-04

- Fixed issue in `nearestGeogPoints()` when polygon lay under centroid
- Fixed bug with "table" call of `getCRS()`
- `getCRS()` (no arguments) now displays `shiny` table of all available CRSs

# enmSdmX 1.0.1

- Fixed bug experienced by some users using `predictEnmSdm()` and `predictMaxNet()` (Thank you, Nikki!)
- Fixed bug in `trainByCrossValid()` using improper call to `evalContBoyce()`
- Fixed `extract()` bug in some examples
- `summarizeByCrossValid()` now summarizes natural spline (NS) models

# enmSdmX 1.0.0

- First release on CRAN
