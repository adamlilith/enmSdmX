# enmSdmX 1.0.3 2023-03-07

- Removed Fedora issue with external link to `raster` package

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
