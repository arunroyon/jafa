# Outputs

## Feature Maps

For each feature, the package writes:

- `<feature>_map.fits`: continuum-subtracted integrated feature map.
- `<feature>_unc.fits`: integrated uncertainty map when an uncertainty cube is available.
- `<feature>_continuum.fits`: integrated fitted continuum over the same feature window.
- `<feature>_diagnostic.png`: feature map, uncertainty map, and median-spectrum continuum diagnostic.
- `<feature>_metadata.yaml`: input cube, feature definition, settings, units, and software version.

## Ratio Maps

For each ratio, the package writes:

- `<feature1>_over_<feature2>_ratio.fits`
- `<feature1>_over_<feature2>_ratio_unc.fits` when both feature uncertainties are available
- `<feature1>_over_<feature2>_ratio.png`
- `<feature1>_over_<feature2>_metadata.yaml`

Intermediate native feature-map products are also written for ratio workflows.

## FITS WCS

The output FITS headers use the celestial WCS from the selected cube. When integer spatial binning is requested, the pixel scale and reference pixels are adjusted so the binned map remains spatially meaningful.

## Metadata

YAML sidecars are intended to make maps reproducible. They record:

- software version
- input cube path and spectral range
- feature definition
- continuum anchors
- map settings
- output units
- reprojection behavior for ratios
