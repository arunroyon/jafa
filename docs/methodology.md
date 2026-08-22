# Methodology

## Cube Selection

For a feature map, the required wavelength span is the union of the feature integration window and all continuum anchor windows. A cube is eligible only if its spectral coverage fully contains that span. If multiple cubes are eligible, the package chooses the one with the smallest excess spectral coverage, with small preference bonuses for matching the feature's preferred instrument or channel.

This prevents hardcoding a MIRI channel name while still favoring the natural cube for a feature.

## Morphology Baseline

The first continuum stage uses `pybaselines.morphological.mor`. The morphology baseline is not the final continuum by itself; it is a morphology-cleaned representation used to suppress narrow emission lines before spline fitting. The default morphology half-window is 10 spectral pixels. Broader residual line structure can be tested with larger values through `--morph-half-window` or config YAML. Continuum fitting is performed over the full wavelength span of the selected cube, while integration remains restricted to the feature window.

NaNs are linearly interpolated for the baseline stage. Spectra with too few finite samples are skipped or assigned NaN outputs.

## Spline Continuum

Continuum anchors can be wavelength windows or discrete anchor points. Every sampled wavelength inside an anchor window contributes an anchor sample from the morphology-cleaned spectrum. Discrete anchor points sample the morphology-cleaned spectrum by interpolation at the requested wavelength when they lie inside the selected cube coverage; outside points are skipped with a warning. The first and last valid wavelength samples of the selected cube are also added as spline edge anchors using the morphology baseline values, which reduces endpoint spikes. If an anchor region has no sampled pixels but its midpoint lies within the cube coverage, the baseline is interpolated at the midpoint. Diagnostic figures show the full selected cube span, shade the feature and anchor windows, mark point anchors, mark the samples used for the spline fit, plot the morphology baseline, and annotate the morphology half-window value.

If too few anchor samples remain, the fitter uses an evenly spaced fallback anchor grid over the finite wavelength span and records a warning. Four or more final anchors use a natural cubic spline; fewer final anchors use linear interpolation.

## Map Integration

The residual spectrum is:

```text
residual = observed spectrum - continuum
```

The feature map value is the trapezoidal integral of the residual over the feature window. By default, negative integrated residuals are clipped to zero. This can be disabled with `--no-clip-negative`.

For cubes with a pixel area, flux-density-per-steradian data are multiplied by the pixel area so the map represents integrated flux per pixel. Integer spatial binning is supported by fitting the mean finite spectrum in each bin, then scaling each wavelength slice by the number of finite contributing pixels so partial-NaN bins do not inflate the flux. The integration uses sampled cube wavelengths inside the feature window and does not currently interpolate residual values onto exact window boundaries.

## Uncertainty

When an uncertainty cube is available, spectra are combined in quadrature for spatial bins. The integrated uncertainty uses trapezoid-equivalent wavelength weights. This is more physically aligned with the residual integration than simply summing channel uncertainties.

## Reprojection Strategy

Ratio maps are computed on a common grid. If the two feature maps are already on the same grid, no reprojection occurs.

If the maps come from different cubes, the target grid is the coarser/lower-resolution cube:

1. Larger pixel scale wins.
2. If pixel scales are similar within 5 percent, smaller spatial dimensions wins.
3. If still tied, the user priority or first cube wins.

The higher-resolution map is reprojected onto the coarser grid. The lower-resolution map is never upsampled to match the higher-resolution cube.

## Spatial Science Mask

For JWST `s3d` cubes, JAFA preserves the `DQ` and `WMAP` extensions. Voxels flagged `DO_NOT_USE` or `NON_SCIENCE`, zero-weight voxels, and non-finite science values are excluded before fitting. Each map pixel must meet the configured finite coverage in the integration window and across the full feature-plus-anchor span. WMAP is normalized by the valid median in each wavelength plane, then the spatial footprint is filtered by `min_relative_weight` and eroded by `edge_erosion_pixels`.

## Ratio Calculation

The ratio is:

```text
ratio = numerator_feature_map / denominator_feature_map
```

Feature science masks are reprojected with nearest-neighbor interpolation and intersected. Pixels are also masked if either map is non-finite or if the denominator is zero or negative. If uncertainties are present, ratio uncertainty is propagated as:

```text
sigma_ratio = ratio * sqrt((sigma_num / num)^2 + (sigma_den / den)^2)
```

Optional SNR thresholds apply to numerator and denominator when uncertainty maps are available.

## Assumptions and Limitations

- The default feature windows are conservative starting points, not universal prescriptions.
- PAH complexes can be blended; users should tune anchors for their target and science question.
- Current continuum mode is morphology plus spline. The module boundaries are designed for future continuum methods.
- Reprojected uncertainties are approximate because variance is interpolated spatially.
- The package expects JWST-like 3D FITS cubes with a usable WCS and, when present, a matching uncertainty extension.
