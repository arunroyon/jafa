![JAFA banner](docs/assets/jafa_banner.png)

<p align="center"><sub>Banner image credit: ChatGPT.</sub></p>

# JAFA — JWST Aromatic Feature Analyzer

JAFA — JWST Aromatic Feature Analyzer creates continuum-subtracted JWST IFU spectral feature maps from MIRI MRS and future NIRSpec IFU cubes. It is designed for PAH and fullerene bands, with feature windows and continuum anchors defined in editable YAML. The Python import package and primary CLI command are both `jafa`.

The package grew out of a working 18.9 um fullerene mapping notebook. The original scientific behavior is preserved in generalized form: morphology-based smoothing/baseline estimation with `pybaselines`, spline continuum fitting through anchor regions, residual integration, FITS output with WCS, diagnostic plots, and ratio maps.

## Example Outputs

TC1 C60 18.9 um continuum diagnostic:

![TC1 C60 18.9 um diagnostic](examples/gallery/tc1_c60_18p9_diagnostic.png)

NGC 7023 PAH 11.0/11.2 ratio map:

![NGC 7023 PAH 11.0 over 11.2 ratio](examples/gallery/ngc7023_pah_11p0_over_11p2_ratio.png)

Scientific validation and citation guidance will be described in Arun et al.,
in preparation. For questions or collaboration, contact Arun Roy at
arunroyon@gmail.com.

## Installation

Core install:

```bash
python -m pip install jafa
```

The core import only requires `numpy`, `scipy`, `astropy`, and `PyYAML`. JWST cube loading, morphology baselines, reprojection, and plotting are optional extras that load only when the relevant functions are called.

For the full JWST science workflow:

```bash
python -m pip install "jafa[jwst,plot]"
```

For a clean conda development environment:

```bash
conda env create -f environment.yml
conda activate jafa
python -m pip install -e .
```

For local development with tests and docs:

```bash
python -m pip install -e ".[jwst,plot,test,docs,dev]"
```

Legacy compatibility: the old `jwst_feature_mapper` import path and
`jwst-feature-map` CLI command are retained as transitional aliases. New
scripts should use `jafa`.

`pybaselines` can trigger numba cache issues in some notebook or packaged environments. The package sets `NUMBA_DISABLE_JIT=1` immediately before lazily importing the morphology baseline routine, matching the original notebook behavior without making `pybaselines` an import-time dependency.

## Command Line

Single feature:

```bash
jafa make-map \
  --cubes cube1.fits cube2.fits \
  --feature PAH_11p2 \
  --outdir results/
```

Two-feature ratio:

```bash
jafa make-ratio \
  --cubes cube1.fits cube2.fits cube3.fits \
  --feature1 PAH_11p2 \
  --feature2 PAH_6p2 \
  --outdir results/
```

Custom wavelength mode:

```bash
jafa make-map \
  --cubes cube1.fits cube2.fits \
  --wavelength 18.9 \
  --feature-window 18.82 19.02 \
  --anchors 18.55 18.70 19.15 19.30 \
  --outdir results/
```

Useful options:

```bash
--config examples/example_config.yaml
--feature-db custom_features.yaml
--snr-threshold 3
--morph-half-window 10
--bin-spatial 2
--output-unit cgs        # default: exact frequency-integrated erg s-1 cm-2 pixel-1
--output-unit native     # legacy wavelength-integrated cube units, e.g. MJy um
--no-stitch-cubes        # require one cube to cover the whole feature
--stitch-gap-tolerance 0.01
--stitch-overlap split   # default: use one cube per side of the overlap midpoint
--allow-noncommon-stitch-footprint
--no-clip-negative
--quiet
--verbose
```

## Python API

```python
from jafa import FeatureMapper

mapper = FeatureMapper(["cube1.fits", "cube2.fits"])
mapper.make_feature_map("C60_18p9", outdir="results/c60")
mapper.make_ratio_map("PAH_11p2", "PAH_6p2", outdir="results/ratio")
```

Lower-level functions are also importable:

```python
from jafa import load_cube, load_feature_database, make_feature_map

cubes = [load_cube("cube_ch3.fits")]
features = load_feature_database()
make_feature_map(cubes, features["C60_18p9"], outdir="results/c60")
```

## Feature Database

Built-in features live in `src/jafa/data/features.yaml`. Each entry has this schema:

```yaml
C60_18p9:
  feature_name: C60_18p9
  central_wavelength: 18.90
  integration_window: [18.75, 19.15]
  continuum_mode: morph_spline
  continuum_anchors: []
  continuum_anchor_points: [17.8, 18.0, 18.2, 18.4, 18.6, 19.2, 19.5, 19.8, 20.2, 20.5]
  notes: Fullerene C60 band.
  preferred_instrument: MIRI MRS
  preferred_channel: Channel 4
```

Included defaults:

`PAH_3p3`, `PAH_6p2`, `PAH_7p7`, `PAH_8p6`, `PAH_11p0`, `PAH_11p2`, `PAH_12p7`, `PAH_16p4`, `C60_7p0`, `C60_8p5`, `C60_17p4`, `C60_18p9`.

Current fullerene defaults include user-tuned TC1 anchors for `C60_7p0`, `C60_8p5`, `C60_17p4`, and `C60_18p9`. Treat them as editable science defaults, not immutable line prescriptions.

Pass `--feature-db custom_features.yaml` to merge or override entries. You can also put custom feature entries under a top-level `features:` block in a pipeline config YAML and pass it with `--config`.

## Workflow

For a requested feature, the software:

1. Loads one or more FITS cubes.
2. Reads the wavelength axis, flux unit, WCS, uncertainty extension when available, and pixel area.
3. Selects the cube whose spectral coverage contains the feature integration window, continuum anchor windows, and discrete anchor points.
4. Fits each spaxel continuum over the full selected cube wavelength span using morphology baseline estimation plus spline anchor fitting.
5. Integrates the continuum-subtracted residual over the feature window.
6. Writes FITS maps, diagnostic PNGs, and metadata YAML.

If no single cube covers the requested feature plus anchors, the mapper can
automatically stitch the narrowest adjacent sequence of cubes. It reprojects
spectral planes onto the coarser/lower-resolution spatial grid, concatenates
and sorts the wavelength axis, splits overlap regions at their midpoint by
default, and then runs the same continuum and integration workflow on the
stitched in-memory cube. This is the intended path for broad PAH complexes
such as the 7.7+8.6 um complex when the feature and continuum anchors span
neighboring MIRI MRS sub-bands. The stitch is recorded in each metadata YAML
under `cube.stitched_from`.

The default PAH 7.7 and 8.6 definitions share a continuum basis over the
stitched 7.7+8.6 complex. Anchor points are `6.66, 6.8, 7.0, 7.15, 9.0, 9.2`
plus the automatic edge anchors. `PAH_7p7` integrates `7.15-8.25 um`;
`PAH_8p6` integrates `8.25-9.0 um`.

Pixel-grid stitching is not full PSF matching. For publication work across
MIRI channels, the sharper/shorter-wavelength cube should ideally also be
convolved to the broader PSF before, or as part of, reprojection. The current
pipeline conservatively matches the spatial pixel grid and avoids duplicated
overlap samples. By default it also masks the stitched cube to the common
spatial footprint of all input channels, so pixels covered by only one channel
become NaN at every stitched wavelength. PSF-kernel matching remains an
explicit validation step.

The default morphology half-window is `10` spectral pixels, matching the original notebook intent: morphology first suppresses narrow emission lines, then the spline continuum is fit to the morphology-cleaned spectrum at anchor samples. Built-in features therefore use `continuum_mode: morph_spline`; use `--continuum-mode spline` only when you intentionally want to skip the morphology stage. The first and last valid wavelength samples of the selected cube are also added as edge anchors to suppress spline endpoint spikes. Increase `--morph-half-window` only when line residuals are broader than the default morphology scale. Diagnostic plots show the full fitted cube span, feature window, anchor windows, anchor samples, morphology baseline, spline continuum, residual, and the morphology half-window value.

By default, feature maps are written in cgs integrated flux units,
`erg s-1 cm-2 pixel-1`. For JWST-style `MJy/sr` cubes, the mapper converts the
surface brightness residual to `erg s-1 cm-2 Hz-1 sr-1`, integrates exactly over
the corresponding frequency grid, and multiplies by the output pixel solid
angle. This implements the cgs conversion that the older plotting notebook
applied only for visualization. Use `--output-unit native` to reproduce the
older wavelength-integrated `MJy um` products.

Integration uses cube samples that fall inside the feature window; it does not
currently interpolate residual values onto exact window boundaries.

## Ratio Maps

If both features are selected from the same cube or already share a grid, the ratio is computed directly.

If the features come from different cubes:

1. Each feature map is created in its native cube.
2. The common grid is chosen from the coarser/lower-resolution cube.
3. Coarser means larger pixel scale; if similar within 5 percent, smaller spatial dimensions; if tied, the first selected cube.
4. The higher-resolution map and uncertainty are reprojected to the coarser grid.
5. The ratio is computed with invalid, zero, and negative denominator pixels masked.

The package deliberately avoids upsampling the lower-resolution map to the higher-resolution cube.

## Outputs

Feature map products:

- `<feature>_map.fits`
- `<feature>_unc.fits` when uncertainty is available
- `<feature>_continuum.fits`
- `<feature>_diagnostic.png`
- `<feature>_metadata.yaml`

The default FITS `BUNIT` for feature maps is `erg s-1 cm-2 pixel-1`. Metadata
records whether the map used exact cgs frequency integration or legacy native
wavelength integration.

Ratio products:

- `<feature1>_over_<feature2>_ratio.fits`
- `<feature1>_over_<feature2>_ratio_unc.fits` when possible
- `<feature1>_over_<feature2>_ratio.png`
- `<feature1>_over_<feature2>_metadata.yaml`

## Troubleshooting

- `No cube covers ...`: provide the cube channel containing both the feature and all continuum anchors, or edit the anchors.
- Stitched broad-band products: inspect the diagnostic spectrum near the channel overlap; adjacent MIRI sub-bands can have small calibration steps.
- Many NaNs in the output: check input masks, wavelength coverage, and SNR threshold.
- Flat or unstable continuum: increase anchor coverage, adjust `morph_half_window`, or inspect the diagnostic PNG.
- Cgs conversion fails: confirm the input cube `BUNIT` is a flux-density surface-brightness unit such as `MJy/sr` and that `PIXAR_SR` or a valid celestial WCS is available.
- Ratio map mostly blank: denominator values may be zero, negative, below SNR threshold, or outside the reprojected footprint.

## Development

Common checks:

```bash
pytest
ruff check src tests
tox -e py312
sphinx-build -b html docs docs/_build/html
python -m build
```

The package is structured in an Astropy-friendly style: src layout, package data for the feature database, FITS/WCS handling through Astropy, Quantity-aware wavelength handling, and optional heavy dependencies behind lazy imports.

## Migration From the Original 18.9 um Notebook

Old notebook concept:

```python
datacube, errorcube, wavelengths, header = load_cube(file_path)
pah_map, pah_unc = integrate_pah(
    datacube,
    errorcube,
    wavelengths,
    pah_range=(18.82, 19.02),
    anchors=[18.55, 18.70, 19.15, 19.30],
    pixel_area_sr=header["PIXAR_SR"],
    morph_half_window=10,
    bin_2x2=True,
)
save_map(pah_map, pah_unc, header, pah_range, outfile, bin_2x2=True)
```

The old maker notebook wrote wavelength-integrated maps in `MJy um`, and the
old plotting notebook multiplied those maps by approximately
`2.99792458e-3 / lambda_um**2` for cgs visualization. The package now performs
the cgs conversion during map creation with the exact frequency grid, so the
saved FITS maps are ready for publication by default.

New package equivalent:

```bash
jafa make-map \
  --cubes file_path.fits \
  --feature C60_18p9 \
  --morph-half-window 10 \
  --bin-spatial 2 \
  --outdir results/c60_18p9
```

Or in Python:

```python
from jafa import FeatureMapper

mapper = FeatureMapper(["file_path.fits"])
mapper.make_feature_map("C60_18p9", outdir="results/c60_18p9")
```
