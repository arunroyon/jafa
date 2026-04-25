# Migration Note

The original notebook workflow was specialized for the C60 18.9 um band. Its core steps map directly onto the package:

| Notebook function or value | Package equivalent |
| --- | --- |
| `load_cube(file_path)` | `jafa.load_cube(file_path)` |
| `continuum_cubic_anchor(...)` | `jafa.continuum.fit_continuum(...)` |
| `mor(..., half_window=10)` | `ContinuumSettings(morph_half_window=10)` |
| `pah_range=(18.82, 19.02)` | built-in `C60_18p9.integration_window` is now `18.75-19.15`; use a config override to reproduce the older narrow window exactly |
| `anchors=[18.55, 18.70, 19.15, 19.30]` | built-in `C60_18p9.continuum_anchor_points`; use a config override to reproduce the older anchor windows exactly |
| `integrate_pah(..., bin_2x2=True)` | `--bin-spatial 2` |
| `save_map(...)` | automatic FITS, PNG, and YAML outputs |

CLI replacement:

```bash
jafa make-map \
  --cubes my_miri_cube.fits \
  --feature C60_18p9 \
  --morph-half-window 10 \
  --bin-spatial 2 \
  --outdir results/c60_18p9
```

Python replacement:

```python
from jafa import FeatureMapper

mapper = FeatureMapper(["my_miri_cube.fits"])
mapper.make_feature_map("C60_18p9", outdir="results/c60_18p9")
```
