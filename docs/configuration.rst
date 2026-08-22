Configuration
=============

Feature database
----------------

Built-in feature definitions live in the packaged data file
``src/jafa/data/features.yaml``. A custom YAML file can override or
extend these definitions:

.. code-block:: yaml

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

Pass the file with:

.. code-block:: bash

   jafa make-map --feature-db custom_features.yaml ...

Pipeline settings
-----------------

Config files can include settings:

.. code-block:: yaml

   settings:
     continuum_mode: morph_spline
     morph_half_window: 10
     min_finite: 10
     min_anchor_points: 4
     fallback_anchor_count: 4
     include_edge_anchors: true
     snr_threshold: null
     clip_negative_residuals: true
     bin_spatial: 1
     min_feature_coverage: 0.95
     min_required_coverage: 0.95
     min_relative_weight: 0.25
     edge_erosion_pixels: 1
     output_unit: cgs
     allow_cube_stitching: true
     stitch_spectral_gap_tolerance_um: 0.0
     stitch_overlap_strategy: split
     stitch_min_spatial_coverage: 0.95
     stitch_require_common_spatial_footprint: true

They can also include feature overrides under a top-level ``features`` block.

``output_unit`` may be ``cgs`` or ``native``. The default ``cgs`` mode writes
pixel-integrated ``erg s-1 cm-2 pixel-1`` maps using exact frequency
integration. ``native`` preserves the older wavelength-integrated cube-unit
behavior.

``allow_cube_stitching`` controls the automatic adjacent-cube fallback used
when a feature and its anchors span neighboring MIRI sub-bands. The default
gap tolerance requires overlapping or touching spectral coverage. The default
overlap strategy splits each spectral overlap at its midpoint; use ``keep``
only for debugging or specialized workflows.

``stitch_require_common_spatial_footprint`` masks stitched products to pixels
covered by all source cubes. This avoids fitting a continuum from a spectrum
that only exists in one channel.

Spatial edge masking uses the feature-window and full anchor-span coverage
fractions together with the JWST ``WMAP`` extension when available. WMAP is
normalized by the valid median in each wavelength plane before applying
``min_relative_weight``. The resulting footprint is eroded by
``edge_erosion_pixels`` and written alongside each map as mask, coverage, and
relative-weight FITS products. Thresholds should be varied in publication
sensitivity checks because WMAP normalization depends on the observing and
cube-build configuration.
