Outputs
=======

Feature maps
------------

For each feature, the writer creates:

- ``<feature>_map.fits``
- ``<feature>_unc.fits`` when uncertainty is available
- ``<feature>_continuum.fits`` when enabled
- ``<feature>_mask.fits`` (1 = retained science pixel)
- ``<feature>_coverage.fits`` (minimum feature/anchor-span coverage fraction)
- ``<feature>_relative_weight.fits`` when WMAP is available
- ``<feature>_diagnostic.png``
- ``<feature>_metadata.yaml``

The FITS products use the selected cube's celestial WCS. Metadata records the
selected cube, feature definition, ignored anchor points, continuum settings,
software version, and output unit.

When a product is made from an adjacent-cube stitch, the metadata field
``cube.stitched_from`` lists the source cubes and the WCS/header correspond to
the chosen common target grid. ``cube.stitch_target`` records the source cube
used as the spatial target grid, while ``cube.stitch_overlap_strategy`` and
``cube.stitch_overlap_ranges_um`` record how overlapping wavelength coverage
was handled. ``cube.stitch_common_footprint_fraction`` records how much of the
target grid remains after requiring spatial coverage from every stitched cube.

By default, feature maps, continuum maps, and uncertainty maps are written as
pixel-integrated cgs fluxes with ``BUNIT = erg s-1 cm-2 pixel-1``. For
``MJy/sr`` cubes, this is computed as an exact trapezoidal
``integral F_nu dnu`` over the feature window multiplied by the spatial pixel
solid angle. Set ``output_unit: native`` or pass ``--output-unit native`` to
write wavelength-integrated products such as ``MJy um`` instead.

Ratio maps
----------

For each ratio, the writer creates:

- ``<feature1>_over_<feature2>_ratio.fits``
- ``<feature1>_over_<feature2>_ratio_unc.fits`` when possible
- ``<feature1>_over_<feature2>_mask.fits``
- ``<feature1>_over_<feature2>_ratio.png``
- ``<feature1>_over_<feature2>_metadata.yaml``

Intermediate feature products are also written in the ratio output directory.
