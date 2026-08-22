Methodology
===========

Feature definitions
-------------------

Each feature is described in YAML with a central wavelength, an integration
window, continuum anchor windows, optional discrete anchor points, and optional
instrument/channel hints. The integration window defines where residual flux is
integrated. Continuum anchors define where the spline continuum is constrained.

Cube selection
--------------

For a single feature, a cube is eligible when its wavelength coverage contains
the feature integration window, anchor windows, and discrete anchor points. The
first and last valid cube wavelengths are added as edge anchors when enabled.
If multiple cubes are eligible, the mapper selects the narrowest coverage with
small preference bonuses for matching feature metadata.

Adjacent-cube stitching
-----------------------

When no single cube covers the full required wavelength span, the mapper can
construct a temporary stitched cube from the narrowest adjacent cube sequence.
This is useful for broad PAH complexes, especially the 7.7+8.6 um complex,
whose feature windows and continuum anchors can cross MIRI MRS sub-band
boundaries.

The stitching procedure is:

#. Select the narrowest adjacent sequence whose combined coverage contains the
   required feature and anchor span without a spectral gap.
#. Choose the common spatial grid using the same conservative rule as ratio
   maps: the coarser/lower-resolution cube is the target.
#. Reproject each spectral plane of the other cube onto that target celestial
   grid.
#. Convert compatible flux units if needed, concatenate spectra along the
   wavelength axis, sort by wavelength, and split overlap regions at their
   midpoint by default so the same spectral interval is not represented twice.
#. Apply the common spatial footprint of all stitched cubes, setting pixels
   covered by only one channel to NaN at every wavelength.
#. Fit the continuum and integrate the feature from this stitched in-memory cube.

The stitch metadata records the source cubes. Overlap regions are not
empirically rescaled by default; this avoids fitting a calibration correction
inside a science feature. Diagnostics should therefore be inspected near the
channel overlap for possible sub-band calibration steps.

For the default PAH 7.7 and 8.6 definitions, both features use the same
stitched continuum basis over the 7.7+8.6 complex. The continuum anchor points
are 6.66, 6.8, 7.0, 7.15, 9.0, and 9.2 um, plus automatic cube edge anchors.
The 7.7 map integrates the residual from 7.15 to 8.25 um; the 8.6 map
integrates 8.25 to 9.0 um.

Spatial reprojection to a common pixel grid is necessary but not identical to
PSF matching. For publication measurements across MIRI channels, the
shorter-wavelength/sharper cube should ideally be convolved to the broader PSF
of the longer-wavelength cube before final science interpretation. The current
stitching machinery records the common grid and overlap policy, but it does not
yet apply a wavelength-dependent PSF kernel.

Morphology baseline
-------------------

The default continuum workflow is morphology plus spline. First,
``pybaselines.morphological.mor`` estimates a morphology-cleaned baseline. This
stage suppresses narrow emission lines before the spline continuum is fit. The
default morphology half-window is 10 spectral pixels and is configurable.

Spline continuum
----------------

Anchor windows contribute every sampled wavelength inside the window. Discrete
anchor points sample the morphology baseline by interpolation. Duplicate anchor
wavelengths are averaged before spline fitting. If too few anchors remain, the
fitter uses an evenly spaced fallback anchor grid over the finite wavelength
span and records a warning. Four or more final anchors use a natural cubic
spline; fewer final anchors use linear interpolation.

Map integration
---------------

The map product first forms the continuum-subtracted residual:

.. code-block:: text

   observed spectrum - spline continuum

For the default ``output_unit="cgs"`` mode, JWST-style surface brightness
spectra such as ``MJy/sr`` are converted to
``erg s-1 cm-2 Hz-1 sr-1`` and integrated over frequency:

.. code-block:: text

   integral F_nu dnu * pixel_area_sr

The frequency grid is computed from each sampled wavelength, so this is an
exact trapezoidal frequency integration rather than a central-wavelength scale
factor. The resulting FITS map unit is ``erg s-1 cm-2 pixel-1``.
The integration uses sampled cube wavelengths inside the feature window and
does not currently interpolate residual values onto exact window boundaries.
For spatially binned products, each wavelength slice is scaled by the number
of finite pixels that actually contributed to that slice, so partial-NaN bins
do not inflate the integrated flux.

The optional ``output_unit="native"`` mode integrates the residual over
wavelength in the cube's native flux-density units, for example
``MJy um`` after pixel-area scaling. Negative integrated values are clipped to
zero by default and can be preserved with ``--no-clip-negative``.

Uncertainty
-----------

When an uncertainty cube is available, spatial bins combine uncertainty in
quadrature. Integrated uncertainty is computed with trapezoid-equivalent
weights in the same coordinate used for the science map: frequency weights for
``cgs`` output and wavelength weights for ``native`` output.

Ratio maps
----------

Ratio maps are built from maps, not full cubes. If two feature maps already
share a grid, the ratio is computed directly. If they come from different
cubes, each map is made in its native cube, the lower-resolution/coarser cube is
chosen as the target grid, and the higher-resolution map is reprojected to that
grid. The lower-resolution map is not upsampled. Because feature maps are in
cgs integrated flux units by default, cross-wavelength ratios compare physical
integrated fluxes rather than ``MJy um`` surface-brightness integrals.

Masking
-------

For JWST ``s3d`` cubes, voxels carrying ``DO_NOT_USE`` or ``NON_SCIENCE`` DQ
bits, zero WMAP weight, or non-finite science values are rejected before the
continuum fit. Each output pixel must meet the configured finite-sample
fractions in both the feature window and the full feature-plus-anchor span.
WMAP values are normalized by the valid median in each wavelength plane; the
median relative weight across the required span must meet
``min_relative_weight``. The resulting footprint is eroded by
``edge_erosion_pixels``.

Ratio masks are the intersection of the two feature science masks. Masks are
reprojected with nearest-neighbor interpolation, while science and variance
continue to use bilinear interpolation. Ratio pixels are also rejected when
the denominator is zero or negative. Optional SNR cuts are applied separately
when uncertainty maps are present.

Limitations
-----------

Default windows and anchors are starting points for scientific inspection, not
universal prescriptions. PAH complexes and fullerene bands can be blended with
ionic lines and plateaus, so diagnostics should be inspected for each target.
Reprojected uncertainties are approximate because variance is interpolated
spatially.
