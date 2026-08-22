import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jafa.feature_db import FeatureDefinition
from jafa.io import CubeData
from jafa.mapping import MapSettings
from jafa.ratios import compute_ratio, make_ratio_map
from jafa.reprojection import ReprojectedMap
from jafa.utils import same_celestial_grid


def test_compute_ratio_masks_bad_denominator_and_propagates_uncertainty():
    numerator = np.array([[4.0, 2.0], [1.0, np.nan]])
    denominator = np.array([[2.0, 0.0], [-1.0, 2.0]])
    nerr = np.full_like(numerator, 0.4)
    derr = np.full_like(denominator, 0.2)

    ratio, ratio_unc = compute_ratio(
        numerator,
        denominator,
        numerator_uncertainty=nerr,
        denominator_uncertainty=derr,
    )

    assert ratio[0, 0] == 2.0
    assert np.isnan(ratio[0, 1])
    assert np.isnan(ratio[1, 0])
    assert np.isnan(ratio[1, 1])
    assert ratio_unc is not None
    assert np.isfinite(ratio_unc[0, 0])


def test_compute_ratio_applies_snr_threshold():
    numerator = np.array([[10.0, 1.0]])
    denominator = np.array([[5.0, 5.0]])
    nerr = np.array([[1.0, 1.0]])
    derr = np.array([[1.0, 1.0]])

    ratio, _ = compute_ratio(
        numerator,
        denominator,
        numerator_uncertainty=nerr,
        denominator_uncertainty=derr,
        snr_threshold=3.0,
    )

    assert ratio[0, 0] == 2.0
    assert np.isnan(ratio[0, 1])


def _celestial_header(*, crval1=10.0, crval2=20.0, cdelt=0.001):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [crval1, crval2]
    wcs.wcs.crpix = [1.0, 1.0]
    wcs.wcs.cdelt = [-cdelt, cdelt]
    return wcs.to_header()


def test_same_celestial_grid_accepts_matching_real_wcs():
    header = _celestial_header()

    assert same_celestial_grid(header, (2, 2), header.copy(), (2, 2))


def test_same_celestial_grid_rejects_different_wcs():
    header1 = _celestial_header(crval1=10.0)
    header2 = _celestial_header(crval1=10.1)

    assert not same_celestial_grid(header1, (2, 2), header2, (2, 2))


def test_same_celestial_grid_rejects_empty_headers_and_different_shapes():
    header = _celestial_header()

    assert not same_celestial_grid(fits.Header(), (2, 2), fits.Header(), (2, 2))
    assert not same_celestial_grid(header, (2, 2), header, (3, 2))


def test_ratio_map_reprojects_same_shape_different_wcs(monkeypatch):
    calls = []

    def fake_reproject_map(data, source_header, target_header, target_shape, method="interp"):
        calls.append((source_header, target_header, target_shape, method))
        return ReprojectedMap(data=np.asarray(data, dtype=float), footprint=np.ones(target_shape, dtype=float))

    def fake_reproject_mask(data, source_header, target_header, target_shape):
        return ReprojectedMap(data=np.asarray(data, dtype=bool), footprint=np.ones(target_shape, dtype=float))

    monkeypatch.setattr("jafa.ratios.reproject_map", fake_reproject_map)
    monkeypatch.setattr("jafa.ratios.reproject_mask", fake_reproject_mask)
    wl1 = np.array([1.0, 2.0, 3.0, 4.0])
    wl2 = np.array([5.0, 6.0, 7.0, 8.0])
    spec1 = np.array([0.0, 1.0, 1.0, 0.0])
    spec2 = np.array([0.0, 2.0, 2.0, 0.0])
    cube1 = CubeData(
        path=None,
        data=np.repeat(spec1[:, None, None], 2, axis=1).repeat(2, axis=2),
        wavelength_um=wl1,
        header=_celestial_header(crval1=10.0),
        flux_unit="MJy/sr",
        pixel_area_sr=1.0,
        name_hint="cube1",
    )
    cube2 = CubeData(
        path=None,
        data=np.repeat(spec2[:, None, None], 2, axis=1).repeat(2, axis=2),
        wavelength_um=wl2,
        header=_celestial_header(crval1=10.1),
        flux_unit="MJy/sr",
        pixel_area_sr=1.0,
        name_hint="cube2",
    )
    feature1 = FeatureDefinition(
        feature_name="feature1",
        central_wavelength=2.5,
        integration_window=(2.0, 3.0),
        continuum_mode="spline",
        continuum_anchor_points=(1.0, 4.0),
    )
    feature2 = FeatureDefinition(
        feature_name="feature2",
        central_wavelength=6.5,
        integration_window=(6.0, 7.0),
        continuum_mode="spline",
        continuum_anchor_points=(5.0, 8.0),
    )
    settings = MapSettings(progress=False, output_unit="native")
    settings.continuum.mode = None
    settings.continuum.min_finite = 2
    settings.continuum.min_anchor_points = 2
    settings.continuum.fallback_anchor_count = 2
    settings.continuum.include_edge_anchors = False
    settings.clip_negative_residuals = False

    result = make_ratio_map([cube1, cube2], feature1, feature2, settings=settings, write_outputs=False)

    assert calls
    assert result.metadata["ratio"]["reprojection"] == "feature2_reprojected_to_feature1_grid"
    assert np.allclose(result.ratio, 0.5)
