import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from jafa.feature_db import FeatureDefinition
from jafa.io import CubeData
from jafa.mapping import MapSettings, make_feature_map
from jafa.ratios import make_ratio_map
from jafa.stitching import _same_grid, select_adjacent_cube_sequence, stitch_adjacent_cubes_for_feature


def _cube(name, wavelengths, value, shape=(3, 3)):
    data = np.full((len(wavelengths), *shape), float(value))
    return CubeData(
        path=None,
        data=data,
        wavelength_um=np.asarray(wavelengths, dtype=float),
        header=_celestial_header(),
        flux_unit="MJy/sr",
        pixel_area_sr=1.0,
        name_hint=name,
    )


def _celestial_header():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.crpix = [1.0, 1.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    return wcs.to_header()


def test_stitch_adjacent_cubes_for_feature_covers_anchor_points():
    blue = _cube("blue", np.linspace(1.0, 2.0, 8), 1.0)
    red = _cube("red", np.linspace(1.8, 3.0, 8), 2.0)
    feature = FeatureDefinition(
        feature_name="joined",
        central_wavelength=2.0,
        integration_window=(1.5, 2.5),
        continuum_anchor_points=(1.0, 3.0),
    )

    stitched = stitch_adjacent_cubes_for_feature([blue, red], feature)

    assert stitched.spectral_range == (1.0, 3.0)
    assert stitched.spatial_shape == blue.spatial_shape
    assert stitched.stitched_from == ("blue", "red")
    assert stitched.stitch_target == "blue"
    assert stitched.stitch_overlap_strategy == "split"
    assert stitched.stitch_overlap_ranges_um == ((1.8, 2.0),)
    assert stitched.stitch_common_footprint_fraction == 1.0
    assert np.all(np.diff(stitched.wavelength_um) > 0)
    overlap = (stitched.wavelength_um >= 1.8) & (stitched.wavelength_um <= 2.0)
    assert np.all(stitched.data[overlap & (stitched.wavelength_um <= 1.9)] == 1.0)
    assert np.all(stitched.data[overlap & (stitched.wavelength_um > 1.9)] == 2.0)


def test_stitching_same_grid_rejects_placeholder_headers():
    assert not _same_grid(fits.Header(), fits.Header(), (2, 2), (2, 2))


def test_make_feature_map_falls_back_to_stitched_cube():
    wl_blue = np.linspace(1.0, 2.0, 24)
    wl_red = np.linspace(1.8, 3.0, 28)
    blue_spectrum = 1.0 + 0.2 * np.exp(-0.5 * ((wl_blue - 1.9) / 0.08) ** 2)
    red_spectrum = 1.0 + 0.2 * np.exp(-0.5 * ((wl_red - 1.9) / 0.08) ** 2)
    blue = CubeData(
        path=None,
        data=np.repeat(blue_spectrum[:, None, None], 4, axis=1).repeat(4, axis=2),
        wavelength_um=wl_blue,
        header=_celestial_header(),
        flux_unit="MJy/sr",
        name_hint="blue",
    )
    red = CubeData(
        path=None,
        data=np.repeat(red_spectrum[:, None, None], 4, axis=1).repeat(4, axis=2),
        wavelength_um=wl_red,
        header=_celestial_header(),
        flux_unit="MJy/sr",
        name_hint="red",
    )
    feature = FeatureDefinition(
        feature_name="joined_map",
        central_wavelength=1.9,
        integration_window=(1.75, 2.1),
        continuum_anchor_points=(1.0, 1.2, 2.8, 3.0),
    )
    settings = MapSettings(progress=False, output_unit="native")
    settings.continuum.mode = "spline"

    result = make_feature_map([blue, red], feature, settings=settings, write_outputs=False)

    assert result.cube.stitched_from == ("blue", "red")
    assert np.nanmedian(result.feature_map) > 0


def test_stitching_masks_to_common_spatial_footprint():
    blue = _cube("blue", np.linspace(1.0, 2.0, 8), 1.0)
    red = _cube("red", np.linspace(1.8, 3.0, 8), 2.0)
    red.data[:, :, -1] = np.nan
    feature = FeatureDefinition(
        feature_name="joined",
        central_wavelength=2.0,
        integration_window=(1.5, 2.5),
        continuum_anchor_points=(1.0, 3.0),
    )

    stitched = stitch_adjacent_cubes_for_feature([blue, red], feature)

    assert stitched.stitch_common_footprint_fraction == 2 / 3
    assert np.all(np.isnan(stitched.data[:, :, -1]))
    assert np.all(np.isfinite(stitched.data[:, :, :-1]))


def test_stitching_preserves_weight_maps_for_feature_masking():
    blue = _cube("blue", np.linspace(1.0, 2.0, 8), 1.0)
    red = _cube("red", np.linspace(1.8, 3.0, 8), 2.0)
    blue.weight = np.ones_like(blue.data)
    red.weight = np.full_like(red.data, 2.0)
    feature = FeatureDefinition(
        feature_name="joined",
        central_wavelength=2.0,
        integration_window=(1.5, 2.5),
        continuum_anchor_points=(1.0, 3.0),
    )

    stitched = stitch_adjacent_cubes_for_feature([blue, red], feature)

    assert stitched.weight is not None
    assert stitched.weight.shape == stitched.data.shape
    assert np.all(np.isfinite(stitched.weight))


def test_select_adjacent_cube_sequence_can_use_three_cubes():
    blue = _cube("blue", np.linspace(6.5, 7.65, 8), 1.0)
    middle = _cube("middle", np.linspace(7.51, 8.77, 8), 2.0)
    red = _cube("red", np.linspace(8.67, 9.3, 8), 3.0)

    selected = select_adjacent_cube_sequence([blue, middle, red], (6.66, 9.2))

    assert selected == (blue, middle, red)


def test_ratio_map_uses_stitched_complex_for_both_features():
    wl_blue = np.linspace(6.5, 7.65, 24)
    wl_middle = np.linspace(7.51, 8.77, 28)
    wl_red = np.linspace(8.67, 9.3, 18)
    anchors = (6.66, 6.8, 7.0, 7.15, 9.0, 9.2)
    feature1 = FeatureDefinition(
        feature_name="PAH_7p7_test",
        central_wavelength=7.7,
        integration_window=(7.15, 8.25),
        continuum_anchor_points=anchors,
    )
    feature2 = FeatureDefinition(
        feature_name="PAH_8p6_test",
        central_wavelength=8.6,
        integration_window=(8.25, 9.0),
        continuum_anchor_points=anchors,
    )

    def spectrum(wavelength):
        return (
            1.0
            + 0.3 * np.exp(-0.5 * ((wavelength - 7.7) / 0.15) ** 2)
            + 0.2 * np.exp(-0.5 * ((wavelength - 8.6) / 0.12) ** 2)
        )

    cubes = [
        CubeData(
            path=None,
            data=np.repeat(spectrum(wl_blue)[:, None, None], 3, axis=1).repeat(3, axis=2),
            wavelength_um=wl_blue,
            header=_celestial_header(),
            flux_unit="MJy/sr",
            name_hint="blue",
        ),
        CubeData(
            path=None,
            data=np.repeat(spectrum(wl_middle)[:, None, None], 3, axis=1).repeat(3, axis=2),
            wavelength_um=wl_middle,
            header=_celestial_header(),
            flux_unit="MJy/sr",
            name_hint="middle",
        ),
        CubeData(
            path=None,
            data=np.repeat(spectrum(wl_red)[:, None, None], 3, axis=1).repeat(3, axis=2),
            wavelength_um=wl_red,
            header=_celestial_header(),
            flux_unit="MJy/sr",
            name_hint="red",
        ),
    ]
    settings = MapSettings(progress=False, output_unit="native")
    settings.continuum.mode = "spline"

    result = make_ratio_map(cubes, feature1, feature2, settings=settings, write_outputs=False)

    assert result.numerator.cube.stitched_from == ("blue", "middle", "red")
    assert result.denominator.cube.stitched_from == ("blue", "middle", "red")
    assert np.nanmedian(result.numerator.feature_map) > 0
    assert np.nanmedian(result.denominator.feature_map) > 0
