import astropy.units as u
import numpy as np
import pytest
from astropy.constants import c
from astropy.io import fits

from jafa.continuum import ContinuumFit
from jafa.feature_db import FeatureDefinition, load_feature_database
from jafa.io import CubeData
from jafa.mapping import MapSettings, _integrate_flux, make_feature_map


def test_make_feature_map_on_synthetic_cube():
    wl = np.linspace(17.8, 21.0, 180)
    baseline = 1.0 + 0.05 * (wl - 18.9)
    feature_profile = 0.6 * np.exp(-0.5 * ((wl - 18.9) / 0.035) ** 2)
    spectrum = baseline + feature_profile
    data = np.repeat(spectrum[:, None, None], 4, axis=1)
    data = np.repeat(data, 4, axis=2)
    unc = np.full_like(data, 0.02)
    cube = CubeData(
        path=None,
        data=data,
        uncertainty=unc,
        wavelength_um=wl * u.micron,
        header=fits.Header(),
        flux_unit="MJy/sr",
        pixel_area_sr=1.0,
    )
    feature = load_feature_database()["C60_18p9"]
    settings = MapSettings(progress=False)

    result = make_feature_map([cube], feature, settings=settings, write_outputs=False)

    assert result.feature_map.shape == (4, 4)
    assert np.nanmedian(result.feature_map) > 0
    assert result.uncertainty is not None
    assert np.all(np.isfinite(result.uncertainty))
    assert result.unit == "erg s-1 cm-2 pixel-1"
    assert result.metadata["settings"]["output_unit_mode"] == "cgs"
    assert result.metadata["settings"]["spectral_integration"] == "frequency_integral_exact"


def test_native_output_unit_preserves_wavelength_integral_behavior():
    wl = np.linspace(17.8, 21.0, 180)
    baseline = np.ones_like(wl)
    feature_profile = 0.6 * np.exp(-0.5 * ((wl - 18.9) / 0.035) ** 2)
    data = np.repeat((baseline + feature_profile)[:, None, None], 2, axis=1)
    data = np.repeat(data, 2, axis=2)
    cube = CubeData(
        path=None,
        data=data,
        wavelength_um=wl * u.micron,
        header=fits.Header(),
        flux_unit="MJy/sr",
        pixel_area_sr=2.0,
    )
    feature = load_feature_database()["C60_18p9"]
    settings = MapSettings(progress=False, output_unit="native")

    result = make_feature_map([cube], feature, settings=settings, write_outputs=False)

    assert result.unit == "MJy um"
    assert result.metadata["settings"]["output_unit_mode"] == "native"
    assert result.metadata["settings"]["spectral_integration"] == "wavelength_integral_native"


def test_cgs_integration_uses_exact_frequency_grid_and_pixel_area():
    wavelengths = np.array([18.75, 19.15])
    values = np.ones_like(wavelengths)
    pixel_area_sr = 2.5e-13
    cube = CubeData(
        path=None,
        data=np.ones((2, 1, 1)),
        wavelength_um=wavelengths * u.micron,
        header=fits.Header(),
        flux_unit="MJy/sr",
        pixel_area_sr=pixel_area_sr,
    )

    value = _integrate_flux(wavelengths, values, cube, npix=1, output_unit="cgs")

    nu = (c / (wavelengths * u.micron)).to_value(u.Hz)
    expected = 1.0e-17 * abs(nu[1] - nu[0]) * pixel_area_sr
    assert value == pytest.approx(expected)


def _binned_test_feature() -> FeatureDefinition:
    return FeatureDefinition(
        feature_name="bin_test",
        central_wavelength=2.5,
        integration_window=(2.0, 3.0),
        continuum_mode="spline",
        continuum_anchor_points=(1.0, 4.0),
    )


def _binned_test_settings() -> MapSettings:
    settings = MapSettings(progress=False, output_unit="native", bin_spatial=2)
    settings.continuum.mode = None
    settings.continuum.min_finite = 2
    settings.continuum.min_anchor_points = 2
    settings.continuum.fallback_anchor_count = 2
    settings.continuum.include_edge_anchors = False
    settings.clip_negative_residuals = False
    return settings


def _binned_test_cube(data: np.ndarray) -> CubeData:
    return CubeData(
        path=None,
        data=data,
        wavelength_um=np.array([1.0, 2.0, 3.0, 4.0]) * u.micron,
        header=fits.Header(),
        flux_unit="MJy/sr",
        pixel_area_sr=1.0,
    )


def test_bin_spatial_all_finite_scales_by_all_contributing_pixels():
    spectrum = np.array([0.0, 1.0, 1.0, 0.0])
    data = np.repeat(spectrum[:, None, None], 2, axis=1).repeat(2, axis=2)

    result = make_feature_map([_binned_test_cube(data)], _binned_test_feature(), settings=_binned_test_settings(), write_outputs=False)

    assert result.feature_map.shape == (1, 1)
    assert result.feature_map[0, 0] == pytest.approx(4.0)


def test_bin_spatial_partial_nans_do_not_inflate_flux():
    spectrum = np.array([0.0, 1.0, 1.0, 0.0])
    data = np.repeat(spectrum[:, None, None], 2, axis=1).repeat(2, axis=2)
    data[:, 0, 0] = np.nan

    result = make_feature_map([_binned_test_cube(data)], _binned_test_feature(), settings=_binned_test_settings(), write_outputs=False)

    assert result.feature_map[0, 0] == pytest.approx(3.0)


def test_bin_spatial_variable_finite_counts_use_per_wavelength_area():
    spectrum = np.array([0.0, 1.0, 1.0, 0.0])
    data = np.repeat(spectrum[:, None, None], 2, axis=1).repeat(2, axis=2)
    data[2, 0, :] = np.nan

    result = make_feature_map([_binned_test_cube(data)], _binned_test_feature(), settings=_binned_test_settings(), write_outputs=False)

    assert result.feature_map[0, 0] == pytest.approx(3.0)


def test_feature_window_integration_uses_sampled_points_without_endpoint_interpolation():
    cube = _binned_test_cube(np.ones((4, 2, 2)))
    value = _integrate_flux(
        np.array([2.0, 3.0]),
        np.array([1.0, 1.0]),
        cube,
        npix=4,
        pixel_counts=np.array([4.0, 4.0]),
        output_unit="native",
    )

    assert value == pytest.approx(4.0)


def test_feature_continuum_mode_controls_runtime_fitting(monkeypatch):
    calls = []

    def fake_fit_continuum(wavelengths, flux, anchors, anchor_points, settings):
        calls.append(settings.mode)
        zeros = np.zeros_like(np.asarray(flux, dtype=float))
        return ContinuumFit(
            continuum=zeros,
            baseline=np.asarray(flux, dtype=float),
            residual=np.asarray(flux, dtype=float),
            anchor_wavelengths=np.asarray(anchor_points, dtype=float),
            anchor_fluxes=np.zeros(len(anchor_points), dtype=float),
        )

    monkeypatch.setattr("jafa.mapping.fit_continuum", fake_fit_continuum)
    cube = _binned_test_cube(np.ones((4, 2, 2)))
    settings = _binned_test_settings()
    settings.continuum.mode = None

    result = make_feature_map([cube], _binned_test_feature(), settings=settings, write_outputs=False)

    assert calls
    assert set(calls) == {"spline"}
    assert result.metadata["settings"]["continuum_mode"] == "spline"
