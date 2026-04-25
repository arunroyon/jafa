import numpy as np

from jafa.continuum import ContinuumSettings, fit_continuum


def test_fit_continuum_tracks_smooth_baseline():
    wl = np.linspace(18.5, 19.3, 160)
    baseline = 2.0 + 0.2 * (wl - 18.9)
    feature = 0.8 * np.exp(-0.5 * ((wl - 18.9) / 0.035) ** 2)
    flux = baseline + feature

    fit = fit_continuum(
        wl,
        flux,
        anchors=((18.55, 18.70), (19.15, 19.30)),
        settings=ContinuumSettings(morph_half_window=8, min_anchor_points=2),
    )

    outside = (wl < 18.75) | (wl > 19.08)
    assert np.nanmedian(np.abs(fit.continuum[outside] - baseline[outside])) < 0.08
    assert np.nanmax(fit.residual) > 0.3


def test_fit_continuum_handles_nans_and_duplicate_point_anchors():
    wl = np.linspace(1.0, 2.0, 40)
    flux = 1.0 + 0.1 * wl
    flux[5:8] = np.nan

    fit = fit_continuum(
        wl,
        flux,
        anchors=(),
        anchor_points=(1.0, 1.0, 1.5, 2.0),
        settings=ContinuumSettings(min_anchor_points=2),
    )

    assert np.all(np.diff(fit.anchor_wavelengths) > 0)
    assert np.isfinite(fit.continuum).all()


def test_fit_continuum_too_few_samples_returns_warning():
    fit = fit_continuum(
        np.array([1.0, 2.0]),
        np.array([1.0, np.nan]),
        anchors=(),
        settings=ContinuumSettings(min_finite=3),
    )
    assert "Too few finite" in fit.warnings[0]
