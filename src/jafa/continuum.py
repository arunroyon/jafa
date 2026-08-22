"""Continuum estimation using morphology smoothing plus spline anchors."""

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import astropy.units as u
import numpy as np
from scipy.interpolate import CubicSpline

from .utils import get_logger, require_dependency

LOG = get_logger(__name__)


@dataclass
class ContinuumSettings:
    """Settings controlling continuum fitting.

    Parameters
    ----------
    mode
        Optional continuum mode override. ``"morph_spline"`` first estimates
        a morphology baseline and then fits a spline through anchor samples on
        that baseline. ``"spline"`` skips the morphology stage and samples the
        input spectrum directly. `None` lets the mapper use the feature
        definition's ``continuum_mode``.
    morph_half_window
        Half-window, in spectral pixels, passed to
        :func:`pybaselines.morphological.mor`.
    min_finite
        Minimum number of finite samples required for a fit.
    min_anchor_points
        Minimum number of anchor wavelengths needed before using the fallback
        anchor grid.
    fallback_anchor_count
        Number of evenly spaced fallback anchors when user anchors are sparse.
    include_edge_anchors
        Add the first and last finite cube wavelengths as spline anchors.
    """

    mode: str | None = None
    morph_half_window: int = 10
    min_finite: int = 10
    min_anchor_points: int = 4
    fallback_anchor_count: int = 4
    include_edge_anchors: bool = True

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> ContinuumSettings:
        """Create continuum settings from a configuration mapping."""

        settings = cls()
        if "continuum_mode" in payload:
            settings.mode = str(payload["continuum_mode"])
        elif "mode" in payload:
            settings.mode = str(payload["mode"])
        if "morph_half_window" in payload:
            settings.morph_half_window = int(payload["morph_half_window"])
        if "min_finite" in payload:
            settings.min_finite = int(payload["min_finite"])
        if "min_anchor_points" in payload:
            settings.min_anchor_points = int(payload["min_anchor_points"])
        if "fallback_anchor_count" in payload:
            settings.fallback_anchor_count = int(payload["fallback_anchor_count"])
        if "include_edge_anchors" in payload:
            settings.include_edge_anchors = bool(payload["include_edge_anchors"])
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate continuum settings."""

        if self.mode is not None and self.mode not in {"morph_spline", "spline"}:
            raise ValueError(f"Unsupported continuum mode: {self.mode!r}.")
        if self.morph_half_window < 1:
            raise ValueError("morph_half_window must be >= 1.")
        if self.min_finite < 1:
            raise ValueError("min_finite must be >= 1.")
        if self.min_anchor_points < 1:
            raise ValueError("min_anchor_points must be >= 1.")
        if self.fallback_anchor_count < 2:
            raise ValueError("fallback_anchor_count must be >= 2.")


@dataclass
class ContinuumFit:
    """Continuum fit result for one spectrum.

    Attributes
    ----------
    continuum
        Fitted spline continuum.
    baseline
        Morphology-cleaned baseline used to sample spline anchors.
    residual
        Original input flux minus fitted continuum.
    anchor_wavelengths
        Wavelengths used in the final spline fit, in microns.
    anchor_fluxes
        Baseline flux values at the final spline anchors.
    warnings
        Non-fatal warnings produced during fitting.
    """

    continuum: np.ndarray
    baseline: np.ndarray
    residual: np.ndarray
    anchor_wavelengths: np.ndarray
    anchor_fluxes: np.ndarray
    warnings: list[str] = field(default_factory=list)


def fit_continuum(
    wavelengths_um: np.ndarray,
    flux: np.ndarray,
    anchors: tuple[tuple[float, float], ...] | list[tuple[float, float]],
    anchor_points: tuple[float, ...] | list[float] | None = None,
    settings: ContinuumSettings | None = None,
) -> ContinuumFit:
    """Fit a continuum with morphology baseline estimation and spline anchors.

    Parameters
    ----------
    wavelengths_um
        Spectral axis in microns. An Astropy `~astropy.units.Quantity` is also
        accepted and converted to microns.
    flux
        One-dimensional flux spectrum.
    anchors
        Continuum anchor windows in microns.
    anchor_points
        Discrete continuum anchor wavelengths in microns.
    settings
        Continuum settings.

    Returns
    -------
    ContinuumFit
        Continuum fit products for the input spectrum.

    The default mode first estimates a morphology-cleaned baseline-like
    spectrum with ``pybaselines.morphological.mor`` to suppress narrow emission
    lines, then fits the final spline continuum to that spectrum at the
    user-defined anchor wavelength samples. The residual is always the input
    spectrum minus the spline continuum.
    """

    settings = settings or ContinuumSettings()
    wl = _wavelength_values_um(wavelengths_um)
    y = np.asarray(flux, dtype=float)
    warnings: list[str] = []
    mode = settings.mode or "morph_spline"
    if mode not in {"morph_spline", "spline"}:
        raise ValueError(f"Unsupported continuum mode: {mode!r}.")

    finite = np.isfinite(wl) & np.isfinite(y)
    if np.count_nonzero(finite) < settings.min_finite:
        msg = "Too few finite spectral samples for continuum fitting."
        warnings.append(msg)
        median = np.nanmedian(y) if np.isfinite(y).any() else np.nan
        baseline = np.full_like(y, median, dtype=float)
        continuum = baseline.copy()
        return ContinuumFit(continuum, baseline, y - continuum, np.array([]), np.array([]), warnings)

    y_filled = _interpolate_nans(wl, y)
    if mode == "morph_spline":
        baseline, baseline_warnings = _morphology_baseline(y_filled, settings.morph_half_window)
        warnings.extend(baseline_warnings)
    else:
        baseline = y_filled.copy()

    anchor_wl, anchor_flux = _anchor_samples(wl, baseline, anchors)
    point_wl, point_flux = _point_anchor_samples(wl, baseline, anchor_points or ())
    anchor_wl = np.concatenate([anchor_wl, point_wl])
    anchor_flux = np.concatenate([anchor_flux, point_flux])
    if settings.include_edge_anchors:
        edge_wl, edge_flux = _edge_anchor_samples(wl, baseline)
        anchor_wl = np.concatenate([anchor_wl, edge_wl])
        anchor_flux = np.concatenate([anchor_flux, edge_flux])
    if anchor_wl.size < settings.min_anchor_points:
        msg = "Continuum fit poorly constrained by anchors; using fallback anchor grid."
        warnings.append(msg)
        anchor_wl = np.linspace(np.nanmin(wl[finite]), np.nanmax(wl[finite]), settings.fallback_anchor_count)
        anchor_flux = np.interp(anchor_wl, wl[finite], baseline[finite])

    anchor_wl, anchor_flux = _sort_unique_anchor_samples(anchor_wl, anchor_flux)

    try:
        if anchor_wl.size >= 4:
            spline = CubicSpline(anchor_wl, anchor_flux, bc_type="natural", extrapolate=True)
            continuum = spline(wl)
        else:
            continuum = np.interp(wl, anchor_wl, anchor_flux, left=anchor_flux[0], right=anchor_flux[-1])
    except Exception as exc:
        msg = f"Spline continuum failed; using linear interpolation: {exc}"
        warnings.append(msg)
        continuum = np.interp(wl, anchor_wl, anchor_flux, left=anchor_flux[0], right=anchor_flux[-1])

    residual = y - continuum
    return ContinuumFit(continuum=continuum, baseline=baseline, residual=residual, anchor_wavelengths=anchor_wl, anchor_fluxes=anchor_flux, warnings=warnings)


def _wavelength_values_um(wavelengths: np.ndarray | u.Quantity) -> np.ndarray:
    if isinstance(wavelengths, u.Quantity):
        return np.asarray(wavelengths.to_value(u.micron), dtype=float)
    return np.asarray(wavelengths, dtype=float)


def _morphology_baseline(flux: np.ndarray, half_window: int) -> tuple[np.ndarray, list[str]]:
    """Return a pybaselines morphology baseline with lazy dependency loading."""

    os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
    morphological = require_dependency(
        "pybaselines.morphological",
        extra="jwst",
        purpose="morphology-based continuum estimation",
    )
    try:
        baseline, _ = morphological.mor(flux, half_window=half_window)
    except Exception as exc:
        msg = f"Morphology baseline failed; using interpolated spectrum as baseline: {exc}"
        return flux.copy(), [msg]
    return np.asarray(baseline, dtype=float), []


def _interpolate_nans(wavelengths_um: np.ndarray, flux: np.ndarray) -> np.ndarray:
    finite = np.isfinite(wavelengths_um) & np.isfinite(flux)
    if np.count_nonzero(finite) == flux.size:
        return flux.astype(float, copy=True)
    filled = flux.astype(float, copy=True)
    filled[~finite] = np.interp(wavelengths_um[~finite], wavelengths_um[finite], flux[finite])
    return filled


def _anchor_samples(
    wavelengths_um: np.ndarray,
    baseline: np.ndarray,
    anchors: tuple[tuple[float, float], ...] | list[tuple[float, float]],
) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(wavelengths_um) & np.isfinite(baseline)
    wl = wavelengths_um[finite]
    y = baseline[finite]
    if wl.size == 0:
        return np.array([]), np.array([])

    anchor_wl: list[float] = []
    anchor_flux: list[float] = []
    for lo, hi in anchors:
        lo, hi = min(float(lo), float(hi)), max(float(lo), float(hi))
        mask = (wl >= lo) & (wl <= hi)
        if np.count_nonzero(mask) > 0:
            anchor_wl.extend(wl[mask].astype(float).tolist())
            anchor_flux.extend(y[mask].astype(float).tolist())
        else:
            mid = 0.5 * (lo + hi)
            if wl.min() <= mid <= wl.max():
                anchor_wl.append(mid)
                anchor_flux.append(float(np.interp(mid, wl, y)))
            else:
                LOG.warning("Continuum anchor %.4f-%.4f um lies outside spectral coverage.", lo, hi)
    anchor_wl_arr = np.asarray(anchor_wl, dtype=float)
    anchor_flux_arr = np.asarray(anchor_flux, dtype=float)
    finite_anchor = np.isfinite(anchor_wl_arr) & np.isfinite(anchor_flux_arr)
    return anchor_wl_arr[finite_anchor], anchor_flux_arr[finite_anchor]


def _point_anchor_samples(
    wavelengths_um: np.ndarray,
    baseline: np.ndarray,
    anchor_points: tuple[float, ...] | list[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Sample the morphology baseline at discrete anchor wavelengths."""

    finite = np.isfinite(wavelengths_um) & np.isfinite(baseline)
    if np.count_nonzero(finite) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    wl = wavelengths_um[finite]
    y = baseline[finite]
    sampled_wl: list[float] = []
    sampled_flux: list[float] = []
    for point in anchor_points:
        x = float(point)
        if wl.min() <= x <= wl.max():
            sampled_wl.append(x)
            sampled_flux.append(float(np.interp(x, wl, y)))
    return np.asarray(sampled_wl, dtype=float), np.asarray(sampled_flux, dtype=float)


def _edge_anchor_samples(wavelengths_um: np.ndarray, baseline: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return first and last finite cube samples as spline edge anchors."""

    finite = np.isfinite(wavelengths_um) & np.isfinite(baseline)
    if np.count_nonzero(finite) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    wl = wavelengths_um[finite]
    y = baseline[finite]
    if wl.size == 1:
        return np.asarray([wl[0]], dtype=float), np.asarray([y[0]], dtype=float)
    return np.asarray([wl[0], wl[-1]], dtype=float), np.asarray([y[0], y[-1]], dtype=float)


def _sort_unique_anchor_samples(anchor_wl: np.ndarray, anchor_flux: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort anchor samples and average repeated wavelengths for spline safety."""

    finite = np.isfinite(anchor_wl) & np.isfinite(anchor_flux)
    if np.count_nonzero(finite) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    x = np.asarray(anchor_wl[finite], dtype=float)
    y = np.asarray(anchor_flux[finite], dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    unique_x: list[float] = []
    unique_y: list[float] = []
    for value in np.unique(x):
        mask = x == value
        unique_x.append(float(value))
        unique_y.append(float(np.nanmean(y[mask])))
    return np.asarray(unique_x, dtype=float), np.asarray(unique_y, dtype=float)
