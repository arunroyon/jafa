"""Feature map creation workflows."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import yaml
from astropy.constants import c
from astropy.io import fits

from .continuum import ContinuumSettings, fit_continuum
from .cube_selection import select_cube_for_feature
from .feature_db import FeatureDefinition, load_feature_database
from .io import CubeData, binned_spatial_header, load_cube, save_fits_product
from .plotting import plot_feature_map
from .stitching import StitchSettings, stitch_adjacent_cubes_for_feature
from .utils import ensure_dir, finite_fraction, get_logger, sanitize_name
from .version import __version__

LOG = get_logger(__name__)


@dataclass
class MapSettings:
    """Settings for feature-map generation.

    Parameters
    ----------
    continuum
        Continuum-fitting settings.
    snr_threshold
        Optional signal-to-noise threshold applied when uncertainties are
        available.
    clip_negative_residuals
        Set negative integrated feature fluxes to zero.
    bin_spatial
        Integer spatial binning factor.
    progress
        Show a progress bar during map creation.
    write_continuum_map
        Also compute and write the continuum integrated over the feature
        window.
    max_nan_fraction_warn
        Warn when the output map has more than this fraction of NaNs.
    output_unit
        Output unit mode. ``"cgs"`` integrates flux density over frequency and
        writes ``erg s-1 cm-2 pixel-1`` maps. ``"native"`` preserves the older
        wavelength-integrated ``<input unit> um`` behavior.
    allow_cube_stitching
        If no single cube covers the full feature plus anchors, try stitching
        an adjacent sequence of cubes on a common spatial grid.
    stitching
        Settings for adjacent-cube stitching.
    """

    continuum: ContinuumSettings = field(default_factory=ContinuumSettings)
    snr_threshold: float | None = None
    clip_negative_residuals: bool = True
    bin_spatial: int = 1
    progress: bool = True
    write_continuum_map: bool = True
    max_nan_fraction_warn: float = 0.5
    output_unit: str = "cgs"
    allow_cube_stitching: bool = True
    stitching: StitchSettings = field(default_factory=StitchSettings)

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> MapSettings:
        """Create map settings from a configuration mapping."""

        allowed = {
            "allow_cube_stitching",
            "bin_spatial",
            "clip_negative_residuals",
            "continuum_mode",
            "fallback_anchor_count",
            "include_edge_anchors",
            "max_nan_fraction_warn",
            "min_anchor_points",
            "min_finite",
            "mode",
            "morph_half_window",
            "output_unit",
            "snr_threshold",
            "stitch_overlap_strategy",
            "stitch_reprojection_method",
            "stitch_require_common_spatial_footprint",
            "stitch_spectral_gap_tolerance_um",
            "stitch_target_grid",
            "write_continuum_map",
        }
        unknown = sorted(set(payload) - allowed)
        if unknown:
            raise ValueError(f"Unknown setting(s): {', '.join(unknown)}")

        settings = cls()
        settings.continuum = ContinuumSettings.from_mapping(payload)
        settings.stitching = StitchSettings.from_mapping(payload)
        if "snr_threshold" in payload:
            settings.snr_threshold = None if payload["snr_threshold"] is None else float(payload["snr_threshold"])
        if "clip_negative_residuals" in payload:
            settings.clip_negative_residuals = bool(payload["clip_negative_residuals"])
        if "bin_spatial" in payload:
            settings.bin_spatial = int(payload["bin_spatial"])
        if "write_continuum_map" in payload:
            settings.write_continuum_map = bool(payload["write_continuum_map"])
        if "max_nan_fraction_warn" in payload:
            settings.max_nan_fraction_warn = float(payload["max_nan_fraction_warn"])
        if "output_unit" in payload:
            settings.output_unit = str(payload["output_unit"])
        if "allow_cube_stitching" in payload:
            settings.allow_cube_stitching = bool(payload["allow_cube_stitching"])
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate map settings."""

        if self.bin_spatial < 1:
            raise ValueError("bin_spatial must be >= 1.")
        if self.output_unit not in {"cgs", "native"}:
            raise ValueError("output_unit must be 'cgs' or 'native'.")
        if self.max_nan_fraction_warn < 0 or self.max_nan_fraction_warn > 1:
            raise ValueError("max_nan_fraction_warn must be between 0 and 1.")
        self.continuum.validate()
        self.stitching.validate()


@dataclass
class FeatureMapResult:
    """Feature-map output in memory.

    Attributes
    ----------
    feature
        Feature definition used for the map.
    cube
        Selected input cube.
    feature_map
        Integrated continuum-subtracted map.
    uncertainty
        Optional propagated uncertainty map.
    continuum_map
        Optional integrated continuum map over the feature window.
    header
        Two-dimensional FITS/WCS header.
    unit
        Output unit string.
    diagnostic
        Continuum diagnostic arrays.
    metadata
        YAML-serializable metadata.
    output_paths
        Paths written by :func:`write_feature_map_outputs`.
    """

    feature: FeatureDefinition
    cube: CubeData
    feature_map: np.ndarray
    uncertainty: np.ndarray | None
    continuum_map: np.ndarray | None
    header: fits.Header
    unit: str
    diagnostic: dict[str, np.ndarray]
    metadata: dict
    output_paths: dict[str, Path] = field(default_factory=dict)


def make_feature_map(
    cubes: list[CubeData],
    feature: FeatureDefinition,
    *,
    outdir: str | Path | None = None,
    settings: MapSettings | None = None,
    write_outputs: bool = True,
) -> FeatureMapResult:
    """Create an integrated continuum-subtracted feature map."""

    settings = settings or MapSettings()
    cube = select_or_stitch_cube(cubes, feature, settings)
    result = make_feature_map_from_cube(cube, feature, settings=settings)
    if write_outputs:
        if outdir is None:
            raise ValueError("outdir is required when write_outputs=True.")
        write_feature_map_outputs(result, outdir)
    return result


def select_or_stitch_cube(cubes: list[CubeData], feature: FeatureDefinition, settings: MapSettings) -> CubeData:
    """Select one cube for a feature, stitching adjacent cubes if needed."""

    try:
        return select_cube_for_feature(cubes, feature).cube
    except ValueError as selection_error:
        if not settings.allow_cube_stitching or not settings.stitching.enabled:
            raise
        try:
            stitched = stitch_adjacent_cubes_for_feature(cubes, feature, settings=settings.stitching)
        except Exception as stitch_error:
            raise ValueError(f"{selection_error} Adjacent-cube stitching also failed: {stitch_error}") from stitch_error
        LOG.warning(
            "No single cube covers %s; using stitched cube from %s.",
            feature.feature_name,
            ", ".join(stitched.stitched_from),
        )
        return stitched


def make_feature_map_from_cube(cube: CubeData, feature: FeatureDefinition, *, settings: MapSettings | None = None) -> FeatureMapResult:
    """Create a feature map from a known cube without running cube selection."""

    settings = settings or MapSettings()
    settings.validate()
    continuum_settings = _continuum_settings_for_feature(settings, feature)

    finite_wave_mask = np.isfinite(cube.wavelength_um)
    wave_mask = _feature_mask(cube.wavelength_um, feature.integration_window)
    if np.count_nonzero(wave_mask) < 2:
        raise ValueError(f"Feature window for {feature.feature_name} has fewer than two spectral samples.")

    wavelengths = cube.wavelength_um[finite_wave_mask]
    data = cube.data[finite_wave_mask, :, :]
    uncertainty = cube.uncertainty[finite_wave_mask, :, :] if cube.uncertainty is not None else None
    local_wave_mask = _feature_mask(wavelengths, feature.integration_window)
    fit_feature, ignored_anchor_points = _feature_with_usable_anchor_points(feature, wavelengths, cube)

    ny, nx = cube.spatial_shape
    bin_size = settings.bin_spatial
    out_ny, out_nx = ny // bin_size, nx // bin_size
    fmap = np.full((out_ny, out_nx), np.nan, dtype=float)
    unc_map = np.full_like(fmap, np.nan) if uncertainty is not None else None
    cont_map = np.full_like(fmap, np.nan) if settings.write_continuum_map else None

    iterator = ((iy, ix) for iy in range(out_ny) for ix in range(out_nx))
    total = out_ny * out_nx

    for iy, ix in _progress_iter(iterator, total, settings.progress):
        ys = slice(iy * bin_size, (iy + 1) * bin_size)
        xs = slice(ix * bin_size, (ix + 1) * bin_size)
        block = data[:, ys, xs]
        finite_counts = np.count_nonzero(np.isfinite(block), axis=(1, 2)).astype(float)
        spectral_sum = np.nansum(block, axis=(1, 2))
        spectrum = np.divide(spectral_sum, finite_counts, out=np.full(data.shape[0], np.nan), where=finite_counts > 0)
        if not np.isfinite(spectrum).any():
            continue
        if finite_fraction(spectrum) < 0.3:
            continue

        fit = fit_continuum(wavelengths, spectrum, fit_feature.continuum_anchors, fit_feature.continuum_anchor_points, continuum_settings)
        residual = fit.residual
        npix = block.shape[1] * block.shape[2]
        # The continuum is fit to the mean spectrum, but the map product is an
        # integrated flux over the finite pixels in each wavelength slice. The
        # count can vary with wavelength near masks/common stitched footprints.
        local_counts = finite_counts[local_wave_mask]
        value = _integrate_flux(
            wavelengths[local_wave_mask],
            residual[local_wave_mask],
            cube,
            npix=npix,
            pixel_counts=local_counts,
            output_unit=settings.output_unit,
        )
        if settings.clip_negative_residuals:
            value = max(value, 0.0)

        fmap[iy, ix] = value
        if cont_map is not None:
            cont_map[iy, ix] = _integrate_flux(
                wavelengths[local_wave_mask],
                fit.continuum[local_wave_mask],
                cube,
                npix=npix,
                pixel_counts=local_counts,
                output_unit=settings.output_unit,
            )

        if uncertainty is not None and unc_map is not None:
            err_block = uncertainty[:, ys, xs]
            valid_err = np.isfinite(err_block) & np.isfinite(block)
            finite_err_counts = np.count_nonzero(valid_err, axis=(1, 2)).astype(float)
            err_sum = np.nansum(np.where(valid_err, err_block, np.nan) ** 2, axis=(1, 2))
            sigma_mean = np.divide(
                np.sqrt(err_sum),
                finite_err_counts,
                out=np.full(uncertainty.shape[0], np.nan),
                where=finite_err_counts > 0,
            )
            unc_map[iy, ix] = _integrated_uncertainty(
                wavelengths[local_wave_mask],
                sigma_mean[local_wave_mask],
                cube,
                npix=npix,
                pixel_counts=finite_err_counts[local_wave_mask],
                output_unit=settings.output_unit,
            )

    if settings.snr_threshold is not None and unc_map is not None:
        snr = np.divide(fmap, unc_map, out=np.full_like(fmap, np.nan), where=unc_map > 0)
        fmap[snr < settings.snr_threshold] = np.nan
        if cont_map is not None:
            cont_map[snr < settings.snr_threshold] = np.nan

    nan_fraction = 1.0 - finite_fraction(fmap)
    if nan_fraction > settings.max_nan_fraction_warn:
        LOG.warning("%.1f percent of %s map pixels are NaN.", 100 * nan_fraction, feature.feature_name)

    header = binned_spatial_header(cube.spatial_header, bin_size)
    unit = _integrated_unit(cube.flux_unit, cube.pixel_area_sr is not None, settings.output_unit)
    diagnostic = _diagnostic_fit(cube, feature, fit_feature, finite_wave_mask, wavelengths, settings, continuum_settings)
    metadata = _metadata(cube, feature, settings, unit, ignored_anchor_points, continuum_mode=continuum_settings.mode)

    return FeatureMapResult(
        feature=feature,
        cube=cube,
        feature_map=fmap,
        uncertainty=unc_map,
        continuum_map=cont_map,
        header=header,
        unit=unit,
        diagnostic=diagnostic,
        metadata=metadata,
    )


def write_feature_map_outputs(result: FeatureMapResult, outdir: str | Path) -> dict[str, Path]:
    """Write FITS, plot, and metadata products for one feature map."""

    outdir = ensure_dir(outdir)
    stem = sanitize_name(result.feature.feature_name)
    paths: dict[str, Path] = {}
    common_meta = {
        "FEATURE": result.feature.feature_name,
        "LAM0": float(result.feature.central_wavelength),
        "LWMIN": float(result.feature.integration_window[0]),
        "LWMAX": float(result.feature.integration_window[1]),
    }
    paths["map"] = save_fits_product(result.feature_map, result.header, outdir / f"{stem}_map.fits", bunit=result.unit, metadata=common_meta)
    if result.uncertainty is not None:
        paths["uncertainty"] = save_fits_product(result.uncertainty, result.header, outdir / f"{stem}_unc.fits", bunit=result.unit, metadata=common_meta)
    if result.continuum_map is not None:
        paths["continuum"] = save_fits_product(result.continuum_map, result.header, outdir / f"{stem}_continuum.fits", bunit=result.unit, metadata=common_meta)
    paths["diagnostic"] = plot_feature_map(
        result.feature_map,
        result.uncertainty,
        outdir / f"{stem}_diagnostic.png",
        title=f"{result.feature.feature_name} integrated map",
        unit_label=result.unit,
        diagnostic=result.diagnostic,
        header=result.header,
    )
    metadata_path = outdir / f"{stem}_metadata.yaml"
    with metadata_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(result.metadata, handle, sort_keys=False)
    paths["metadata"] = metadata_path
    result.output_paths.update(paths)
    return paths


class FeatureMapper:
    """High-level interface for creating feature and ratio maps."""

    def __init__(self, cubes: Iterable[str | Path | CubeData], *, feature_db_path: str | Path | None = None, settings: MapSettings | None = None):
        self.cubes = [cube if isinstance(cube, CubeData) else load_cube(cube) for cube in cubes]
        self.features = load_feature_database(feature_db_path)
        self.settings = settings or MapSettings()

    def make_feature_map(self, feature_name: str, *, outdir: str | Path | None = None, write_outputs: bool = True) -> FeatureMapResult:
        """Create a map for a named feature from the feature database."""

        feature = self.features[feature_name]
        return make_feature_map(self.cubes, feature, outdir=outdir, settings=self.settings, write_outputs=write_outputs)

    def make_ratio_map(self, feature1: str, feature2: str, *, outdir: str | Path, write_outputs: bool = True):
        """Create a ratio map for two named features."""

        from .ratios import make_ratio_map

        return make_ratio_map(self.cubes, self.features[feature1], self.features[feature2], outdir=outdir, settings=self.settings, write_outputs=write_outputs)


def _feature_mask(wavelengths: np.ndarray, window: tuple[float, float]) -> np.ndarray:
    wavelengths = _wavelength_values_um(wavelengths)
    lo, hi = min(window), max(window)
    return (wavelengths >= lo) & (wavelengths <= hi)


def _integrate(wavelengths: np.ndarray, values: np.ndarray) -> float:
    wavelengths = _wavelength_values_um(wavelengths)
    finite = np.isfinite(wavelengths) & np.isfinite(values)
    if np.count_nonzero(finite) < 2:
        return np.nan
    return float(np.trapz(values[finite], wavelengths[finite]))


def _integrate_scaled(wavelengths: np.ndarray, values: np.ndarray, scale: np.ndarray) -> float:
    wavelengths = _wavelength_values_um(wavelengths)
    y = np.asarray(values, dtype=float)
    s = np.asarray(scale, dtype=float)
    finite = np.isfinite(wavelengths) & np.isfinite(y) & np.isfinite(s)
    if np.count_nonzero(finite) < 2:
        return np.nan
    x = wavelengths[finite]
    integrand = y[finite] * s[finite]
    order = np.argsort(x)
    return float(np.trapz(integrand[order], x[order]))


def _integrate_flux(
    wavelengths: np.ndarray,
    values: np.ndarray,
    cube: CubeData,
    *,
    npix: int,
    pixel_counts: np.ndarray | None = None,
    output_unit: str,
) -> float:
    if output_unit == "native":
        scale = _area_factors(cube, npix, pixel_counts, len(_wavelength_values_um(wavelengths)))
        return _integrate_scaled(wavelengths, values, scale)
    return _integrate_cgs(wavelengths, values, cube, npix=npix, pixel_counts=pixel_counts)


def _integrate_cgs(wavelengths: np.ndarray, values: np.ndarray, cube: CubeData, *, npix: int, pixel_counts: np.ndarray | None = None) -> float:
    """Integrate Fnu over frequency and return cgs flux per output pixel."""

    wl = _wavelength_values_um(wavelengths)
    y = np.asarray(values, dtype=float)
    area_factors = _area_factors(cube, npix, pixel_counts, wl.size)
    finite = np.isfinite(wl) & np.isfinite(y) & np.isfinite(area_factors)
    if np.count_nonzero(finite) < 2:
        return np.nan
    if cube.pixel_area_sr is None:
        raise ValueError("Cgs output requires cube.pixel_area_sr so surface brightness can become flux per pixel.")

    wl = wl[finite]
    y = y[finite]
    area_factors = area_factors[finite]
    order = np.argsort(wl)
    wl = wl[order]
    y = y[order]
    area_factors = area_factors[order]

    unit = _surface_brightness_unit(cube.flux_unit)
    y_cgs = (y * unit).to_value(u.erg / u.s / u.cm**2 / u.Hz / u.sr)
    frequency_hz = (c / (wl * u.micron)).to_value(u.Hz)
    return float(np.trapz((y_cgs * area_factors)[::-1], frequency_hz[::-1]))


def _integrated_uncertainty(
    wavelengths: np.ndarray,
    sigma: np.ndarray,
    cube: CubeData,
    *,
    npix: int,
    pixel_counts: np.ndarray | None = None,
    output_unit: str,
) -> float:
    wavelengths = _wavelength_values_um(wavelengths)
    area_factors = _area_factors(cube, npix, pixel_counts, wavelengths.size)
    finite = np.isfinite(wavelengths) & np.isfinite(sigma) & np.isfinite(area_factors)
    x = wavelengths[finite]
    s = sigma[finite]
    area_factors = area_factors[finite]
    if x.size < 2:
        return np.nan

    if output_unit == "cgs":
        if cube.pixel_area_sr is None:
            raise ValueError("Cgs uncertainty output requires cube.pixel_area_sr.")
        unit = _surface_brightness_unit(cube.flux_unit)
        s = (s * unit).to_value(u.erg / u.s / u.cm**2 / u.Hz / u.sr)
        freq = (c / (x * u.micron)).to_value(u.Hz)
        order = np.argsort(freq)
        x = freq[order]
        s = s[order]
        area_factors = area_factors[order]
    else:
        order = np.argsort(x)
        x = x[order]
        s = s[order]
        area_factors = area_factors[order]

    weights = np.empty_like(x)
    weights[0] = 0.5 * (x[1] - x[0])
    weights[-1] = 0.5 * (x[-1] - x[-2])
    if x.size > 2:
        weights[1:-1] = 0.5 * (x[2:] - x[:-2])
    return float(np.sqrt(np.nansum((s * weights * area_factors) ** 2)))


def _area_factor(cube: CubeData, npix: int) -> float:
    if cube.pixel_area_sr is None:
        return 1.0
    return float(cube.pixel_area_sr * npix)


def _area_factors(cube: CubeData, npix: int, pixel_counts: np.ndarray | None, size: int) -> np.ndarray:
    if cube.pixel_area_sr is None:
        return np.ones(size, dtype=float)
    if pixel_counts is None:
        return np.full(size, _area_factor(cube, npix), dtype=float)
    counts = np.asarray(pixel_counts, dtype=float)
    if counts.size != size:
        raise ValueError("pixel_counts must have the same length as wavelengths.")
    return counts * float(cube.pixel_area_sr)


def _integrated_unit(flux_unit: str, area_scaled: bool, output_unit: str) -> str:
    if output_unit == "cgs":
        return "erg s-1 cm-2 pixel-1"
    unit = flux_unit.strip() or "flux"
    if area_scaled and "/sr" in unit:
        unit = unit.replace("/sr", "").replace("sr-1", "").strip()
    return f"{unit} um"


def _surface_brightness_unit(flux_unit: str) -> u.UnitBase:
    try:
        unit = u.Unit(flux_unit)
    except Exception as exc:
        raise ValueError(f"Cannot parse input flux unit {flux_unit!r} for cgs conversion.") from exc
    target = u.erg / u.s / u.cm**2 / u.Hz / u.sr
    if not unit.is_equivalent(target):
        raise ValueError(
            f"Cgs output requires a flux-density surface-brightness unit equivalent to {target}; "
            f"got {flux_unit!r}."
        )
    return unit


def _wavelength_values_um(wavelengths: np.ndarray | u.Quantity) -> np.ndarray:
    if isinstance(wavelengths, u.Quantity):
        return np.asarray(wavelengths.to_value(u.micron), dtype=float)
    return np.asarray(wavelengths, dtype=float)


def _progress_iter(iterator, total: int, enabled: bool):
    if not enabled:
        yield from iterator
        return
    from astropy.utils.console import ProgressBar

    with ProgressBar(total) as bar:
        for item in iterator:
            yield item
            bar.update()


def _continuum_settings_for_feature(settings: MapSettings, feature: FeatureDefinition) -> ContinuumSettings:
    mode = settings.continuum.mode or feature.continuum_mode
    return replace(settings.continuum, mode=mode)


def _feature_with_usable_anchor_points(
    feature: FeatureDefinition,
    wavelengths: np.ndarray,
    cube: CubeData,
) -> tuple[FeatureDefinition, tuple[float, ...]]:
    """Return a feature definition whose point anchors lie inside this cube."""

    finite = wavelengths[np.isfinite(wavelengths)]
    if finite.size == 0 or not feature.continuum_anchor_points:
        return feature, ()

    lo = float(np.nanmin(finite))
    hi = float(np.nanmax(finite))
    usable: list[float] = []
    ignored: list[float] = []
    for point in feature.continuum_anchor_points:
        x = float(point)
        if lo <= x <= hi:
            usable.append(x)
        else:
            ignored.append(x)

    if ignored:
        LOG.warning(
            "Ignoring %s anchor point(s) outside %s coverage %.4f-%.4f um for %s: %s. "
            "Cube edge anchors are still included when enabled.",
            len(ignored),
            cube.name,
            lo,
            hi,
            feature.feature_name,
            ", ".join(f"{x:.4f}" for x in ignored),
        )
    if len(usable) == len(feature.continuum_anchor_points):
        return feature, ()
    return replace(feature, continuum_anchor_points=tuple(usable)), tuple(ignored)


def _diagnostic_fit(
    cube: CubeData,
    feature: FeatureDefinition,
    fit_feature: FeatureDefinition,
    work_mask: np.ndarray,
    wavelengths: np.ndarray,
    settings: MapSettings,
    continuum_settings: ContinuumSettings,
) -> dict[str, np.ndarray]:
    data = cube.data[work_mask, :, :]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spectrum = np.nanmedian(data, axis=(1, 2))
    fit = fit_continuum(wavelengths, spectrum, fit_feature.continuum_anchors, fit_feature.continuum_anchor_points, continuum_settings)
    return {
        "wavelength_um": wavelengths,
        "spectrum": spectrum,
        "baseline": fit.baseline,
        "continuum": fit.continuum,
        "residual": fit.baseline - fit.continuum,
        "anchor_wavelengths": fit.anchor_wavelengths,
        "anchor_fluxes": fit.anchor_fluxes,
        "continuum_anchors": np.asarray(feature.continuum_anchors, dtype=float),
        "continuum_anchor_points": np.asarray(fit_feature.continuum_anchor_points, dtype=float),
        "requested_continuum_anchor_points": np.asarray(feature.continuum_anchor_points, dtype=float),
        "feature_window": np.asarray(feature.integration_window, dtype=float),
        "morph_half_window": np.asarray([settings.continuum.morph_half_window], dtype=float),
    }


def _metadata(
    cube: CubeData,
    feature: FeatureDefinition,
    settings: MapSettings,
    unit: str,
    ignored_anchor_points: tuple[float, ...] = (),
    *,
    continuum_mode: str | None = None,
) -> dict:
    return {
        "software": {"name": "JAFA", "full_name": "JWST Aromatic Feature Analyzer", "package": "jafa", "version": __version__},
        "cube": {
            "path": str(cube.path) if cube.path else None,
            "name": cube.name,
            "spectral_range_um": list(cube.spectral_range),
            "flux_unit": cube.flux_unit,
            "pixel_area_sr": cube.pixel_area_sr,
            "instrument": cube.instrument,
            "channel": cube.channel,
            "stitched_from": list(cube.stitched_from),
            "stitch_target": cube.stitch_target,
            "stitch_overlap_strategy": cube.stitch_overlap_strategy,
            "stitch_overlap_ranges_um": [list(item) for item in cube.stitch_overlap_ranges_um],
            "stitch_common_footprint_fraction": cube.stitch_common_footprint_fraction,
        },
        "feature": feature.to_metadata(),
        "ignored_anchor_points_um": list(ignored_anchor_points),
        "settings": {
            "morph_half_window": settings.continuum.morph_half_window,
            "continuum_mode": continuum_mode or settings.continuum.mode or feature.continuum_mode,
            "min_finite": settings.continuum.min_finite,
            "min_anchor_points": settings.continuum.min_anchor_points,
            "fallback_anchor_count": settings.continuum.fallback_anchor_count,
            "include_edge_anchors": settings.continuum.include_edge_anchors,
            "snr_threshold": settings.snr_threshold,
            "clip_negative_residuals": settings.clip_negative_residuals,
            "bin_spatial": settings.bin_spatial,
            "continuum_fit_span": "full_cube",
            "output_unit_mode": settings.output_unit,
            "spectral_integration": "frequency_integral_exact" if settings.output_unit == "cgs" else "wavelength_integral_native",
            "allow_cube_stitching": settings.allow_cube_stitching,
            "stitch_target_grid": settings.stitching.target_grid,
            "stitch_spectral_gap_tolerance_um": settings.stitching.spectral_gap_tolerance_um,
            "stitch_overlap_strategy": settings.stitching.overlap_strategy,
            "stitch_require_common_spatial_footprint": settings.stitching.require_common_spatial_footprint,
        },
        "output_unit": unit,
    }
