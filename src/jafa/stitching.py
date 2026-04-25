"""Spectral stitching utilities for adjacent IFU cubes."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits

from .cube_selection import choose_coarser_cube
from .feature_db import FeatureDefinition
from .io import CubeData
from .reprojection import reproject_map, reproject_uncertainty
from .utils import get_logger, same_celestial_grid

LOG = get_logger(__name__)


@dataclass
class StitchSettings:
    """Settings controlling adjacent-cube stitching.

    Parameters
    ----------
    enabled
        Allow automatic stitching when no single cube covers the requested
        feature and its continuum anchors.
    spectral_gap_tolerance_um
        Maximum allowed wavelength gap between adjacent cube coverages. The
        default requires overlapping or touching coverage.
    target_grid
        Spatial grid for the stitched cube. ``"coarser"`` uses the
        lower-resolution/lower-pixel-count cube, matching ratio-map behavior.
        ``"first"`` uses the first cube in the selected adjacent pair.
    overlap_strategy
        Spectral overlap handling. ``"split"`` uses one cube on each side of
        the midpoint of every overlap, preventing duplicated overlap samples.
        ``"keep"`` keeps every sample from every cube.
    require_common_spatial_footprint
        Set pixels outside the common spatial footprint of all stitched cubes
        to NaN across every wavelength plane.
    reprojection_method
        Reprojection method passed to :func:`jafa.reprojection.reproject_map`.
    """

    enabled: bool = True
    spectral_gap_tolerance_um: float = 0.0
    target_grid: str = "coarser"
    overlap_strategy: str = "split"
    require_common_spatial_footprint: bool = True
    reprojection_method: str = "interp"

    @classmethod
    def from_mapping(cls, payload: Mapping[str, Any]) -> StitchSettings:
        """Create stitching settings from a configuration mapping."""

        settings = cls()
        if "stitch_spectral_gap_tolerance_um" in payload:
            settings.spectral_gap_tolerance_um = float(payload["stitch_spectral_gap_tolerance_um"])
        if "stitch_target_grid" in payload:
            settings.target_grid = str(payload["stitch_target_grid"])
        if "stitch_overlap_strategy" in payload:
            settings.overlap_strategy = str(payload["stitch_overlap_strategy"])
        if "stitch_require_common_spatial_footprint" in payload:
            settings.require_common_spatial_footprint = bool(payload["stitch_require_common_spatial_footprint"])
        if "stitch_reprojection_method" in payload:
            settings.reprojection_method = str(payload["stitch_reprojection_method"])
        settings.validate()
        return settings

    def validate(self) -> None:
        """Validate stitching settings."""

        if self.spectral_gap_tolerance_um < 0:
            raise ValueError("stitch_spectral_gap_tolerance_um must be >= 0.")
        if self.target_grid not in {"coarser", "first"}:
            raise ValueError("stitch_target_grid must be 'coarser' or 'first'.")
        if self.overlap_strategy not in {"split", "keep"}:
            raise ValueError("stitch_overlap_strategy must be 'split' or 'keep'.")
        if self.reprojection_method != "interp":
            raise ValueError("Only interpolation reprojection is currently implemented.")


def stitch_adjacent_cubes_for_feature(
    cubes: Sequence[CubeData],
    feature: FeatureDefinition,
    *,
    settings: StitchSettings | None = None,
) -> CubeData:
    """Create an in-memory stitched cube covering a requested feature.

    The selected adjacent cube sequence must collectively cover the feature integration
    window, continuum anchor windows, and discrete anchor points. Each spectral
    plane is reprojected to a common celestial grid before the wavelength axis
    is concatenated, sorted, and duplicate wavelength samples are averaged.

    Parameters
    ----------
    cubes
        Candidate input cubes.
    feature
        Feature definition whose required wavelength span should be covered.
    settings
        Stitching settings.

    Returns
    -------
    CubeData
        Stitched cube ordered as ``(wavelength, y, x)``.

    Raises
    ------
    ValueError
        If no adjacent cube sequence can cover the required wavelength span.
    """

    settings = settings or StitchSettings()
    sequence = select_adjacent_cube_sequence(cubes, feature.required_window, settings=settings)
    return stitch_cubes(sequence, settings=settings, name=f"stitched_{feature.feature_name}")


def select_adjacent_cube_sequence(
    cubes: Sequence[CubeData],
    required_window: tuple[float, float],
    *,
    settings: StitchSettings | None = None,
) -> tuple[CubeData, ...]:
    """Select the narrowest adjacent cube sequence covering ``required_window``."""

    settings = settings or StitchSettings()
    if len(cubes) < 2:
        raise ValueError("At least two cubes are required for adjacent-cube stitching.")

    required_lo, required_hi = min(required_window), max(required_window)
    ordered_cubes = tuple(sorted((cube for cube in cubes if _finite_range(cube)), key=lambda cube: cube.spectral_range[0]))
    candidates: list[tuple[float, int, tuple[CubeData, ...]]] = []
    for start in range(len(ordered_cubes)):
        current: list[CubeData] = []
        max_gap = 0.0
        for stop in range(start, len(ordered_cubes)):
            cube = ordered_cubes[stop]
            if current:
                previous_hi = current[-1].spectral_range[1]
                gap = cube.spectral_range[0] - previous_hi
                if gap > settings.spectral_gap_tolerance_um:
                    break
                max_gap = max(max_gap, gap)
            current.append(cube)
            if len(current) < 2:
                continue
            combined_lo = min(item.spectral_range[0] for item in current)
            combined_hi = max(item.spectral_range[1] for item in current)
            if combined_lo <= required_lo and combined_hi >= required_hi:
                excess_width = (combined_hi - combined_lo) - (required_hi - required_lo)
                candidates.append((excess_width + max(max_gap, 0.0), len(current), tuple(current)))
                break

    if not candidates:
        coverage = ", ".join(f"{cube.name}: {cube.spectral_range[0]:.3f}-{cube.spectral_range[1]:.3f} um" for cube in cubes)
        raise ValueError(
            "No adjacent cube sequence covers "
            f"{required_lo:.3f}-{required_hi:.3f} um with gap tolerance "
            f"{settings.spectral_gap_tolerance_um:.4g} um. Available coverage: {coverage}"
        )

    candidates.sort(key=lambda item: (item[0], item[1]))
    selected = candidates[0][2]
    LOG.info(
        "Selected adjacent cubes %s for stitched coverage %.3f-%.3f um.",
        " + ".join(cube.name for cube in selected),
        required_lo,
        required_hi,
    )
    return selected


def select_adjacent_cube_pair(
    cubes: Sequence[CubeData],
    required_window: tuple[float, float],
    *,
    settings: StitchSettings | None = None,
) -> tuple[CubeData, CubeData]:
    """Select an adjacent two-cube sequence covering ``required_window``."""

    selected = select_adjacent_cube_sequence(cubes, required_window, settings=settings)
    if len(selected) != 2:
        raise ValueError(f"Selected adjacent sequence contains {len(selected)} cubes, not two.")
    return selected[0], selected[1]


def _finite_range(cube: CubeData) -> bool:
    low, high = cube.spectral_range
    return bool(np.isfinite(low) and np.isfinite(high))


def stitch_cubes(
    cubes: Sequence[CubeData],
    *,
    settings: StitchSettings | None = None,
    name: str | None = None,
) -> CubeData:
    """Stitch adjacent cubes on a common spatial grid.

    Parameters
    ----------
    cubes
        Adjacent cubes to stitch. The current implementation accepts a sequence
        of neighboring MIRI MRS sub-band cubes with compatible flux units.
    settings
        Stitching settings.
    name
        Optional display name for the output cube.

    Returns
    -------
    CubeData
        Stitched in-memory cube.
    """

    settings = settings or StitchSettings()
    if len(cubes) < 2:
        raise ValueError("At least two cubes are required for stitching.")
    ordered = tuple(sorted(cubes, key=lambda cube: cube.spectral_range[0]))
    target = _target_cube(ordered, settings)
    target_header = target.spatial_header
    target_shape = target.spatial_shape
    target_flux_unit = _reference_flux_unit(ordered)
    keep_masks, overlap_ranges = _spectral_keep_masks(ordered, settings.overlap_strategy)

    wavelengths: list[np.ndarray] = []
    data_blocks: list[np.ndarray] = []
    uncertainty_blocks: list[np.ndarray] = []
    spatial_footprints: list[np.ndarray] = []
    have_all_uncertainties = all(cube.uncertainty is not None for cube in ordered)

    for cube, keep_mask in zip(ordered, keep_masks, strict=True):
        scale = _unit_scale(cube.flux_unit, target_flux_unit)
        data, uncertainty, spatial_footprint = _cube_on_target_grid(cube, target, target_header, target_shape, settings)
        wavelengths.append(np.asarray(cube.wavelength_um, dtype=float)[keep_mask])
        data_blocks.append(data[keep_mask] * scale)
        spatial_footprints.append(spatial_footprint)
        if have_all_uncertainties and uncertainty is not None:
            uncertainty_blocks.append(uncertainty[keep_mask] * scale)
        elif cube.uncertainty is not None:
            LOG.warning("Dropping stitched uncertainty because not all source cubes provide uncertainty.")

    common_footprint_fraction = None
    if settings.require_common_spatial_footprint:
        common_spatial_footprint = np.logical_and.reduce(spatial_footprints)
        if not np.any(common_spatial_footprint):
            raise ValueError("Stitched cubes have no common spatial footprint.")
        common_footprint_fraction = float(np.count_nonzero(common_spatial_footprint) / common_spatial_footprint.size)
        for block in data_blocks:
            block[:, ~common_spatial_footprint] = np.nan
        for block in uncertainty_blocks:
            block[:, ~common_spatial_footprint] = np.nan
        LOG.info("Retained %.1f percent of target pixels in stitched common footprint.", 100.0 * common_footprint_fraction)

    wavelength_um = np.concatenate(wavelengths)
    data = np.concatenate(data_blocks, axis=0)
    uncertainty = np.concatenate(uncertainty_blocks, axis=0) if have_all_uncertainties else None
    wavelength_um, data, uncertainty = _sort_and_merge_duplicate_wavelengths(wavelength_um, data, uncertainty)

    source_names = tuple(cube.name for cube in ordered)
    LOG.info(
        "Created stitched cube %s from %s with coverage %.3f-%.3f um.",
        name or "+".join(source_names),
        ", ".join(source_names),
        float(np.nanmin(wavelength_um)),
        float(np.nanmax(wavelength_um)),
    )

    return CubeData(
        path=None,
        data=data,
        uncertainty=uncertainty,
        wavelength_um=wavelength_um,
        wcs_2d=target.wcs_2d,
        header=target.header,
        flux_unit=target_flux_unit,
        uncertainty_unit=target.uncertainty_unit,
        pixel_area_sr=target.pixel_area_sr,
        instrument=target.instrument,
        channel="+".join(str(cube.channel or cube.name) for cube in ordered),
        name_hint=name or f"stitched_{'_'.join(source_names)}",
        stitched_from=source_names,
        stitch_target=target.name,
        stitch_overlap_strategy=settings.overlap_strategy,
        stitch_overlap_ranges_um=overlap_ranges,
        stitch_common_footprint_fraction=common_footprint_fraction,
    )


def _target_cube(cubes: Sequence[CubeData], settings: StitchSettings) -> CubeData:
    if settings.target_grid == "coarser":
        return choose_coarser_cube(list(cubes), priority=list(cubes))
    if settings.target_grid == "first":
        return cubes[0]
    raise ValueError("target_grid must be 'coarser' or 'first'.")


def _reference_flux_unit(cubes: Sequence[CubeData]) -> str:
    for cube in cubes:
        if cube.flux_unit:
            return cube.flux_unit
    return ""


def _spectral_keep_masks(
    cubes: Sequence[CubeData],
    overlap_strategy: str,
) -> tuple[tuple[np.ndarray, ...], tuple[tuple[float, float], ...]]:
    if overlap_strategy not in {"split", "keep"}:
        raise ValueError("overlap_strategy must be 'split' or 'keep'.")

    keep_masks = [np.ones(cube.wavelength_um.size, dtype=bool) for cube in cubes]
    overlap_ranges: list[tuple[float, float]] = []
    if overlap_strategy == "keep":
        return tuple(keep_masks), tuple(overlap_ranges)

    for idx in range(len(cubes) - 1):
        blue = cubes[idx]
        red = cubes[idx + 1]
        blue_lo, blue_hi = blue.spectral_range
        red_lo, red_hi = red.spectral_range
        overlap_lo = max(blue_lo, red_lo)
        overlap_hi = min(blue_hi, red_hi)
        if overlap_hi <= overlap_lo:
            continue
        split = 0.5 * (overlap_lo + overlap_hi)
        keep_masks[idx] &= blue.wavelength_um <= split
        keep_masks[idx + 1] &= red.wavelength_um > split
        overlap_ranges.append((float(overlap_lo), float(overlap_hi)))

    for cube, keep_mask in zip(cubes, keep_masks, strict=True):
        if np.count_nonzero(keep_mask) < 2:
            raise ValueError(f"Overlap strategy {overlap_strategy!r} removed too many spectral samples from {cube.name}.")
    return tuple(keep_masks), tuple(overlap_ranges)


def _unit_scale(source_unit_text: str, target_unit_text: str) -> float:
    if not source_unit_text or not target_unit_text or source_unit_text == target_unit_text:
        return 1.0
    try:
        source_unit = u.Unit(source_unit_text)
        target_unit = u.Unit(target_unit_text)
    except Exception as exc:
        raise ValueError(f"Cannot stitch cubes with incompatible flux units {source_unit_text!r} and {target_unit_text!r}.") from exc
    if not source_unit.is_equivalent(target_unit):
        raise ValueError(f"Cannot stitch cubes with incompatible flux units {source_unit_text!r} and {target_unit_text!r}.")
    return float((1.0 * source_unit).to_value(target_unit))


def _cube_on_target_grid(
    cube: CubeData,
    target: CubeData,
    target_header: fits.Header,
    target_shape: tuple[int, int],
    settings: StitchSettings,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
    if cube is target:
        return cube.data, cube.uncertainty, _native_spatial_footprint(cube.data)
    if _same_grid(cube.spatial_header, target_header, cube.spatial_shape, target_shape):
        return cube.data, cube.uncertainty, _native_spatial_footprint(cube.data)

    source_header = cube.spatial_header
    data = np.full((cube.data.shape[0], *target_shape), np.nan, dtype=float)
    footprints = np.full_like(data, np.nan)
    for idx, plane in enumerate(cube.data):
        reprojected = reproject_map(plane, source_header, target_header, target_shape, method=settings.reprojection_method)
        data[idx] = reprojected.data
        footprints[idx] = reprojected.footprint

    uncertainty = None
    if cube.uncertainty is not None:
        uncertainty = np.full_like(data, np.nan)
        for idx, plane in enumerate(cube.uncertainty):
            uncertainty[idx] = reproject_uncertainty(plane, source_header, target_header, target_shape, method=settings.reprojection_method).data
    return data, uncertainty, np.nanmedian(footprints, axis=0) > 0


def _native_spatial_footprint(data: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return np.any(np.isfinite(data), axis=0)


def _same_grid(
    first_header: fits.Header,
    second_header: fits.Header,
    first_shape: tuple[int, int],
    second_shape: tuple[int, int],
) -> bool:
    return same_celestial_grid(first_header, first_shape, second_header, second_shape)


def _sort_and_merge_duplicate_wavelengths(
    wavelength_um: np.ndarray,
    data: np.ndarray,
    uncertainty: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    order = np.argsort(wavelength_um)
    wavelength_um = np.asarray(wavelength_um[order], dtype=float)
    data = data[order]
    if uncertainty is not None:
        uncertainty = uncertainty[order]

    unique, inverse, counts = np.unique(wavelength_um, return_inverse=True, return_counts=True)
    if np.all(counts == 1):
        return wavelength_um, data, uncertainty

    merged_data = np.full((unique.size, *data.shape[1:]), np.nan, dtype=float)
    merged_uncertainty = np.full_like(merged_data, np.nan) if uncertainty is not None else None
    for idx in range(unique.size):
        mask = inverse == idx
        with np.errstate(invalid="ignore"):
            merged_data[idx] = np.nanmean(data[mask], axis=0)
        if uncertainty is not None and merged_uncertainty is not None:
            merged_uncertainty[idx] = np.sqrt(np.nansum(uncertainty[mask] ** 2, axis=0)) / np.count_nonzero(mask)
    return unique, merged_data, merged_uncertainty
