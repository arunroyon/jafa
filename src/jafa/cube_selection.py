"""Select the cube best suited for a requested spectral feature."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .feature_db import FeatureDefinition
from .io import CubeData
from .utils import get_logger

LOG = get_logger(__name__)


@dataclass(frozen=True)
class CubeSelection:
    """Result of selecting a cube for a feature."""

    cube: CubeData
    feature: FeatureDefinition
    required_window: tuple[float, float]
    spectral_margin_um: float


def select_cube_for_feature(cubes: list[CubeData], feature: FeatureDefinition) -> CubeSelection:
    """Select the narrowest cube that fully covers a feature and its anchors.

    Raises
    ------
    ValueError
        If no cube covers the feature integration window plus continuum anchor
        windows and discrete anchor points.
    """

    if not cubes:
        raise ValueError("At least one cube is required.")

    required = feature.required_window
    candidates: list[tuple[float, CubeData]] = []
    for cube in cubes:
        low, high = cube.spectral_range
        if np.isfinite(low) and np.isfinite(high) and low <= required[0] and high >= required[1]:
            coverage_width = high - low
            needed_width = required[1] - required[0]
            margin = coverage_width - needed_width
            channel_bonus = _preference_bonus(cube, feature)
            candidates.append((margin - channel_bonus, cube))
        else:
            LOG.debug(
                "Cube %s rejected for %s: coverage %.4g-%.4g um does not cover %.4g-%.4g um",
                cube.name,
                feature.feature_name,
                low,
                high,
                required[0],
                required[1],
            )

    if not candidates:
        coverage = ", ".join(f"{cube.name}: {cube.spectral_range[0]:.3f}-{cube.spectral_range[1]:.3f} um" for cube in cubes)
        raise ValueError(
            f"No cube covers {feature.feature_name} plus continuum anchors "
            f"({required[0]:.3f}-{required[1]:.3f} um). Available coverage: {coverage}"
        )

    candidates.sort(key=lambda item: item[0])
    cube = candidates[0][1]
    margin = float(cube.spectral_range[1] - cube.spectral_range[0] - (required[1] - required[0]))
    LOG.info(
        "Selected cube %s for %s (coverage %.3f-%.3f um)",
        cube.name,
        feature.feature_name,
        cube.spectral_range[0],
        cube.spectral_range[1],
    )
    return CubeSelection(cube=cube, feature=feature, required_window=required, spectral_margin_um=margin)


def choose_coarser_cube(cubes: list[CubeData], priority: list[CubeData] | None = None) -> CubeData:
    """Choose the lower-resolution/coarser cube for a common ratio grid.

    Rule order:
    1. Larger pixel scale wins.
    2. If pixel scales are similar within 5 percent, smaller spatial dimensions win.
    3. If still tied, use the optional priority order or the first cube.
    """

    if not cubes:
        raise ValueError("At least one cube is required.")
    priority = priority or cubes
    priority_index = {id(cube): idx for idx, cube in enumerate(priority)}

    def key(cube: CubeData) -> tuple[float, int, int]:
        pixel_scale = cube.pixel_scale_arcsec
        scale = pixel_scale if pixel_scale is not None and np.isfinite(pixel_scale) else -1.0
        ny, nx = cube.spatial_shape
        return (scale, -(ny * nx), -priority_index.get(id(cube), len(priority)))

    ranked = sorted(cubes, key=key, reverse=True)
    first = ranked[0]
    if len(ranked) > 1:
        second = ranked[1]
        s1, s2 = first.pixel_scale_arcsec, second.pixel_scale_arcsec
        if s1 and s2 and abs(s1 - s2) / max(s1, s2) <= 0.05:
            ranked = sorted(cubes, key=lambda c: (-(c.spatial_shape[0] * c.spatial_shape[1]), -priority_index.get(id(c), len(priority))), reverse=True)
            first = ranked[0]
    LOG.info("Selected common ratio grid from cube %s", first.name)
    return first


def _preference_bonus(cube: CubeData, feature: FeatureDefinition) -> float:
    bonus = 0.0
    if feature.preferred_instrument and cube.instrument:
        if feature.preferred_instrument.lower().split()[0] in cube.instrument.lower():
            bonus += 0.01
    if feature.preferred_channel and cube.channel:
        if str(feature.preferred_channel).lower() in str(cube.channel).lower():
            bonus += 0.01
    return bonus
