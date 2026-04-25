"""Ratio-map workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import yaml
from astropy.io import fits

from .cube_selection import choose_coarser_cube
from .feature_db import FeatureDefinition
from .io import CubeData, save_fits_product
from .mapping import (
    FeatureMapResult,
    MapSettings,
    make_feature_map_from_cube,
    select_or_stitch_cube,
    write_feature_map_outputs,
)
from .plotting import plot_ratio_map
from .reprojection import reproject_map, reproject_uncertainty
from .utils import ensure_dir, get_logger, same_celestial_grid, sanitize_name
from .version import __version__

LOG = get_logger(__name__)


@dataclass
class RatioMapResult:
    """Ratio-map output in memory."""

    feature1: FeatureDefinition
    feature2: FeatureDefinition
    numerator: FeatureMapResult
    denominator: FeatureMapResult
    ratio: np.ndarray
    ratio_uncertainty: np.ndarray | None
    header: fits.Header
    metadata: dict
    output_paths: dict[str, Path] = field(default_factory=dict)


def make_ratio_map(
    cubes: list[CubeData],
    feature1: FeatureDefinition,
    feature2: FeatureDefinition,
    *,
    outdir: str | Path | None = None,
    settings: MapSettings | None = None,
    write_outputs: bool = True,
) -> RatioMapResult:
    """Create a ratio map for two requested spectral features."""

    settings = settings or MapSettings()
    selected1 = select_or_stitch_cube(cubes, feature1, settings)
    selected2 = select_or_stitch_cube(cubes, feature2, settings)

    num = make_feature_map_from_cube(selected1, feature1, settings=settings)
    den = make_feature_map_from_cube(selected2, feature2, settings=settings)

    same_grid = selected1 is selected2 or same_celestial_grid(num.header, num.feature_map.shape, den.header, den.feature_map.shape)

    reprojection_note = "same_cube_or_grid"
    if same_grid:
        num_map, den_map = num.feature_map, den.feature_map
        num_unc, den_unc = num.uncertainty, den.uncertainty
        header = num.header
    else:
        target_cube = choose_coarser_cube([selected1, selected2], priority=[selected1, selected2])
        if target_cube is selected1:
            header = num.header
            target_shape = num.feature_map.shape
            num_map, num_unc = num.feature_map, num.uncertainty
            den_map = reproject_map(den.feature_map, den.header, header, target_shape).data
            den_unc = reproject_uncertainty(den.uncertainty, den.header, header, target_shape).data if den.uncertainty is not None else None
            reprojection_note = f"{feature2.feature_name}_reprojected_to_{feature1.feature_name}_grid"
        else:
            header = den.header
            target_shape = den.feature_map.shape
            den_map, den_unc = den.feature_map, den.uncertainty
            num_map = reproject_map(num.feature_map, num.header, header, target_shape).data
            num_unc = reproject_uncertainty(num.uncertainty, num.header, header, target_shape).data if num.uncertainty is not None else None
            reprojection_note = f"{feature1.feature_name}_reprojected_to_{feature2.feature_name}_grid"

    ratio, ratio_unc = compute_ratio(
        num_map,
        den_map,
        numerator_uncertainty=num_unc,
        denominator_uncertainty=den_unc,
        snr_threshold=settings.snr_threshold,
    )
    metadata = _ratio_metadata(num, den, reprojection_note)
    result = RatioMapResult(feature1, feature2, num, den, ratio, ratio_unc, header, metadata)
    if write_outputs:
        if outdir is None:
            raise ValueError("outdir is required when write_outputs=True.")
        write_ratio_outputs(result, outdir, write_intermediate=True)
    return result


def compute_ratio(
    numerator: np.ndarray,
    denominator: np.ndarray,
    *,
    numerator_uncertainty: np.ndarray | None = None,
    denominator_uncertainty: np.ndarray | None = None,
    snr_threshold: float | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Compute a masked ratio and propagate uncertainty when available."""

    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    valid = np.isfinite(num) & np.isfinite(den) & (den > 0)

    if snr_threshold is not None:
        if numerator_uncertainty is not None:
            num_snr = np.divide(num, numerator_uncertainty, out=np.full_like(num, np.nan), where=numerator_uncertainty > 0)
            valid &= num_snr >= snr_threshold
        if denominator_uncertainty is not None:
            den_snr = np.divide(den, denominator_uncertainty, out=np.full_like(den, np.nan), where=denominator_uncertainty > 0)
            valid &= den_snr >= snr_threshold

    ratio = np.full_like(num, np.nan, dtype=float)
    ratio[valid] = num[valid] / den[valid]

    ratio_unc = None
    if numerator_uncertainty is not None and denominator_uncertainty is not None:
        nerr = np.asarray(numerator_uncertainty, dtype=float)
        derr = np.asarray(denominator_uncertainty, dtype=float)
        valid_unc = valid & (num > 0) & (nerr >= 0) & (derr >= 0)
        ratio_unc = np.full_like(num, np.nan, dtype=float)
        rel = np.sqrt((nerr[valid_unc] / num[valid_unc]) ** 2 + (derr[valid_unc] / den[valid_unc]) ** 2)
        ratio_unc[valid_unc] = ratio[valid_unc] * rel

    invalid_count = int(np.size(ratio) - np.count_nonzero(np.isfinite(ratio)))
    if invalid_count:
        LOG.warning("Masked %d invalid ratio pixels.", invalid_count)
    return ratio, ratio_unc


def write_ratio_outputs(result: RatioMapResult, outdir: str | Path, *, write_intermediate: bool = True) -> dict[str, Path]:
    """Write ratio FITS, plot, metadata, and optional intermediate feature maps."""

    outdir = ensure_dir(outdir)
    if write_intermediate:
        write_feature_map_outputs(result.numerator, outdir)
        write_feature_map_outputs(result.denominator, outdir)

    stem = f"{sanitize_name(result.feature1.feature_name)}_over_{sanitize_name(result.feature2.feature_name)}"
    paths: dict[str, Path] = {}
    paths["ratio"] = save_fits_product(result.ratio, result.header, outdir / f"{stem}_ratio.fits", bunit="ratio", metadata={"RATIO": stem})
    if result.ratio_uncertainty is not None:
        paths["ratio_uncertainty"] = save_fits_product(result.ratio_uncertainty, result.header, outdir / f"{stem}_ratio_unc.fits", bunit="ratio", metadata={"RATIO": stem})
    paths["figure"] = plot_ratio_map(result.ratio, outdir / f"{stem}_ratio.png", title=stem.replace("_", " "), header=result.header)
    metadata_path = outdir / f"{stem}_metadata.yaml"
    with metadata_path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(result.metadata, handle, sort_keys=False)
    paths["metadata"] = metadata_path
    result.output_paths.update(paths)
    return paths


def _ratio_metadata(num: FeatureMapResult, den: FeatureMapResult, reprojection_note: str) -> dict:
    return {
        "software": {"name": "JAFA", "full_name": "JWST Aromatic Feature Analyzer", "package": "jafa", "version": __version__},
        "ratio": {
            "numerator": num.feature.feature_name,
            "denominator": den.feature.feature_name,
            "reprojection": reprojection_note,
            "masking": "denominator must be finite and positive; SNR thresholds applied when uncertainties are available",
        },
        "numerator_metadata": num.metadata,
        "denominator_metadata": den.metadata,
    }
