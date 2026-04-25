"""Shared utility helpers."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from importlib import import_module
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .exceptions import OptionalDependencyError

LOGGER_NAME = "jafa"


def configure_logging(verbose: bool = False, quiet: bool = False) -> logging.Logger:
    """Configure package logging and return the package logger."""

    level = logging.WARNING if quiet else logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    return logger


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a child logger for this package."""

    if name:
        return logging.getLogger(f"{LOGGER_NAME}.{name}")
    return logging.getLogger(LOGGER_NAME)


def ensure_dir(path: str | Path) -> Path:
    """Create a directory if needed and return it as a ``Path``."""

    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def sanitize_name(name: str) -> str:
    """Return a filesystem-friendly feature or product name."""

    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return clean.strip("_") or "product"


def finite_fraction(values: np.ndarray) -> float:
    """Return the fraction of finite values in an array."""

    if values.size == 0:
        return 0.0
    return float(np.count_nonzero(np.isfinite(values)) / values.size)


def pairwise(values: Iterable[float]) -> list[tuple[float, float]]:
    """Convert an even-length sequence into wavelength pairs."""

    seq = [float(v) for v in values]
    if len(seq) % 2:
        raise ValueError("Anchor values must be supplied as wavelength pairs.")
    return [(seq[i], seq[i + 1]) for i in range(0, len(seq), 2)]


def require_dependency(module_name: str, *, extra: str | None = None, purpose: str | None = None):
    """Import an optional dependency or raise a helpful installation error.

    Parameters
    ----------
    module_name
        Importable module name.
    extra
        Optional package extra that provides the dependency.
    purpose
        Short description of why the dependency is needed.

    Returns
    -------
    module
        Imported Python module.

    Raises
    ------
    OptionalDependencyError
        If the dependency cannot be imported.
    """

    try:
        return import_module(module_name)
    except ImportError as exc:
        install = f"jafa[{extra}]" if extra else "the required optional extra"
        why = f" for {purpose}" if purpose else ""
        raise OptionalDependencyError(
            f"{module_name!r} is required{why}. Install it with `pip install {install}`."
        ) from exc


def robust_percentile_limits(values: np.ndarray, lower: float = 2, upper: float = 98) -> tuple[float, float]:
    """Return finite percentile limits with a fallback for nearly constant arrays."""

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("No finite values are available for plotting.")
    vmin, vmax = np.nanpercentile(finite, (lower, upper))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(finite))
        vmax = float(np.nanmax(finite))
        if vmin == vmax:
            eps = 1e-12 if vmin == 0 else abs(vmin) * 1e-6
            vmin -= eps
            vmax += eps
    return float(vmin), float(vmax)


def same_celestial_grid(
    header1: fits.Header,
    shape1: tuple[int, int],
    header2: fits.Header,
    shape2: tuple[int, int],
    *,
    tolerance_deg: float = 1e-7,
) -> bool:
    """Return whether two map products share the same celestial pixel grid.

    Empty or placeholder headers are treated conservatively as unknown, not as
    proof of a shared grid.
    """

    if tuple(shape1) != tuple(shape2):
        return False
    wcs1 = _validated_celestial_wcs(header1)
    wcs2 = _validated_celestial_wcs(header2)
    if wcs1 is None or wcs2 is None:
        return False
    if tuple(list(wcs1.wcs.ctype)[:2]) != tuple(list(wcs2.wcs.ctype)[:2]):
        return False

    ny, nx = int(shape1[0]), int(shape1[1])
    x = np.asarray([0.0, max(nx - 1, 0), 0.0, max(nx - 1, 0), 0.5 * max(nx - 1, 0)])
    y = np.asarray([0.0, 0.0, max(ny - 1, 0), max(ny - 1, 0), 0.5 * max(ny - 1, 0)])
    world1 = np.asarray(wcs1.pixel_to_world_values(x, y), dtype=float)
    world2 = np.asarray(wcs2.pixel_to_world_values(x, y), dtype=float)
    if world1.shape != world2.shape or not np.all(np.isfinite(world1)) or not np.all(np.isfinite(world2)):
        return False

    lon_delta = np.abs((world1[0] - world2[0] + 180.0) % 360.0 - 180.0)
    lat_delta = np.abs(world1[1] - world2[1])
    return bool(np.all(lon_delta <= tolerance_deg) and np.all(lat_delta <= tolerance_deg))


def _validated_celestial_wcs(header: fits.Header | None) -> WCS | None:
    if not header:
        return None
    required = ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2")
    if any(key not in header for key in required):
        return None
    try:
        wcs = WCS(header).celestial
    except Exception:
        return None
    if wcs.pixel_n_dim != 2 or wcs.world_n_dim != 2:
        return None
    if any(not item for item in list(wcs.wcs.ctype)[:2]):
        return None
    return wcs
