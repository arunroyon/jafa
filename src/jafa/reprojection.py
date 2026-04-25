"""Map reprojection utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.io import fits

from .utils import require_dependency


@dataclass(frozen=True)
class ReprojectedMap:
    """A reprojected image and its footprint.

    Parameters
    ----------
    data
        Reprojected map with pixels outside the reproject footprint set to
        NaN.
    footprint
        Reprojection footprint returned by ``reproject``.
    """

    data: np.ndarray
    footprint: np.ndarray


def reproject_map(
    data: np.ndarray,
    source_header: fits.Header,
    target_header: fits.Header,
    target_shape: tuple[int, int],
    *,
    method: str = "interp",
) -> ReprojectedMap:
    """Reproject a 2D map to a target WCS grid.

    Parameters
    ----------
    data
        Two-dimensional source map.
    source_header
        FITS/WCS header for ``data``.
    target_header
        FITS/WCS header defining the output celestial grid.
    target_shape
        Output shape as ``(ny, nx)``.
    method
        Reprojection method. Only ``"interp"`` is currently implemented.

    Returns
    -------
    ReprojectedMap
        Reprojected map and footprint.
    """

    if method != "interp":
        raise ValueError("Only interpolation reprojection is currently implemented.")
    reproject = require_dependency("reproject", extra="jwst", purpose="map reprojection")
    reproject_interp = reproject.reproject_interp
    result, footprint = reproject_interp((data, source_header), target_header, shape_out=target_shape)
    result = np.asarray(result, dtype=float)
    result[footprint <= 0] = np.nan
    return ReprojectedMap(data=result, footprint=np.asarray(footprint, dtype=float))


def reproject_uncertainty(
    uncertainty: np.ndarray,
    source_header: fits.Header,
    target_header: fits.Header,
    target_shape: tuple[int, int],
    *,
    method: str = "interp",
) -> ReprojectedMap:
    """Reproject an uncertainty map by interpolating variance and taking sqrt."""

    variance = np.asarray(uncertainty, dtype=float) ** 2
    rep = reproject_map(variance, source_header, target_header, target_shape, method=method)
    return ReprojectedMap(data=np.sqrt(rep.data), footprint=rep.footprint)
