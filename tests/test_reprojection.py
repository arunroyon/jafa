import numpy as np
import pytest
from astropy.io import fits

from jafa.reprojection import reproject_map, reproject_mask, reproject_uncertainty

pytest.importorskip("reproject")


def _header():
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 4
    header["NAXIS2"] = 4
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CRPIX1"] = 2.0
    header["CRPIX2"] = 2.0
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CDELT1"] = -1 / 3600
    header["CDELT2"] = 1 / 3600
    return header


def test_reproject_map_masks_by_footprint():
    data = np.arange(16, dtype=float).reshape(4, 4)
    header = _header()
    rep = reproject_map(data, header, header, data.shape)
    assert rep.data.shape == data.shape
    assert np.nanmedian(rep.footprint) > 0
    assert np.all(np.isfinite(rep.data[rep.footprint > 0]))


def test_reproject_uncertainty_uses_variance():
    uncertainty = np.full((4, 4), 2.0)
    header = _header()
    rep = reproject_uncertainty(uncertainty, header, header, uncertainty.shape)
    assert np.nanmedian(rep.data) == pytest.approx(2.0)


def test_reproject_mask_uses_nearest_neighbor_and_stays_boolean():
    mask = np.ones((4, 4), dtype=bool)
    mask[0, :] = False
    header = _header()

    rep = reproject_mask(mask, header, header, mask.shape)

    assert rep.data.dtype == bool
    assert np.array_equal(rep.data, mask)
