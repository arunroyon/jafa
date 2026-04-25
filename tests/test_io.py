import numpy as np
from astropy.io import fits

from jafa.io import binned_spatial_header, load_cube, save_fits_product


def _simple_cube_header():
    header = fits.Header()
    header["BUNIT"] = "MJy/sr"
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CTYPE3"] = "WAVE"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CUNIT3"] = "um"
    header["CRPIX1"] = 1.0
    header["CRPIX2"] = 1.0
    header["CRPIX3"] = 1.0
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CRVAL3"] = 5.0
    header["CDELT1"] = -1 / 3600
    header["CDELT2"] = 1 / 3600
    header["CDELT3"] = 0.1
    return header


def test_load_cube_astropy_fallback_reads_wavelength_quantity(tmp_path):
    path = tmp_path / "cube.fits"
    data = np.ones((4, 3, 2), dtype=float)
    fits.PrimaryHDU(data=data, header=_simple_cube_header()).writeto(path)

    cube = load_cube(path, use_spectral_cube=False)

    assert cube.data.shape == (4, 3, 2)
    assert np.allclose(cube.wavelength_um, [5.0, 5.1, 5.2, 5.3])
    assert str(cube.wavelength.unit) == "micron"
    assert cube.flux_unit_object is not None


def test_save_fits_product_and_binned_header(tmp_path):
    header = fits.Header()
    header["CRPIX1"] = 3.0
    header["CRPIX2"] = 5.0
    header["CDELT1"] = -0.1
    header["CDELT2"] = 0.1
    binned = binned_spatial_header(header, 2)
    assert binned["CDELT1"] == -0.2
    assert binned["CDELT2"] == 0.2

    path = save_fits_product(np.ones((2, 2)), binned, tmp_path / "map.fits", bunit="MJy um")
    with fits.open(path) as hdul:
        assert hdul[0].header["BUNIT"] == "MJy um"
        assert hdul[0].data.shape == (2, 2)

