"""FITS cube I/O and FITS product writing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area, proj_plane_pixel_scales

from .utils import get_logger, require_dependency

LOG = get_logger(__name__)

JWST_DO_NOT_USE = np.uint32(1)
JWST_NON_SCIENCE = np.uint32(512)
DEFAULT_BAD_DQ_BITS = JWST_DO_NOT_USE | JWST_NON_SCIENCE


@dataclass
class CubeData:
    """Loaded spectral cube with spectral axis first: ``(wavelength, y, x)``.

    Parameters
    ----------
    path
        Original FITS file path, or `None` for synthetic/in-memory cubes.
    data
        Flux data ordered as ``(spectral, y, x)``.
    wavelength_um
        Spectral axis in microns. A `~astropy.units.Quantity` is accepted and
        converted to microns.
    uncertainty
        Optional uncertainty cube ordered like ``data``.
    dq
        Optional bit-encoded data-quality cube ordered like ``data``.
    weight
        Optional cube-build weight map (the JWST ``WMAP`` extension) ordered
        like ``data``.
    wcs_3d, wcs_2d, header
        Full cube WCS, celestial WCS, and science FITS header.
    flux_unit, uncertainty_unit
        FITS unit strings.
    pixel_area_sr
        Spatial pixel solid angle in steradians, when available.
    instrument, channel
        Instrument metadata used as cube-selection hints.
    name_hint
        Optional display name for in-memory products such as stitched cubes.
    stitched_from
        Names of source cubes used to construct this cube.
    stitch_target
        Source cube whose spatial grid defines a stitched product.
    stitch_overlap_strategy
        Overlap handling used for a stitched product.
    stitch_overlap_ranges_um
        Wavelength overlap ranges encountered while stitching.
    stitch_common_footprint_fraction
        Fraction of target-grid pixels retained after requiring common spatial
        coverage from every stitched source cube.
    """

    path: Path | None
    data: np.ndarray
    wavelength_um: np.ndarray
    uncertainty: np.ndarray | None = None
    dq: np.ndarray | None = None
    weight: np.ndarray | None = None
    wcs_3d: WCS | None = None
    wcs_2d: WCS | None = None
    header: fits.Header | None = None
    flux_unit: str = ""
    uncertainty_unit: str = ""
    pixel_area_sr: float | None = None
    instrument: str | None = None
    channel: str | None = None
    name_hint: str | None = None
    stitched_from: tuple[str, ...] = ()
    stitch_target: str | None = None
    stitch_overlap_strategy: str | None = None
    stitch_overlap_ranges_um: tuple[tuple[float, float], ...] = ()
    stitch_common_footprint_fraction: float | None = None

    def __post_init__(self) -> None:
        if isinstance(self.wavelength_um, u.Quantity):
            self.wavelength_um = np.asarray(self.wavelength_um.to_value(u.micron), dtype=float)
        else:
            self.wavelength_um = np.asarray(self.wavelength_um, dtype=float)
        self.data = np.asarray(self.data, dtype=float)
        if self.uncertainty is not None:
            self.uncertainty = np.asarray(self.uncertainty, dtype=float)
        if self.dq is not None:
            self.dq = np.asarray(self.dq, dtype=np.uint32)
        if self.weight is not None:
            self.weight = np.asarray(self.weight, dtype=float)
        if self.data.ndim != 3:
            raise ValueError("Cube data must have shape (spectral, y, x).")
        if self.wavelength_um.size != self.data.shape[0]:
            raise ValueError("Wavelength axis length must match the cube spectral dimension.")
        for name, values in (("uncertainty", self.uncertainty), ("dq", self.dq), ("weight", self.weight)):
            if values is not None and values.shape != self.data.shape:
                raise ValueError(f"Cube {name} shape {values.shape} does not match science shape {self.data.shape}.")

    def valid_voxel_mask(self, *, bad_dq_bits: int = int(DEFAULT_BAD_DQ_BITS)) -> np.ndarray:
        """Return voxels with finite science data and usable cube-build quality."""

        valid = np.isfinite(self.data)
        if self.dq is not None:
            valid &= (self.dq & np.uint32(bad_dq_bits)) == 0
        if self.weight is not None:
            valid &= np.isfinite(self.weight) & (self.weight > 0)
        return valid

    @property
    def wavelength(self) -> u.Quantity:
        """Return the spectral axis as a micron Quantity."""

        return self.wavelength_um * u.micron

    @property
    def flux_unit_object(self) -> u.UnitBase | None:
        """Return ``flux_unit`` as an Astropy unit when parseable."""

        return _unit_or_none(self.flux_unit)

    @property
    def uncertainty_unit_object(self) -> u.UnitBase | None:
        """Return ``uncertainty_unit`` as an Astropy unit when parseable."""

        return _unit_or_none(self.uncertainty_unit)

    @property
    def name(self) -> str:
        """Return a human-readable cube name."""

        if self.name_hint:
            return self.name_hint
        return self.path.name if self.path else "in_memory_cube"

    @property
    def spectral_range(self) -> tuple[float, float]:
        """Return the wavelength range in microns."""

        finite = self.wavelength_um[np.isfinite(self.wavelength_um)]
        if finite.size == 0:
            return (np.nan, np.nan)
        return float(np.nanmin(finite)), float(np.nanmax(finite))

    @property
    def spatial_shape(self) -> tuple[int, int]:
        """Return spatial shape as ``(ny, nx)``."""

        return int(self.data.shape[1]), int(self.data.shape[2])

    @property
    def spatial_header(self) -> fits.Header:
        """Return a 2D WCS header for map products."""

        if self.wcs_2d is not None:
            return self.wcs_2d.to_header()
        if self.header is not None:
            try:
                return WCS(self.header).celestial.to_header()
            except Exception:
                return fits.Header()
        return fits.Header()

    @property
    def pixel_scale_arcsec(self) -> float | None:
        """Return an approximate spatial pixel scale in arcsec."""

        wcs = self.wcs_2d
        if wcs is None and self.header is not None:
            try:
                wcs = WCS(self.header).celestial
            except Exception:
                wcs = None
        if wcs is None:
            return None
        try:
            scales = proj_plane_pixel_scales(wcs) * u.deg
            return float(np.nanmean(np.abs(scales.to_value(u.arcsec))))
        except Exception:
            return None


def load_cube(path: str | Path, *, use_spectral_cube: bool = True) -> CubeData:
    """Load a JWST-style IFU cube from FITS.

    The returned data are always sorted by increasing wavelength and arranged as
    ``(spectral, y, x)``. JWST ``SCI``, ``ERR``, ``DQ``, and ``WMAP``
    extensions are preserved when present, while the loader falls back to the
    first 3D image extension for simple cubes.

    Parameters
    ----------
    path
        Input FITS cube path.
    use_spectral_cube
        Use :mod:`spectral_cube` for robust cube/WCS interpretation. When
        `False`, a pure-Astropy fallback is used for simple FITS cubes.

    Returns
    -------
    CubeData
        Loaded cube data and metadata.
    """

    path = Path(path)
    with fits.open(path, memmap=False) as hdul:
        sci_index = _find_image_hdu(hdul, preferred_names=("SCI", "FLUX", "DATA"))
        sci_hdu = hdul[sci_index]
        header = sci_hdu.header.copy()
        wcs_3d = WCS(header)

        if use_spectral_cube:
            spectral_cube = require_dependency("spectral_cube", extra="jwst", purpose="loading JWST spectral cubes")
            cube = spectral_cube.SpectralCube.read(sci_hdu)
            data = np.asarray(cube.unmasked_data[:].value, dtype=float)
            wavelength_um = cube.spectral_axis.to(u.micron).value.astype(float)
        else:
            data, wavelength_um = _read_cube_astropy(sci_hdu, wcs_3d)

        err_index = _find_image_hdu(hdul, preferred_names=("ERR", "ERROR", "UNCERTAINTY"), ndim=3, required=False)
        uncertainty = None
        uncertainty_unit = ""
        if err_index is not None:
            uncertainty = _match_shape(np.asarray(hdul[err_index].data, dtype=float), data.shape)
            uncertainty_unit = str(hdul[err_index].header.get("BUNIT", header.get("BUNIT", "")))

        dq_index = _find_image_hdu(hdul, preferred_names=("DQ",), ndim=3, required=False)
        dq = None
        if dq_index is not None:
            dq = _match_shape(np.asarray(hdul[dq_index].data, dtype=np.uint32), data.shape)

        weight_index = _find_image_hdu(hdul, preferred_names=("WMAP", "WEIGHT", "WHT"), ndim=3, required=False)
        weight = None
        if weight_index is not None:
            weight = _match_shape(np.asarray(hdul[weight_index].data, dtype=float), data.shape)

        if not np.all(np.diff(wavelength_um) > 0):
            order = np.argsort(wavelength_um)
            wavelength_um = wavelength_um[order]
            data = data[order, :, :]
            if uncertainty is not None:
                uncertainty = uncertainty[order, :, :]
            if dq is not None:
                dq = dq[order, :, :]
            if weight is not None:
                weight = weight[order, :, :]

        pixel_area_sr = _pixel_area_sr(header, wcs_3d)

        return CubeData(
            path=path,
            data=data,
            uncertainty=uncertainty,
            dq=dq,
            weight=weight,
            wavelength_um=wavelength_um,
            wcs_3d=wcs_3d,
            wcs_2d=wcs_3d.celestial,
            header=header,
            flux_unit=str(header.get("BUNIT", "")),
            uncertainty_unit=uncertainty_unit,
            pixel_area_sr=pixel_area_sr,
            instrument=header.get("INSTRUME"),
            channel=str(header.get("CHANNEL", header.get("BAND", "")) or ""),
        )


def save_fits_product(
    data: np.ndarray,
    header: fits.Header,
    output_path: str | Path,
    *,
    bunit: str | None = None,
    overwrite: bool = True,
    metadata: dict[str, Any] | None = None,
) -> Path:
    """Save a 2D FITS product with a WCS header and lightweight metadata."""

    output_path = Path(output_path)
    out_header = header.copy()
    if bunit:
        out_header["BUNIT"] = bunit
    if metadata:
        for key, value in metadata.items():
            fits_key = str(key).upper()[:8]
            if fits_key in out_header:
                continue
            try:
                out_header[fits_key] = value
            except Exception:
                continue
    fits.PrimaryHDU(data=np.asarray(data), header=out_header).writeto(output_path, overwrite=overwrite)
    LOG.info("Saved FITS product: %s", output_path)
    return output_path


def binned_spatial_header(header: fits.Header, bin_spatial: int) -> fits.Header:
    """Return a spatial header adjusted for integer spatial binning."""

    out = header.copy()
    if bin_spatial <= 1:
        return out
    if "CRPIX1" in out:
        out["CRPIX1"] = (out["CRPIX1"] + 0.5 * (bin_spatial - 1)) / bin_spatial
    if "CRPIX2" in out:
        out["CRPIX2"] = (out["CRPIX2"] + 0.5 * (bin_spatial - 1)) / bin_spatial
    if "CDELT1" in out:
        out["CDELT1"] *= bin_spatial
    if "CDELT2" in out:
        out["CDELT2"] *= bin_spatial
    if "CD1_1" in out:
        out["CD1_1"] *= bin_spatial
    if "CD1_2" in out:
        out["CD1_2"] *= bin_spatial
    if "CD2_1" in out:
        out["CD2_1"] *= bin_spatial
    if "CD2_2" in out:
        out["CD2_2"] *= bin_spatial
    return out


def _find_image_hdu(
    hdul: fits.HDUList,
    *,
    preferred_names: tuple[str, ...],
    ndim: int = 3,
    required: bool = True,
) -> int | None:
    preferred = {name.upper() for name in preferred_names}
    candidates: list[int] = []
    for idx, hdu in enumerate(hdul):
        data = getattr(hdu, "data", None)
        if data is None or np.ndim(data) != ndim:
            continue
        candidates.append(idx)
        extname = str(hdu.header.get("EXTNAME", "")).upper()
        if extname in preferred:
            return idx
    if candidates and required:
        return candidates[0]
    if required:
        raise ValueError(f"No {ndim}D image extension found in FITS file.")
    return None


def _match_shape(values: np.ndarray, target_shape: tuple[int, ...]) -> np.ndarray:
    if values.shape == target_shape:
        return values
    for axes in ((2, 0, 1), (1, 2, 0), (0, 2, 1), (2, 1, 0), (1, 0, 2)):
        if values.ndim == 3 and np.transpose(values, axes).shape == target_shape:
            return np.transpose(values, axes)
    raise ValueError(f"Uncertainty cube shape {values.shape} cannot be aligned to science shape {target_shape}.")


def _pixel_area_sr(header: fits.Header, wcs_3d: WCS) -> float | None:
    if "PIXAR_SR" in header:
        try:
            return float(header["PIXAR_SR"])
        except Exception:
            pass
    try:
        area_deg2 = abs(float(proj_plane_pixel_area(wcs_3d.celestial)))
        return float((area_deg2 * (u.deg**2)).to_value(u.sr))
    except Exception:
        LOG.warning("Could not determine pixel area; map units will omit pixel-area scaling.")
        return None


def _read_cube_astropy(sci_hdu: fits.ImageHDU | fits.PrimaryHDU, wcs_3d: WCS) -> tuple[np.ndarray, np.ndarray]:
    """Read a simple cube using Astropy WCS only."""

    data = np.asarray(sci_hdu.data, dtype=float)
    if data.ndim != 3:
        raise ValueError("Astropy fallback loader requires a 3D image cube.")
    spectral_numpy_axis = _spectral_numpy_axis(wcs_3d, data.ndim)
    if spectral_numpy_axis != 0:
        data = np.moveaxis(data, spectral_numpy_axis, 0)
    wavelength_um = _spectral_axis_from_wcs(wcs_3d, data.shape[0])
    return data, wavelength_um


def _spectral_numpy_axis(wcs: WCS, ndim: int) -> int:
    physical_types = [str(item or "").lower() for item in wcs.world_axis_physical_types]
    spectral_world_axis = None
    for idx, physical_type in enumerate(physical_types):
        if "em." in physical_type or "spect" in physical_type or "freq" in physical_type:
            spectral_world_axis = idx
            break
    if spectral_world_axis is None:
        ctype = [str(item).upper() for item in wcs.wcs.ctype]
        for idx, item in enumerate(ctype):
            if any(token in item for token in ("WAVE", "FREQ", "VRAD", "VELO")):
                spectral_world_axis = idx
                break
    if spectral_world_axis is None:
        return 0
    correlated = np.where(wcs.axis_correlation_matrix[spectral_world_axis])[0]
    if correlated.size != 1:
        return 0
    fits_pixel_axis = int(correlated[0])
    return ndim - 1 - fits_pixel_axis


def _spectral_axis_from_wcs(wcs: WCS, nspec: int) -> np.ndarray:
    spectral_numpy_axis = _spectral_numpy_axis(wcs, 3)
    fits_pixel_axis = 3 - 1 - spectral_numpy_axis
    pixels = [np.zeros(nspec), np.zeros(nspec), np.zeros(nspec)]
    pixels[fits_pixel_axis] = np.arange(nspec, dtype=float)
    world_values = wcs.pixel_to_world_values(*pixels)
    spectral_world_axis = None
    for idx, physical_type in enumerate(wcs.world_axis_physical_types):
        text = str(physical_type or "").lower()
        if "em." in text or "spect" in text or "freq" in text:
            spectral_world_axis = idx
            break
    if spectral_world_axis is None:
        spectral_world_axis = fits_pixel_axis
    values = np.asarray(world_values[spectral_world_axis], dtype=float)
    unit_name = wcs.world_axis_units[spectral_world_axis] or "m"
    return (values * u.Unit(unit_name)).to_value(u.micron)


def _unit_or_none(unit_text: str | None) -> u.UnitBase | None:
    if not unit_text:
        return None
    try:
        return u.Unit(str(unit_text))
    except Exception:
        return None
