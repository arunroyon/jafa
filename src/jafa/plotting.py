"""Publication-oriented plotting helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.wcs import WCS

from .utils import require_dependency, robust_percentile_limits


def plot_feature_map(
    feature_map: np.ndarray,
    uncertainty: np.ndarray | None,
    output_path: str | Path,
    *,
    title: str,
    unit_label: str,
    diagnostic: dict[str, np.ndarray] | None = None,
    header=None,
    cmap_name: str = "rainbow",
) -> Path:
    """Save a feature-map diagnostic figure.

    Parameters
    ----------
    feature_map
        Integrated feature map.
    uncertainty
        Optional uncertainty map with the same shape as ``feature_map``.
    output_path
        Destination PNG/PDF path.
    title
        Map panel title.
    unit_label
        Colorbar unit label.
    diagnostic
        Optional continuum diagnostic payload produced by the mapper.
    header
        Optional two-dimensional WCS header.
    cmap_name
        Matplotlib colormap name.

    Returns
    -------
    pathlib.Path
        Written figure path.
    """

    plt = require_dependency("matplotlib.pyplot", extra="plot", purpose="plotting feature maps")
    output_path = Path(output_path)
    ncols = 3 if diagnostic else 2 if uncertainty is not None else 1
    fig = plt.figure(figsize=(5.2 * ncols, 4.8), constrained_layout=True)
    axes = []
    wcs = _wcs_from_header(header)
    for idx in range(ncols):
        if idx < (2 if uncertainty is not None else 1) and wcs is not None:
            axes.append(fig.add_subplot(1, ncols, idx + 1, projection=wcs))
        else:
            axes.append(fig.add_subplot(1, ncols, idx + 1))
    axes = np.asarray(axes)

    vmin, vmax = robust_percentile_limits(feature_map)
    im = axes[0].imshow(feature_map, origin="lower", interpolation="nearest", cmap=cmap_name, vmin=vmin, vmax=vmax)
    axes[0].set_title(title)
    _style_map_axis(axes[0])
    cbar = fig.colorbar(im, ax=axes[0], pad=0.02, fraction=0.046)
    cbar.set_label(unit_label, fontsize=11, fontweight="bold")
    cbar.ax.tick_params(labelsize=9)

    cursor = 1
    if uncertainty is not None:
        finite_unc = np.isfinite(uncertainty)
        if np.any(finite_unc):
            uvmin, uvmax = robust_percentile_limits(uncertainty)
        else:
            uvmin, uvmax = 0.0, 1.0
        im_unc = axes[cursor].imshow(uncertainty, origin="lower", interpolation="nearest", cmap=cmap_name, vmin=uvmin, vmax=uvmax)
        axes[cursor].set_title("Uncertainty")
        _style_map_axis(axes[cursor])
        cbar_unc = fig.colorbar(im_unc, ax=axes[cursor], pad=0.02, fraction=0.046)
        cbar_unc.set_label(unit_label, fontsize=11, fontweight="bold")
        cbar_unc.ax.tick_params(labelsize=9)
        cursor += 1

    if diagnostic:
        ax = axes[cursor]
        wl = diagnostic.get("wavelength_um")
        spectrum = diagnostic.get("spectrum")
        baseline = diagnostic.get("baseline")
        continuum = diagnostic.get("continuum")
        residual = diagnostic.get("residual")
        anchor_wavelengths = diagnostic.get("anchor_wavelengths")
        anchor_fluxes = diagnostic.get("anchor_fluxes")
        feature_window = diagnostic.get("feature_window")
        continuum_anchors = diagnostic.get("continuum_anchors")
        continuum_anchor_points = diagnostic.get("continuum_anchor_points")
        morph_half_window = diagnostic.get("morph_half_window")
        if wl is not None and spectrum is not None:
            ax.plot(wl, spectrum, color="0.15", lw=1.2, label="Median spectrum")
        if wl is not None and baseline is not None:
            ax.plot(wl, baseline, color="#2ca02c", lw=1.0, ls="--", label="Morph baseline")
        if wl is not None and continuum is not None:
            ax.plot(wl, continuum, color="#1f77b4", lw=1.4, label="Continuum")
        if wl is not None and residual is not None:
            ax.plot(wl, residual, color="#d62728", lw=1.0, label="Morph baseline - continuum")
        if continuum_anchors is not None:
            for lo, hi in continuum_anchors:
                ax.axvspan(lo, hi, color="#1f77b4", alpha=0.12)
        if continuum_anchor_points is not None:
            for point in continuum_anchor_points:
                ax.axvline(point, color="#1f77b4", alpha=0.35, lw=0.8, ls=":")
        if feature_window is not None:
            ax.axvspan(feature_window[0], feature_window[1], color="#d62728", alpha=0.10)
        if anchor_wavelengths is not None and anchor_fluxes is not None and len(anchor_wavelengths) > 0:
            ax.scatter(anchor_wavelengths, anchor_fluxes, marker="o", s=18, facecolor="white", edgecolor="#1f77b4", linewidth=0.9, label="Spline anchor samples", zorder=4)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_title("Continuum diagnostic")
        ax.set_xlabel("Wavelength (um)")
        ax.set_ylabel("Flux")
        if morph_half_window is not None and len(morph_half_window) > 0:
            ax.text(
                0.02,
                0.96,
                f"morph_half_window = {int(morph_half_window[0])}\nedge anchors = on",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.85, "boxstyle": "round,pad=0.25"},
            )
        ax.legend(loc="best", fontsize=8)

    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def plot_ratio_map(ratio: np.ndarray, output_path: str | Path, *, title: str, header=None, cmap_name: str = "rainbow") -> Path:
    """Save a ratio-map figure."""

    plt = require_dependency("matplotlib.pyplot", extra="plot", purpose="plotting ratio maps")
    output_path = Path(output_path)
    wcs = _wcs_from_header(header)
    fig = plt.figure(figsize=(8, 6), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=wcs) if wcs is not None else fig.add_subplot(1, 1, 1)
    vmin, vmax = robust_percentile_limits(ratio)
    im = ax.imshow(ratio, origin="lower", interpolation="nearest", cmap=cmap_name, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    _style_map_axis(ax)
    cbar = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label("Ratio", fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def _wcs_from_header(header):
    if header is None:
        return None
    try:
        return WCS(header)
    except Exception:
        return None


def _style_map_axis(ax) -> None:
    if hasattr(ax, "coords"):
        ax.coords.grid(color="0.85", ls=":", lw=0.8)
        ax.coords[0].set_axislabel("R.A (H:M:S)", fontsize=12, fontweight="bold")
        ax.coords[1].set_axislabel("Dec (D:M:S)", fontsize=12, fontweight="bold")
        ax.tick_params(axis="both", which="major", direction="in", labelsize=9, length=6, width=1.2)
    else:
        ax.set_xlabel("x pixel")
        ax.set_ylabel("y pixel")
