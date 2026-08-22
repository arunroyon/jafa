#!/usr/bin/env python3
"""Render the Tc 1 and NGC 7023 publication gallery from JAFA FITS products.

The input products are generated with the edge-mask-aware JAFA pipeline. This
renderer applies an additional signal-to-noise threshold for display only and
writes both 300 dpi PNG and vector PDF figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize
from astropy.wcs import WCS
from matplotlib.colors import Normalize

HD_200775 = SkyCoord("21h01m36.92s", "+68d09m47.8s")


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "font.weight": "normal",
            "axes.labelsize": 9,
            "axes.labelweight": "bold",
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
        }
    )


def _load_fits(path: Path) -> tuple[np.ndarray, fits.Header]:
    with fits.open(path, memmap=False) as hdul:
        return np.asarray(hdul[0].data, dtype=float), hdul[0].header.copy()


def _display_data(data: np.ndarray, uncertainty: np.ndarray, mask: np.ndarray, snr: float) -> np.ndarray:
    valid = np.asarray(mask, dtype=bool)
    valid &= np.isfinite(data) & np.isfinite(uncertainty) & (uncertainty > 0)
    valid &= data / uncertainty >= snr
    return np.where(valid, data, np.nan)


def _turbo_with_masked_background():
    cmap = plt.get_cmap("turbo").copy()
    cmap.set_bad("#eeeeee")
    return cmap


def _style_wcs_axis(ax, *, show_ylabel: bool = True) -> None:
    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_axislabel("Right ascension (J2000)", minpad=0.8, fontweight="bold")
    dec.set_axislabel("Declination (J2000)" if show_ylabel else "", minpad=0.8, fontweight="bold")
    ra.set_major_formatter("hh:mm:ss.s")
    dec.set_major_formatter("dd:mm:ss")
    ra.set_ticklabel(size=8, weight="bold", exclude_overlapping=True)
    dec.set_ticklabel(size=8, weight="bold", exclude_overlapping=True)
    if not show_ylabel:
        dec.set_ticklabel_visible(False)
    ax.coords.grid(color="white", alpha=0.42, linestyle=":", linewidth=0.7)
    ax.tick_params(axis="both", which="major", direction="in", length=5, width=1.0)


def _add_panel_label(ax, label: str) -> None:
    ax.text(
        0.035,
        0.955,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="white",
        path_effects=[path_effects.withStroke(linewidth=1.8, foreground="black")],
    )


def _save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", dpi=300)
    fig.savefig(output_dir / f"{stem}.pdf")
    plt.close(fig)


def render_tc1(tc1_dir: Path, output_dir: Path, *, snr: float) -> None:
    feature_map, header = _load_fits(tc1_dir / "C60_18p9_map.fits")
    uncertainty, _ = _load_fits(tc1_dir / "C60_18p9_unc.fits")
    mask, _ = _load_fits(tc1_dir / "C60_18p9_mask.fits")
    display_map = _display_data(feature_map, uncertainty, mask, snr)
    display_uncertainty = np.where(np.isfinite(display_map), uncertainty, np.nan)

    map_scale = 1.0e-14
    uncertainty_scale = 1.0e-17
    scaled_map = display_map / map_scale
    scaled_uncertainty = display_uncertainty / uncertainty_scale
    map_vmax = float(np.nanpercentile(scaled_map, 99.5))
    unc_vmax = float(np.nanpercentile(scaled_uncertainty, 99.5))
    map_norm = ImageNormalize(vmin=0.0, vmax=map_vmax, stretch=AsinhStretch(0.25))
    unc_norm = Normalize(vmin=0.0, vmax=unc_vmax)
    cmap = _turbo_with_masked_background()
    wcs = WCS(header).celestial

    fig = plt.figure(figsize=(7.2, 3.25), constrained_layout=True)
    axes = [fig.add_subplot(1, 2, index, projection=wcs) for index in (1, 2)]
    images = [
        axes[0].imshow(scaled_map, origin="lower", interpolation="nearest", cmap=cmap, norm=map_norm),
        axes[1].imshow(scaled_uncertainty, origin="lower", interpolation="nearest", cmap=cmap, norm=unc_norm),
    ]
    axes[0].set_title(r"C$_{60}$ 18.9 $\mu$m integrated flux")
    axes[1].set_title(r"Statistical uncertainty ($1\sigma$)")
    _style_wcs_axis(axes[0])
    _style_wcs_axis(axes[1], show_ylabel=False)
    _add_panel_label(axes[0], "(a)")
    _add_panel_label(axes[1], "(b)")

    units = (
        r"$10^{-14}$ erg s$^{-1}$ cm$^{-2}$ pixel$^{-1}$",
        r"$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ pixel$^{-1}$",
    )
    for ax, image, unit in zip(axes, images, units, strict=True):
        colorbar = fig.colorbar(image, ax=ax, pad=0.025, fraction=0.048)
        colorbar.set_label(unit, fontsize=7.5, fontweight="bold", labelpad=3)
        colorbar.ax.tick_params(labelsize=7.5, width=1.0)
        for tick in colorbar.ax.get_yticklabels():
            tick.set_fontweight("bold")

    fig.suptitle(r"Tc 1: spatially resolved C$_{60}$ emission", fontsize=11, fontweight="bold")
    _save_figure(fig, output_dir, "tc1_c60_18p9_diagnostic")


def _add_star_direction(ax, wcs: WCS, star: SkyCoord) -> None:
    ny, nx = ax.images[0].get_array().shape
    center = wcs.pixel_to_world((nx - 1) / 2.0, (ny - 1) / 2.0)
    position_angle = center.position_angle(star).to_value(u.rad)
    dx = -np.sin(position_angle)
    dy = np.cos(position_angle)
    origin = np.array([0.86, 0.88])
    endpoint = origin + 0.115 * np.array([dx, dy])
    ax.annotate(
        "",
        xy=endpoint,
        xytext=origin,
        xycoords="axes fraction",
        arrowprops={
            "arrowstyle": "-|>",
            "facecolor": "white",
            "edgecolor": "black",
            "lw": 1.2,
            "mutation_scale": 13,
        },
    )
    ax.text(
        0.965,
        0.955,
        "toward HD 200775",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        fontweight="bold",
        color="white",
        path_effects=[path_effects.withStroke(linewidth=2.0, foreground="black")],
    )


def render_ngc7023(ngc7023_dir: Path, output_dir: Path, *, snr: float) -> None:
    stem = "PAH_11p0_over_PAH_11p2"
    ratio, header = _load_fits(ngc7023_dir / f"{stem}_ratio.fits")
    uncertainty, _ = _load_fits(ngc7023_dir / f"{stem}_ratio_unc.fits")
    mask, _ = _load_fits(ngc7023_dir / f"{stem}_mask.fits")
    display_ratio = _display_data(ratio, uncertainty, mask, snr)
    vmin, vmax = np.nanpercentile(display_ratio, [2.0, 98.0])
    norm = Normalize(vmin=max(0.0, float(vmin)), vmax=float(vmax))
    cmap = _turbo_with_masked_background()
    wcs = WCS(header).celestial

    fig = plt.figure(figsize=(4.35, 5.0), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    image = ax.imshow(display_ratio, origin="lower", interpolation="nearest", cmap=cmap, norm=norm)
    _style_wcs_axis(ax)
    _add_panel_label(ax, "NGC 7023 NW")
    _add_star_direction(ax, wcs, HD_200775)
    ax.text(
        0.04,
        0.045,
        rf"S/N $\geq$ {snr:g}",
        transform=ax.transAxes,
        fontsize=7.5,
        fontweight="bold",
        color="white",
        path_effects=[path_effects.withStroke(linewidth=2.0, foreground="black")],
    )
    ax.set_title(r"PAH ionization proxy: 11.0 $\mu$m / 11.2 $\mu$m", fontsize=10.5)
    colorbar = fig.colorbar(image, ax=ax, pad=0.025, fraction=0.05)
    colorbar.set_label(r"PAH 11.0 $\mu$m / 11.2 $\mu$m", fontsize=9, fontweight="bold")
    colorbar.ax.tick_params(labelsize=8, width=1.0)
    for tick in colorbar.ax.get_yticklabels():
        tick.set_fontweight("bold")
    _save_figure(fig, output_dir, "ngc7023_pah_11p0_over_11p2_ratio")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tc1-dir", type=Path, default=Path("outputs/publication_gallery/tc1"))
    parser.add_argument("--ngc7023-dir", type=Path, default=Path("outputs/publication_gallery/ngc7023"))
    parser.add_argument("--output-dir", type=Path, default=Path("examples/gallery"))
    parser.add_argument("--snr", type=float, default=3.0, help="Display-only minimum signal-to-noise ratio.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    _configure_style()
    render_tc1(args.tc1_dir, args.output_dir, snr=args.snr)
    render_ngc7023(args.ngc7023_dir, args.output_dir, snr=args.snr)


if __name__ == "__main__":
    main()
