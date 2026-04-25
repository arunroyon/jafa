"""Command-line interface for JAFA."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from .exceptions import OptionalDependencyError
from .feature_db import FeatureDefinition, custom_feature_definition, load_feature_database
from .io import load_cube
from .mapping import MapSettings, make_feature_map
from .ratios import make_ratio_map
from .utils import configure_logging
from .version import __version__


def build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser."""

    parser = argparse.ArgumentParser(
        prog="jafa",
        description="JAFA — JWST Aromatic Feature Analyzer: create JWST IFU spectral feature and ratio maps.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    make_map = sub.add_parser("make-map", help="Create a continuum-subtracted feature map.")
    _add_common_args(make_map)
    make_map.add_argument("--feature", help="Feature name from the database.")
    make_map.add_argument("--wavelength", type=float, help="Custom central wavelength in microns.")
    make_map.add_argument("--feature-window", nargs=2, type=float, metavar=("MIN", "MAX"), help="Custom integration window in microns.")
    make_map.add_argument("--anchors", nargs="+", type=float, help="Custom continuum anchor pairs: lo1 hi1 lo2 hi2 ...")

    make_ratio = sub.add_parser("make-ratio", help="Create a ratio map from two named features.")
    _add_common_args(make_ratio)
    make_ratio.add_argument("--feature1", required=True, help="Numerator feature name.")
    make_ratio.add_argument("--feature2", required=True, help="Denominator feature name.")

    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the command-line interface."""

    parser = build_parser()
    args = parser.parse_args(argv)
    logger = configure_logging(verbose=args.verbose, quiet=args.quiet)
    try:
        config = _load_config(args.config)
        settings = _settings_from_args(args, config)
        feature_db = _feature_db_from_args(args.feature_db, config)

        if args.command == "make-map":
            feature = _feature_from_args(args, feature_db)
            cubes = [load_cube(path) for path in args.cubes]
            logger.info("Loaded %d cube(s).", len(cubes))
            result = make_feature_map(cubes, feature, outdir=args.outdir, settings=settings, write_outputs=True)
            logger.info("Selected cube: %s", result.cube.name)
            return 0

        if args.command == "make-ratio":
            feature1 = _feature_by_name(args.feature1, feature_db)
            feature2 = _feature_by_name(args.feature2, feature_db)
            cubes = [load_cube(path) for path in args.cubes]
            logger.info("Loaded %d cube(s).", len(cubes))
            result = make_ratio_map(cubes, feature1, feature2, outdir=args.outdir, settings=settings, write_outputs=True)
            logger.info("Numerator cube: %s", result.numerator.cube.name)
            logger.info("Denominator cube: %s", result.denominator.cube.name)
            return 0
    except (FileNotFoundError, KeyError, OSError, OptionalDependencyError, ValueError, yaml.YAMLError) as exc:
        if args.verbose:
            raise
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    parser.error(f"Unknown command: {args.command}")
    return 2


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cubes", nargs="+", required=True, help="Input FITS spectral cubes.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--config", help="Optional pipeline YAML configuration.")
    parser.add_argument("--feature-db", help="Optional YAML feature database to merge with built-ins.")
    parser.add_argument("--snr-threshold", type=float, help="Mask pixels below this SNR when uncertainties are available.")
    parser.add_argument("--morph-half-window", type=int, help="Morphology baseline half-window in spectral pixels.")
    parser.add_argument("--bin-spatial", type=int, default=None, help="Integer spatial binning factor.")
    parser.add_argument(
        "--output-unit",
        choices=("cgs", "native"),
        help="Integrated map unit mode: cgs flux per pixel or native wavelength-integrated cube units.",
    )
    parser.add_argument(
        "--continuum-mode",
        choices=("morph_spline", "spline"),
        help="Override feature continuum mode.",
    )
    parser.add_argument("--no-stitch-cubes", action="store_true", help="Disable automatic adjacent-cube stitching.")
    parser.add_argument(
        "--stitch-gap-tolerance",
        type=float,
        help="Maximum wavelength gap, in microns, allowed when stitching adjacent cubes.",
    )
    parser.add_argument(
        "--stitch-overlap",
        choices=("split", "keep"),
        help="How to handle spectral overlap between stitched cubes.",
    )
    parser.add_argument(
        "--allow-noncommon-stitch-footprint",
        action="store_true",
        help="Keep pixels that are only covered by one source cube in a stitched product.",
    )
    parser.add_argument("--no-clip-negative", action="store_true", help="Allow negative integrated residual fluxes.")
    parser.add_argument("--quiet", action="store_true", help="Show warnings and errors only.")
    parser.add_argument("--verbose", action="store_true", help="Show detailed log messages.")


def _load_config(path: str | None) -> dict:
    if not path:
        return {}
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _settings_from_args(args: argparse.Namespace, config: dict) -> MapSettings:
    settings_payload = config.get("settings") or {}
    settings = MapSettings.from_mapping(settings_payload)
    if args.snr_threshold is not None:
        settings.snr_threshold = args.snr_threshold
    if args.no_clip_negative:
        settings.clip_negative_residuals = False
    if args.bin_spatial is not None:
        settings.bin_spatial = args.bin_spatial
    if args.output_unit is not None:
        settings.output_unit = args.output_unit
    if args.no_stitch_cubes:
        settings.allow_cube_stitching = False
    if args.stitch_gap_tolerance is not None:
        settings.stitching.spectral_gap_tolerance_um = args.stitch_gap_tolerance
    if args.stitch_overlap is not None:
        settings.stitching.overlap_strategy = args.stitch_overlap
    if args.allow_noncommon_stitch_footprint:
        settings.stitching.require_common_spatial_footprint = False
    if args.continuum_mode is not None:
        settings.continuum.mode = args.continuum_mode
    if args.morph_half_window is not None:
        settings.continuum.morph_half_window = args.morph_half_window
    settings.progress = not args.quiet
    settings.validate()
    return settings


def _feature_from_args(args: argparse.Namespace, feature_db: dict):
    if args.feature:
        return _feature_by_name(args.feature, feature_db)
    if args.wavelength is None or args.feature_window is None or args.anchors is None:
        raise ValueError("Either --feature or all custom wavelength options are required.")
    return custom_feature_definition(args.wavelength, tuple(args.feature_window), args.anchors)


def _feature_by_name(name: str, feature_db: dict):
    if name not in feature_db:
        raise ValueError(f"Unknown feature {name!r}.")
    return feature_db[name]


def _feature_db_from_args(path: str | None, config: dict) -> dict:
    feature_db = load_feature_database(path)
    for name, payload in (config.get("features") or {}).items():
        feature_db[name] = FeatureDefinition.from_mapping(name, payload)
    return feature_db


if __name__ == "__main__":
    raise SystemExit(main())
