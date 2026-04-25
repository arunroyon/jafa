"""JAFA — JWST Aromatic Feature Analyzer.

Tools for making JWST spectral feature maps from IFU cubes.

The package keeps optional JWST, reprojection, and plotting dependencies out of
import time. Public objects are loaded lazily when first accessed.
"""

from importlib import import_module

from .version import __version__

_PUBLIC_IMPORTS = {
    "CubeData": ".io",
    "FeatureDefinition": ".feature_db",
    "FeatureMapResult": ".mapping",
    "FeatureMapper": ".mapping",
    "RatioMapResult": ".ratios",
    "load_cube": ".io",
    "load_feature_database": ".feature_db",
    "make_feature_map": ".mapping",
    "make_ratio_map": ".ratios",
    "save_fits_product": ".io",
    "select_cube_for_feature": ".cube_selection",
    "select_or_stitch_cube": ".mapping",
    "stitch_adjacent_cubes_for_feature": ".stitching",
    "stitch_cubes": ".stitching",
}

__all__ = [
    "CubeData",
    "FeatureDefinition",
    "FeatureMapResult",
    "FeatureMapper",
    "RatioMapResult",
    "__version__",
    "load_cube",
    "load_feature_database",
    "make_feature_map",
    "make_ratio_map",
    "save_fits_product",
    "select_cube_for_feature",
    "select_or_stitch_cube",
    "stitch_adjacent_cubes_for_feature",
    "stitch_cubes",
]


def __getattr__(name: str):
    """Load public API objects on first use."""

    if name in _PUBLIC_IMPORTS:
        module = import_module(_PUBLIC_IMPORTS[name], __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
