"""Legacy compatibility namespace for JAFA.

New code should import :mod:`jafa`. This module keeps the old
``jwst_feature_mapper`` top-level import working during the rename transition.
"""

from __future__ import annotations

from importlib import import_module

_jafa = import_module("jafa")
__path__ = _jafa.__path__
__all__ = list(getattr(_jafa, "__all__", ()))
__version__ = _jafa.__version__


def __getattr__(name: str):
    """Delegate legacy top-level attributes to :mod:`jafa`."""

    return getattr(_jafa, name)
