"""Package-specific exceptions."""

from __future__ import annotations


class JWSTFeatureMapperError(Exception):
    """Base exception for ``jafa`` errors."""


class OptionalDependencyError(ImportError, JWSTFeatureMapperError):
    """Raised when an optional dependency is required but unavailable."""

