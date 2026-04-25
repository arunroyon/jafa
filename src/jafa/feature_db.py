"""Feature database loading and validation."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from .utils import pairwise


@dataclass(frozen=True)
class FeatureDefinition:
    """Definition of a spectral feature and its continuum anchors.

    Parameters
    ----------
    feature_name
        Stable feature identifier, for example ``"C60_18p9"``.
    central_wavelength
        Approximate central wavelength in microns.
    integration_window
        Integration bounds in microns.
    continuum_mode
        Continuum mode used by the mapper unless runtime settings explicitly
        override it.
    continuum_anchors
        Anchor wavelength windows in microns.
    continuum_anchor_points
        Discrete anchor wavelengths in microns.
    notes
        Optional scientific notes.
    preferred_instrument, preferred_channel
        Optional cube-selection hints.
    """

    feature_name: str
    central_wavelength: float
    integration_window: tuple[float, float]
    continuum_mode: str = "morph_spline"
    continuum_anchors: tuple[tuple[float, float], ...] = field(default_factory=tuple)
    continuum_anchor_points: tuple[float, ...] = field(default_factory=tuple)
    notes: str | None = None
    preferred_instrument: str | None = None
    preferred_channel: str | None = None

    def __post_init__(self) -> None:
        if not str(self.feature_name).strip():
            raise ValueError("feature_name must be non-empty.")
        if not np.isfinite(float(self.central_wavelength)):
            raise ValueError(f"central_wavelength must be finite for {self.feature_name!r}.")
        object.__setattr__(self, "central_wavelength", float(self.central_wavelength))
        object.__setattr__(self, "integration_window", _as_pair(self.integration_window, "integration_window"))
        object.__setattr__(
            self,
            "continuum_anchors",
            tuple(_as_pair(item, "continuum_anchors") for item in self.continuum_anchors),
        )
        anchor_points = tuple(float(item) for item in self.continuum_anchor_points)
        if any(not np.isfinite(item) for item in anchor_points):
            raise ValueError(f"continuum_anchor_points must be finite for {self.feature_name!r}.")
        object.__setattr__(self, "continuum_anchor_points", anchor_points)
        if self.continuum_mode not in {"spline", "morph_spline"}:
            raise ValueError(f"Unsupported continuum_mode for {self.feature_name!r}: {self.continuum_mode!r}.")

    @property
    def required_window(self) -> tuple[float, float]:
        """Return the full wavelength span needed for feature plus anchors."""

        lows = [self.integration_window[0]]
        highs = [self.integration_window[1]]
        for lo, hi in self.continuum_anchors:
            lows.append(min(lo, hi))
            highs.append(max(lo, hi))
        if self.continuum_anchor_points:
            lows.extend(self.continuum_anchor_points)
            highs.extend(self.continuum_anchor_points)
        return min(lows), max(highs)

    @classmethod
    def from_mapping(cls, name: str, payload: Mapping[str, Any]) -> FeatureDefinition:
        """Build a feature definition from YAML/dict data."""

        feature_name = str(payload.get("feature_name", name))
        window = _as_pair(payload["integration_window"], "integration_window")
        anchors = tuple(_as_pair(item, "continuum_anchors") for item in payload.get("continuum_anchors", ()))
        anchor_points = tuple(float(item) for item in payload.get("continuum_anchor_points", ()))
        return cls(
            feature_name=feature_name,
            central_wavelength=float(payload["central_wavelength"]),
            integration_window=window,
            continuum_mode=str(payload.get("continuum_mode", "morph_spline")),
            continuum_anchors=anchors,
            continuum_anchor_points=anchor_points,
            notes=payload.get("notes"),
            preferred_instrument=payload.get("preferred_instrument"),
            preferred_channel=payload.get("preferred_channel"),
        )

    def to_metadata(self) -> dict[str, Any]:
        """Return a YAML-serializable representation."""

        payload: dict[str, Any] = {
            "feature_name": self.feature_name,
            "central_wavelength": self.central_wavelength,
            "integration_window": list(self.integration_window),
            "continuum_mode": self.continuum_mode,
            "continuum_anchors": [list(a) for a in self.continuum_anchors],
        }
        if self.continuum_anchor_points:
            payload["continuum_anchor_points"] = list(self.continuum_anchor_points)
        if self.notes:
            payload["notes"] = self.notes
        if self.preferred_instrument:
            payload["preferred_instrument"] = self.preferred_instrument
        if self.preferred_channel:
            payload["preferred_channel"] = self.preferred_channel
        return payload


def _as_pair(value: Any, field_name: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError(f"{field_name} entries must be [min, max] pairs.")
    lo, hi = float(value[0]), float(value[1])
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise ValueError(f"{field_name} pair must contain finite values: {value!r}")
    if lo == hi:
        raise ValueError(f"{field_name} pair has zero width: {value!r}")
    return (min(lo, hi), max(lo, hi))


def load_feature_database(path: str | Path | None = None) -> dict[str, FeatureDefinition]:
    """Load the built-in feature database, optionally merged with a user YAML file.

    User entries with the same key replace built-in entries.
    """

    with resources.files("jafa.data").joinpath("features.yaml").open("r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    if path is not None:
        with Path(path).open("r", encoding="utf-8") as handle:
            custom = yaml.safe_load(handle) or {}
        raw.update(custom)

    return {name: FeatureDefinition.from_mapping(name, payload) for name, payload in raw.items()}


def custom_feature_definition(
    wavelength: float,
    feature_window: tuple[float, float],
    anchors: list[float] | tuple[float, ...],
    name: str | None = None,
) -> FeatureDefinition:
    """Create a feature definition from CLI wavelength/window/anchor values."""

    anchor_pairs = tuple(pairwise(anchors))
    return FeatureDefinition(
        feature_name=name or f"custom_{float(wavelength):.3f}um".replace(".", "p"),
        central_wavelength=float(wavelength),
        integration_window=(min(feature_window), max(feature_window)),
        continuum_mode="morph_spline",
        continuum_anchors=anchor_pairs,
        notes="User-supplied custom wavelength feature.",
    )
