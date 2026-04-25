"""Sphinx configuration for JAFA."""

from __future__ import annotations

project = "JAFA — JWST Aromatic Feature Analyzer"
author = "Arun Roy"
copyright = "2026, Arun Roy"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "alabaster"
