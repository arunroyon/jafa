import numpy as np

from jafa.cube_selection import choose_coarser_cube, select_cube_for_feature
from jafa.feature_db import FeatureDefinition
from jafa.io import CubeData


def _cube(name, wavelengths, shape=(2, 2)):
    data = np.zeros((len(wavelengths), *shape), dtype=float)
    return CubeData(path=None, data=data, wavelength_um=np.asarray(wavelengths), channel=name)


def test_select_cube_for_feature_requires_anchors():
    feature = FeatureDefinition(
        feature_name="X",
        central_wavelength=10.0,
        integration_window=(9.9, 10.1),
        continuum_anchors=((9.5, 9.7), (10.3, 10.5)),
    )
    short = _cube("short", np.linspace(9.8, 10.2, 20))
    full = _cube("full", np.linspace(9.4, 10.6, 20))

    selected = select_cube_for_feature([short, full], feature)
    assert selected.cube is full


def test_choose_coarser_cube_uses_smaller_shape_when_scale_unknown():
    high_res = _cube("high", np.linspace(1, 2, 4), shape=(20, 20))
    coarse = _cube("coarse", np.linspace(1, 2, 4), shape=(10, 10))
    assert choose_coarser_cube([high_res, coarse]) is coarse
