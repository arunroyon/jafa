import pytest

from jafa.feature_db import FeatureDefinition, custom_feature_definition, load_feature_database


def test_builtin_feature_database_contains_required_features():
    features = load_feature_database()
    for name in [
        "PAH_3p3",
        "PAH_6p2",
        "PAH_7p7",
        "PAH_8p6",
        "PAH_11p0",
        "PAH_11p2",
        "PAH_12p7",
        "C60_17p4",
        "PAH_16p4",
        "C60_18p9",
    ]:
        assert name in features

    c60 = features["C60_18p9"]
    assert c60.integration_window == (18.75, 19.15)
    assert c60.continuum_anchor_points[-1] == 20.5
    assert c60.continuum_mode == "morph_spline"
    assert c60.preferred_channel == "Channel 4"

    pah = features["PAH_16p4"]
    assert pah.integration_window == (16.2, 16.625)


def test_custom_feature_definition_from_cli_pairs():
    feature = custom_feature_definition(
        18.9,
        (18.82, 19.02),
        [18.55, 18.70, 19.15, 19.30],
    )
    assert feature.feature_name == "custom_18p900um"
    assert feature.continuum_anchors == ((18.55, 18.70), (19.15, 19.30))


def test_invalid_feature_definition_is_rejected():
    with pytest.raises(ValueError, match="zero width"):
        FeatureDefinition(
            feature_name="bad",
            central_wavelength=10.0,
            integration_window=(10.0, 10.0),
        )

    with pytest.raises(ValueError, match="Unsupported continuum_mode"):
        FeatureDefinition(
            feature_name="bad_mode",
            central_wavelength=10.0,
            integration_window=(9.9, 10.1),
            continuum_mode="not_a_mode",
        )


def test_required_window_includes_discrete_anchor_points():
    feature = FeatureDefinition(
        feature_name="anchored",
        central_wavelength=2.0,
        integration_window=(1.8, 2.2),
        continuum_anchor_points=(1.2, 2.8),
    )

    assert feature.required_window == (1.2, 2.8)


def test_package_data_loads_from_installed_resource():
    features = load_feature_database()
    assert features["C60_7p0"].continuum_anchor_points == (6.48, 6.81, 6.87, 7.25, 7.4, 7.55, 7.61)
