import pytest

from jafa.cli import _settings_from_args, build_parser, main


def test_make_map_parser_accepts_custom_feature_args():
    parser = build_parser()
    args = parser.parse_args(
        [
            "make-map",
            "--cubes",
            "cube.fits",
            "--wavelength",
            "18.9",
            "--feature-window",
            "18.75",
            "19.15",
            "--anchors",
            "18.0",
            "18.2",
            "19.2",
            "19.5",
            "--outdir",
            "results",
            "--output-unit",
            "native",
            "--continuum-mode",
            "spline",
            "--no-stitch-cubes",
            "--stitch-overlap",
            "keep",
            "--allow-noncommon-stitch-footprint",
        ]
    )
    assert args.command == "make-map"
    assert args.wavelength == 18.9
    assert args.feature_window == [18.75, 19.15]
    assert args.output_unit == "native"
    assert args.continuum_mode == "spline"
    assert args.no_stitch_cubes is True
    assert args.stitch_overlap == "keep"
    assert args.allow_noncommon_stitch_footprint is True


def test_make_ratio_parser_accepts_named_features():
    parser = build_parser()
    args = parser.parse_args(
        [
            "make-ratio",
            "--cubes",
            "cube1.fits",
            "cube2.fits",
            "--feature1",
            "C60_7p0",
            "--feature2",
            "C60_18p9",
            "--outdir",
            "results",
            "--quiet",
        ]
    )
    assert args.command == "make-ratio"
    assert args.quiet is True


def test_cli_help_uses_jafa_branding(capsys):
    parser = build_parser()

    with pytest.raises(SystemExit) as excinfo:
        parser.parse_args(["--help"])

    captured = capsys.readouterr()
    assert excinfo.value.code == 0
    assert "jafa" in captured.out
    assert "JAFA" in captured.out
    assert "JWST Aromatic Feature Analyzer" in captured.out


def test_settings_from_config_accepts_documented_settings():
    parser = build_parser()
    args = parser.parse_args(
        [
            "make-map",
            "--cubes",
            "cube.fits",
            "--feature",
            "C60_18p9",
            "--outdir",
            "results",
        ]
    )

    settings = _settings_from_args(
        args,
        {
            "settings": {
                "continuum_mode": "morph_spline",
                "morph_half_window": 12,
                "min_finite": 8,
                "min_anchor_points": 3,
                "fallback_anchor_count": 5,
                "include_edge_anchors": False,
                "min_feature_coverage": 0.9,
                "min_required_coverage": 0.85,
                "min_relative_weight": 0.2,
                "edge_erosion_pixels": 2,
                "output_unit": "native",
                "stitch_spectral_gap_tolerance_um": 0.01,
                "stitch_overlap_strategy": "keep",
                "stitch_min_spatial_coverage": 0.9,
                "stitch_require_common_spatial_footprint": False,
            }
        },
    )

    assert settings.output_unit == "native"
    assert settings.continuum.mode == "morph_spline"
    assert settings.continuum.morph_half_window == 12
    assert settings.continuum.min_finite == 8
    assert settings.continuum.min_anchor_points == 3
    assert settings.continuum.fallback_anchor_count == 5
    assert settings.continuum.include_edge_anchors is False
    assert settings.min_feature_coverage == 0.9
    assert settings.min_required_coverage == 0.85
    assert settings.min_relative_weight == 0.2
    assert settings.edge_erosion_pixels == 2
    assert settings.stitching.spectral_gap_tolerance_um == 0.01
    assert settings.stitching.overlap_strategy == "keep"
    assert settings.stitching.min_spatial_coverage == 0.9
    assert settings.stitching.require_common_spatial_footprint is False


def test_cli_args_override_config_settings():
    parser = build_parser()
    args = parser.parse_args(
        [
            "make-map",
            "--cubes",
            "cube.fits",
            "--feature",
            "C60_18p9",
            "--outdir",
            "results",
            "--continuum-mode",
            "spline",
            "--morph-half-window",
            "13",
        ]
    )

    settings = _settings_from_args(
        args,
        {
            "settings": {
                "continuum_mode": "morph_spline",
                "morph_half_window": 3,
            }
        },
    )

    assert settings.continuum.mode == "spline"
    assert settings.continuum.morph_half_window == 13


def test_unknown_config_setting_fails_cleanly():
    parser = build_parser()
    args = parser.parse_args(
        [
            "make-map",
            "--cubes",
            "cube.fits",
            "--feature",
            "C60_18p9",
            "--outdir",
            "results",
        ]
    )

    with pytest.raises(ValueError, match="Unknown setting"):
        _settings_from_args(args, {"settings": {"not_a_real_setting": True}})


def test_main_unknown_feature_reports_clean_error(capsys):
    status = main(
        [
            "make-map",
            "--cubes",
            "missing.fits",
            "--feature",
            "NOPE",
            "--outdir",
            "results",
        ]
    )

    captured = capsys.readouterr()
    assert status == 1
    assert "ERROR: Unknown feature 'NOPE'." in captured.err
    assert "Traceback" not in captured.err


def test_main_bad_custom_anchors_reports_clean_error(capsys):
    status = main(
        [
            "make-map",
            "--cubes",
            "missing.fits",
            "--wavelength",
            "18.9",
            "--feature-window",
            "18.0",
            "19.0",
            "--anchors",
            "18.1",
            "--outdir",
            "results",
        ]
    )

    captured = capsys.readouterr()
    assert status == 1
    assert "ERROR: Anchor values must be supplied as wavelength pairs." in captured.err
    assert "Traceback" not in captured.err


def test_main_missing_file_reports_clean_error(capsys):
    status = main(
        [
            "make-map",
            "--cubes",
            "does_not_exist.fits",
            "--feature",
            "C60_18p9",
            "--outdir",
            "results",
        ]
    )

    captured = capsys.readouterr()
    assert status == 1
    assert "ERROR:" in captured.err
    assert "does_not_exist.fits" in captured.err
    assert "Traceback" not in captured.err
