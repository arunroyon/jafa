import subprocess
import sys

import pytest

from jafa.exceptions import OptionalDependencyError
from jafa.utils import require_dependency


def test_top_level_import_does_not_load_heavy_optional_modules():
    code = (
        "import sys; "
        "sys.path.insert(0, 'src'); "
        "import jafa; "
        "forbidden={'spectral_cube','pybaselines','reproject','matplotlib'}; "
        "loaded=sorted(forbidden.intersection(sys.modules)); "
        "assert not loaded, loaded; "
        "print(jafa.__version__)"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=".", text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr


def test_legacy_import_alias_still_loads_feature_database():
    code = (
        "import sys; "
        "sys.path.insert(0, 'src'); "
        "import jwst_feature_mapper; "
        "from jwst_feature_mapper.feature_db import load_feature_database; "
        "features=load_feature_database(); "
        "assert 'C60_18p9' in features; "
        "print(jwst_feature_mapper.__version__)"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=".", text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr


def test_optional_dependency_error_is_informative(monkeypatch):
    def missing_import(name):
        raise ImportError(name)

    monkeypatch.setattr("jafa.utils.import_module", missing_import)
    with pytest.raises(OptionalDependencyError, match=r"pip install jafa\[plot\]"):
        require_dependency("not_real", extra="plot", purpose="testing optional errors")
