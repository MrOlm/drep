import os
import pytest
import tests.test_utils as test_utils

def pytest_configure(config):
    test_utils.fix_bdb_paths(test_utils.load_solutions_wd())
    config.addinivalue_line("markers", "requires_nsimscan: requires nsimscan (goANI/gANI), not available in CI")

def pytest_collection_modifyitems(config, items):
    if os.environ.get("CI"):
        skip = pytest.mark.skip(reason="nsimscan not available in CI (https://github.com/abadona/qsimscan)")
        for item in items:
            if "requires_nsimscan" in item.keywords:
                item.add_marker(skip)
