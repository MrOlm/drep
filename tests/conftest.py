import os
import pytest
import tests.test_utils as test_utils

def pytest_configure(config):
    test_utils.fix_bdb_paths(test_utils.load_solutions_wd())
    config.addinivalue_line("markers", "requires_nsimscan: requires nsimscan (goANI/gANI), not available in CI")
    config.addinivalue_line("markers", "requires_checkm: requires checkm database, not available in CI")
    config.addinivalue_line("markers", "requires_mummer: mummer version-sensitive test, skipped in CI")

def pytest_collection_modifyitems(config, items):
    if os.environ.get("CI"):
        skip_nsimscan = pytest.mark.skip(reason="nsimscan not available in CI (https://github.com/abadona/qsimscan)")
        skip_checkm = pytest.mark.skip(reason="checkm database not set up in CI")
        skip_mummer = pytest.mark.skip(reason="mummer version-sensitive, skipped in CI")
        for item in items:
            if "requires_nsimscan" in item.keywords:
                item.add_marker(skip_nsimscan)
            if "requires_checkm" in item.keywords:
                item.add_marker(skip_checkm)
            if "requires_mummer" in item.keywords:
                item.add_marker(skip_mummer)
