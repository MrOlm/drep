import tests.test_utils as test_utils

def pytest_configure(config):
    test_utils.fix_bdb_paths(test_utils.load_solutions_wd())
