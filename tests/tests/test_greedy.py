import os
import glob
import shutil
import pandas as pd
import importlib
import logging
import subprocess
import pytest

import tests.test_utils as test_utils

@pytest.fixture()
def self():
    self = test_utils.load_common_self()
    yield self
    self.teardown()

def test_multiround_primary_clustering_1(self):
    assert True