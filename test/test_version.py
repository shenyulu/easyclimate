"""
pytest for plot.taylor_diagrams.py
"""

import pytest
import easyclimate as ecl


def test_version():
    result_data = ecl.show_versions()
    refer_data = ecl.__version__
    assert result_data == refer_data
