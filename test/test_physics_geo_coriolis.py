"""
pytest for physics.geo.coriolis.py
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr

from easyclimate.physics import get_coriolis_parameter
from .util import round_sf_np


def test_get_coriolis_parameter():
    latdata = np.array([30, 60])
    x = get_coriolis_parameter(latdata, omega=7.292e-5)

    result_data = round_sf_np(x)
    refer_data = np.array([7.292e-05, 1.263e-04])
    assert np.isclose(result_data, refer_data).all()
