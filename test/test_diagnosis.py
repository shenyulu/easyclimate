import pytest

import easyclimate as ecl
import numpy as np
from .util import round_sf_np

def test_get_coriolis_parameter():
    latdata = np.array([30, 60])
    x = ecl.get_coriolis_parameter(latdata, omega = 7.292e-5)
    x = round_sf_np(x)
    assert (x == np.array([7.292e-05, 1.263e-04])).all()
