"""
pytest for physics/geo/linrood_latwgt.py
"""

import pytest
import numpy as np
from easyclimate.physics.geo.linrood_latwgt import calc_lat_weight_lin_rood


class TestCalcLatWeightLinRood:

    def test_calc_lat_weight_lin_rood1(self):
        lat, wgt = calc_lat_weight_lin_rood(8)
        refer_data1 = np.array(
            [
                -90.0,
                -64.28571429,
                -38.57142857,
                -12.85714286,
                12.85714286,
                38.57142857,
                64.28571429,
                90.0,
            ]
        )
        refer_data2 = np.array(
            [
                0.02507209,
                0.19309643,
                0.34794774,
                0.43388374,
                0.43388374,
                0.34794774,
                0.19309643,
                0.02507209,
            ]
        )
        assert np.isclose(lat, refer_data1).all()
        assert np.isclose(wgt, refer_data2).all()

    def test_basic_output_structure(self):
        """Test that the function returns a tuple of two numpy arrays."""
        nlat = 5
        lat, weight = calc_lat_weight_lin_rood(nlat)

        assert isinstance(lat, np.ndarray)
        assert isinstance(weight, np.ndarray)
        assert lat.shape == (nlat,)
        assert weight.shape == (nlat,)

    def test_pole_latitudes(self):
        """Test that the first and last latitudes are -90 and 90 degrees."""
        nlat = 10
        lat, _ = calc_lat_weight_lin_rood(nlat)

        assert lat[0] == -90.0
        assert lat[-1] == 90.0

    def test_uniform_latitude_spacing(self):
        """Test that latitudes are uniformly spaced."""
        nlat = 5
        lat, _ = calc_lat_weight_lin_rood(nlat)

        expected_lats = np.linspace(-90, 90, nlat)
        np.testing.assert_allclose(lat, expected_lats, rtol=1e-10)

    def test_weight_sum(self):
        """Test that the sum of weights equals 2 (integral over sphere)."""
        nlat = 20
        _, weight = calc_lat_weight_lin_rood(nlat)

        assert np.isclose(np.sum(weight), 2.0, rtol=1e-10)

    def test_minimum_nlat(self):
        """Test that nlat < 2 raises ValueError."""
        with pytest.raises(ValueError):
            calc_lat_weight_lin_rood(1)

    def test_against_known_values(self):
        """Test against precomputed values for nlat=3."""
        nlat = 3
        lat, weight = calc_lat_weight_lin_rood(nlat)

        # Expected latitudes: -90, 0, 90
        expected_lat = np.array([-90.0, 0.0, 90.0])
        np.testing.assert_allclose(lat, expected_lat)

        # Expected weights: [1 + sin(-45°), sin(45°) - sin(-45°), 1 - sin(45°)]
        sqrt2 = np.sqrt(2)
        expected_weight = np.array(
            [1.0 + (-sqrt2 / 2), sqrt2 / 2 - (-sqrt2 / 2), 1.0 - sqrt2 / 2]
        )
        np.testing.assert_allclose(weight, expected_weight, rtol=1e-10)

    def test_weight_properties(self):
        """Test that weights are positive and symmetric for even nlat."""
        nlat = 10
        _, weight = calc_lat_weight_lin_rood(nlat)

        assert np.all(weight > 0)
        np.testing.assert_allclose(weight, weight[::-1], rtol=1e-10)  # Symmetry check
