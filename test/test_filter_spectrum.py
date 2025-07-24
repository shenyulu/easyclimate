"""
pytest for filter/spectrum.py
"""

import numpy as np
import xarray as xr
import pytest

from easyclimate.filter.spectrum import calc_time_spectrum, calc_mean_fourier_amplitude


class TestSpectrumFunctions:
    @pytest.fixture
    def sample_data(self):
        time = np.arange(100)
        data = np.sin(2 * np.pi * time / 10) + 0.5 * np.random.randn(100)
        da = xr.DataArray(data, dims="time", coords={"time": time}, name="test_signal")
        return da

    def test_calc_time_spectrum_output(self, sample_data):
        result = calc_time_spectrum(sample_data)
        assert "freq" in result
        assert "period" in result
        assert "spectrum_freq" in result
        assert "spectrum_period" in result

        spectrum = result["spectrum_freq"]
        freq = result["freq"]

        # Check if the spectrum values are non-negative
        assert np.all(spectrum.values >= 0)

        # Check if frequency and spectrum dimensions match
        assert spectrum.dims[-1] == "freq"
        assert freq.dims == ("freq",)

    def test_calc_mean_fourier_amplitude_range(self, sample_data):
        mean_amp = calc_mean_fourier_amplitude(sample_data, lower=8, upper=12)

        # Output should not have time dimension
        assert "time" not in mean_amp.dims

        # Check name and attributes
        assert mean_amp.name == "test_signal"
        assert isinstance(mean_amp.values.item(), float)

    def test_calc_mean_fourier_amplitude_zero_range(self, sample_data):
        # Use a very narrow period range that likely excludes everything
        mean_amp = calc_mean_fourier_amplitude(sample_data, lower=1000, upper=2000)
        assert mean_amp.values.item() == pytest.approx(0.0, abs=1e-6)
