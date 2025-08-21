"""
pytest for filter/spectrum.py
"""

import numpy as np
import xarray as xr
import pytest

from easyclimate.filter.spectrum import (
    calc_time_spectrum,
    calc_mean_fourier_amplitude,
    filter_fourier_harmonic_analysis,
)


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


class TestFilterFourierHarmonicAnalysis:
    @pytest.fixture
    def sample_data(self):
        """Create a sample DataArray for testing."""
        np.random.seed(42)
        return xr.DataArray(
            np.random.randn(56, 90, 180),
            dims=["time", "lat", "lon"],
            coords={
                "time": np.arange(1948, 2004),
                "lat": np.linspace(-90, 90, 90),
                "lon": np.linspace(0, 360, 180, endpoint=False),
            },
            name="z200",
        )

    def test_invalid_time_dim(self, sample_data):
        """Test error for invalid time dimension."""
        with pytest.raises(ValueError, match="Time dimension 'invalid_dim' not found"):
            filter_fourier_harmonic_analysis(sample_data, time_dim="invalid_dim")

    def test_invalid_filter_type(self, sample_data):
        """Test error for invalid filter type."""
        with pytest.raises(ValueError, match="Invalid filter_type: invalid"):
            filter_fourier_harmonic_analysis(sample_data, filter_type="invalid")

    def test_highpass_missing_max_period(self, sample_data):
        """Test error for highpass filter with missing max_period."""
        with pytest.raises(ValueError, match="High-pass filter requires max_period"):
            filter_fourier_harmonic_analysis(
                sample_data, filter_type="highpass", period_bounds=(None, None)
            )

    def test_lowpass_missing_min_period(self, sample_data):
        """Test error for lowpass filter with missing min_period."""
        with pytest.raises(ValueError, match="Low-pass filter requires min_period"):
            filter_fourier_harmonic_analysis(
                sample_data, filter_type="lowpass", period_bounds=(None, None)
            )

    def test_bandpass_missing_periods(self, sample_data):
        """Test error for bandpass filter with missing period bounds."""
        with pytest.raises(
            ValueError, match="Bandpass filter requires both min_period and max_period"
        ):
            filter_fourier_harmonic_analysis(
                sample_data, filter_type="bandpass", period_bounds=(None, None)
            )

    def test_invalid_period_bounds(self, sample_data):
        """Test error for min_period >= max_period."""
        with pytest.raises(
            ValueError, match="min_period.*must be less than max_period"
        ):
            filter_fourier_harmonic_analysis(sample_data, period_bounds=(10.0, 5.0))

    def test_output_dimensions(self, sample_data):
        """Test that output dimensions match input dimensions."""
        result = filter_fourier_harmonic_analysis(
            sample_data, period_bounds=(2.0, 10.0), filter_type="bandpass"
        )
        assert result.dims == sample_data.dims
        assert result.shape == sample_data.shape
        assert result.coords.equals(sample_data.coords)

    def test_metadata_attributes(self, sample_data):
        """Test that filter metadata is added to output attributes."""
        result = filter_fourier_harmonic_analysis(
            sample_data,
            period_bounds=(2.0, 10.0),
            filter_type="bandpass",
            sampling_interval=1.0,
            apply_window=True,
        )
        assert result.attrs["filter_type"] == "bandpass"
        assert result.attrs["period_bounds"] == (2.0, 10.0)
        assert result.attrs["sampling_interval"] == 1.0
        assert result.attrs["window_applied"] is True

    def test_window_application(self, sample_data):
        """Test that Hann window is applied when requested."""
        result_with_window = filter_fourier_harmonic_analysis(
            sample_data,
            period_bounds=(2.0, 10.0),
            filter_type="bandpass",
            apply_window=True,
        )
        result_without_window = filter_fourier_harmonic_analysis(
            sample_data,
            period_bounds=(2.0, 10.0),
            filter_type="bandpass",
            apply_window=False,
        )
        assert not np.allclose(result_with_window.values, result_without_window.values)

    def test_filter_types(self, sample_data):
        """Test different filter types produce different results."""
        highpass = filter_fourier_harmonic_analysis(
            sample_data, period_bounds=(None, 5.0), filter_type="highpass"
        )
        lowpass = filter_fourier_harmonic_analysis(
            sample_data, period_bounds=(5.0, None), filter_type="lowpass"
        )
        bandpass = filter_fourier_harmonic_analysis(
            sample_data, period_bounds=(2.0, 10.0), filter_type="bandpass"
        )
        assert not np.allclose(highpass.values, lowpass.values)
        assert not np.allclose(highpass.values, bandpass.values)
        assert not np.allclose(lowpass.values, bandpass.values)
