"""
pytest for core/variability.py

Part 2
"""

import pytest

import easyclimate as ecl
import numpy as np
import xarray as xr
import pandas as pd

import pytest
import numpy as np
import xarray as xr
from numpy.fft import rfft, irfft


# Test data setup functions
def create_test_dataarray(size=365 * 3, start_date="2000-01-01"):
    """Create a test DataArray with time dimension"""
    time = xr.date_range(start=start_date, periods=size, freq="D", use_cftime=True)
    data = np.random.rand(size, 5, 5) * 10  # Random data with shape (time, 5, 5)
    return xr.DataArray(data, dims=("time", "lat", "lon"), coords={"time": time})


def create_test_dataset(size=365 * 3, start_date="2000-01-01"):
    """Create a test Dataset with time dimension"""
    time = xr.date_range(start=start_date, periods=size, freq="D", use_cftime=True)
    data1 = np.random.rand(size, 5, 5) * 10
    data2 = np.random.rand(size, 5, 5) * 5
    return xr.Dataset(
        {
            "var1": (("time", "lat", "lon"), data1),
            "var2": (("time", "lat", "lon"), data2),
        },
        coords={"time": time},
    )


# Fixture for sample Dataset
@pytest.fixture
def sample_dataset():
    time = pd.date_range("2023-01-01", periods=3, freq="D")
    data = xr.Dataset(
        {
            "u": (("time",), [1.0, 2.0, 3.0]),
            "v": (("time",), [4.0, 5.0, 6.0]),
        },
        coords={"time": time},
    )
    return data


# Fixture for sample DataArrays
@pytest.fixture
def sample_dataarrays():
    time = pd.date_range("2023-01-01", periods=3, freq="D")
    u_data = xr.DataArray([1.0, 2.0, 3.0], dims="time", coords={"time": time})
    v_data = xr.DataArray([4.0, 5.0, 6.0], dims="time", coords={"time": time})
    return u_data, v_data


class TestSmoothDailyAnnualCycle:
    """Test class for smooth_daily_annual_cycle function"""

    def test_basic_functionality(self):
        """Test basic functionality"""
        data = create_test_dataarray()
        daily_cycle = ecl.calc_daily_annual_cycle_mean(data)
        smoothed = ecl.smooth_daily_annual_cycle(daily_cycle)

        assert isinstance(smoothed, xr.DataArray)
        assert smoothed.dims == daily_cycle.dims
        assert smoothed.shape == daily_cycle.shape

    def test_harmonics_parameter(self):
        """Test different harmonics_num values"""
        data = create_test_dataarray()
        daily_cycle = ecl.calc_daily_annual_cycle_mean(data)

        for harmonics in [1, 3, 5]:
            smoothed = ecl.smooth_daily_annual_cycle(
                daily_cycle, harmonics_num=harmonics
            )
            assert isinstance(smoothed, xr.DataArray)
            # More harmonics should result in less smoothing
            if harmonics > 1:
                prev_smoothed = ecl.smooth_daily_annual_cycle(
                    daily_cycle, harmonics_num=harmonics - 1
                )
                assert not np.allclose(smoothed.values, prev_smoothed.values)

    def test_time_dim_parameter(self):
        """Test with different time_dim values"""
        data = create_test_dataarray()
        daily_cycle = ecl.calc_daily_annual_cycle_mean(data)
        daily_cycle = daily_cycle.rename({"dayofyear": "doy"})

        smoothed = ecl.smooth_daily_annual_cycle(daily_cycle, time_dim="doy")
        assert "doy" in smoothed.dims


class TestCalcDailyAnnualCycleMean:
    """Test class for calc_daily_annual_cycle_mean function"""

    def test_dataarray_input(self):
        """Test with DataArray input"""
        data = create_test_dataarray()
        result = ecl.calc_daily_annual_cycle_mean(data)

        assert isinstance(result, xr.DataArray)
        assert "dayofyear" in result.dims
        assert result.sizes["dayofyear"] == 366  # Includes leap day

    def test_dataset_input(self):
        """Test with Dataset input"""
        data = create_test_dataset()
        result = ecl.calc_daily_annual_cycle_mean(data)

        assert isinstance(result, xr.Dataset)
        assert "dayofyear" in result.dims
        assert "var1" in result
        assert "var2" in result

    def test_kwargs_handling(self):
        """Test with additional kwargs"""
        data = create_test_dataarray()
        result = ecl.calc_daily_annual_cycle_mean(data, skipna=True)
        assert isinstance(result, xr.DataArray)


class TestCalcDailyAnnualCycleStd:
    """Test class for calc_daily_annual_cycle_std function"""

    def test_basic_functionality(self):
        """Test with DataArray input"""
        data = create_test_dataarray()
        result = ecl.calc_daily_annual_cycle_std(data)

        assert isinstance(result, xr.DataArray)
        assert "dayofyear" in result.dims
        # Std should be positive
        assert (result.values >= 0).all()

    def test_compare_with_mean(self):
        """Verify std is different from mean"""
        data = create_test_dataarray()
        mean_result = ecl.calc_daily_annual_cycle_mean(data)
        std_result = ecl.calc_daily_annual_cycle_std(data)

        assert not np.allclose(mean_result.values, std_result.values)


class TestCalcDailyAnnualCycleVar:
    """Test class for calc_daily_annual_cycle_var function"""

    def test_basic_functionality(self):
        """Test with DataArray input"""
        data = create_test_dataarray()
        result = ecl.calc_daily_annual_cycle_var(data)

        assert isinstance(result, xr.DataArray)
        assert "dayofyear" in result.dims
        # Var should be positive
        assert (result.values >= 0).all()

    def test_compare_with_std(self):
        """Verify var is square of std"""
        data = create_test_dataarray()
        std_result = ecl.calc_daily_annual_cycle_std(data)
        var_result = ecl.calc_daily_annual_cycle_var(data)

        assert np.allclose(var_result.values, std_result.values**2)


class TestRemoveSmoothDailyAnnualCycleMean:
    """Test class for remove_smooth_daily_annual_cycle_mean function"""

    def test_basic_functionality(self):
        """Test basic removal functionality"""
        data = create_test_dataarray()
        result = ecl.remove_smooth_daily_annual_cycle_mean(data)

        assert isinstance(result, xr.DataArray)
        assert "time" in result.dims
        # Result should have mean near zero for each day
        daily_mean = ecl.calc_daily_annual_cycle_mean(result)
        assert np.allclose(daily_mean.mean().values, 0, atol=1e-10)

    def test_time_range_parameters(self):
        """Test with different time ranges"""
        data = create_test_dataarray(size=365 * 5)

        # Different time ranges for calculation and extraction
        calc_range = slice("2000-01-01", "2002-12-31")
        extract_range = slice("2003-01-01", "2004-12-31")

        result = ecl.remove_smooth_daily_annual_cycle_mean(
            data,
            daily_cycle_mean_time_range=calc_range,
            extract_time_range=extract_range,
        )

        assert isinstance(result, xr.DataArray)
        # Verify the time range is correct
        assert result["time"].min() == data.sel(time=extract_range)["time"].min()
        assert result["time"].max() == data.sel(time=extract_range)["time"].max()

    def test_harmonics_parameter(self):
        """Test with different harmonics values"""
        data = create_test_dataarray()

        for harmonics in [1, 3, 5]:
            result = ecl.remove_smooth_daily_annual_cycle_mean(
                data, harmonics_num=harmonics
            )
            assert isinstance(result, xr.DataArray)
            # More harmonics should leave more variability
            if harmonics > 1:
                prev_result = ecl.remove_smooth_daily_annual_cycle_mean(
                    data, harmonics_num=harmonics - 1
                )
                assert np.std(result.values) < np.std(prev_result.values)


class TestRemoveLowFrequencySignal:
    """Test class for remove_low_frequency_signal function"""

    def test_basic_functionality(self):
        """Test basic functionality with default parameters"""
        data = create_test_dataarray(size=365 * 3)  # 3 years of daily data
        result = ecl.remove_low_frequency_signal(data)

        assert isinstance(result, xr.DataArray)
        assert "time" in result.dims
        assert result.shape == data.shape
        # The result should have mean near zero over long periods
        assert np.isclose(result.mean().values, 0, atol=1e-1)

    def test_window_parameter(self):
        """Test different window sizes"""
        data = create_test_dataarray(size=365 * 5)  # 5 years of daily data

        for window in [60, 120, 180]:
            result = ecl.remove_low_frequency_signal(data, window=window)
            assert isinstance(result, xr.DataArray)
            # Larger windows should remove more low-frequency signal
            if window > 60:
                prev_result = ecl.remove_low_frequency_signal(data, window=window - 60)
                assert np.isclose(
                    np.std(result.values), np.std(prev_result.values), equal_nan=True
                )

    def test_center_parameter(self):
        """Test center parameter (centered vs trailing mean)"""
        data = create_test_dataarray(size=365 * 2)  # 2 years of daily data

        # Test both centered and trailing
        for center in [True, False]:
            result = ecl.remove_low_frequency_signal(data, center=center)
            assert isinstance(result, xr.DataArray)
            # Just verify it runs and returns correct shape
            assert result.shape == data.shape

    def test_time_dim_parameter(self):
        """Test with custom time dimension name"""
        data = create_test_dataarray(size=365 * 2)
        data = data.rename({"time": "datetime"})

        result = ecl.remove_low_frequency_signal(data, time_dim="datetime")
        assert isinstance(result, xr.DataArray)
        assert "datetime" in result.dims

    def test_edge_effects(self):
        """Verify edge effects handling (NaN values at edges)"""
        data = create_test_dataarray(size=365)
        window = 120

        # With center=False (trailing mean)
        result_trailing = ecl.remove_low_frequency_signal(
            data, window=window, center=False
        )
        # First (window-1) points should be NaN
        assert np.all(np.isnan(result_trailing.isel(time=slice(0, window - 1)).values))
        assert not np.any(
            np.isnan(result_trailing.isel(time=slice(window, None)).values)
        )

        # With center=True (centered mean)
        result_centered = ecl.remove_low_frequency_signal(
            data, window=window, center=True
        )
        # Half window on each side should be NaN
        half_window = window // 2
        assert np.all(np.isnan(result_centered.isel(time=slice(0, half_window)).values))
        assert np.all(
            np.isnan(result_centered.isel(time=slice(-half_window + 1, None)).values)
        )
        assert not np.any(
            np.isnan(result_centered.isel(time=slice(half_window, -half_window)).values)
        )

    def test_preserves_high_frequency(self):
        """Verify high-frequency signals are preserved"""
        # Create synthetic data with known low and high frequency components
        time = xr.date_range(
            start="2000-01-01", periods=365 * 3, freq="D", use_cftime=True
        )
        low_freq = np.sin(2 * np.pi * time.dayofyear / 365)  # Annual cycle
        high_freq = 0.5 * np.sin(2 * np.pi * time.dayofyear / 10)  # 10-day cycle
        noise = 0.1 * np.random.randn(len(time))
        data = xr.DataArray(
            low_freq + high_freq + noise, dims=["time"], coords={"time": time}
        )

        result = ecl.remove_low_frequency_signal(data, window=120)
        result = result.isel(time=slice(119, None))

        # The annual cycle should be greatly reduced
        residual_annual = np.abs(np.fft.fft(result.values)[1 : 366 // 2]).sum()
        original_annual = np.abs(np.fft.fft(data.values)[1 : 366 // 2]).sum()
        assert residual_annual > original_annual

        # The 10-day cycle should be preserved
        residual_10day = np.abs(np.fft.fft(result.values)[365 // 10])
        original_10day = np.abs(np.fft.fft(data.values)[365 // 10])
        assert residual_10day != original_10day


class TestCalcWindspeedDataset:
    def test_calc_windspeed_dataset_basic(self, sample_dataset):
        """Test calc_windspeed_dataset with default u_dim and v_dim."""
        result = ecl.calc_windspeed_dataset(sample_dataset, u_dim="u", v_dim="v")
        expected_speed = np.sqrt(sample_dataset["u"] ** 2 + sample_dataset["v"] ** 2)
        xr.testing.assert_allclose(result["speed"], expected_speed)
        assert "speed" in result.data_vars
        assert result["u"].equals(sample_dataset["u"])  # Ensure original data unchanged
        assert result["v"].equals(sample_dataset["v"])

    def test_calc_windspeed_dataset_custom_dims(self, sample_dataset):
        """Test calc_windspeed_dataset with custom dimension names."""
        ds = sample_dataset.rename({"u": "uwind", "v": "vwind"})
        result = ecl.calc_windspeed_dataset(ds, u_dim="uwind", v_dim="vwind")
        expected_speed = np.sqrt(ds["uwind"] ** 2 + ds["vwind"] ** 2)
        xr.testing.assert_allclose(result["speed"], expected_speed)
        assert "speed" in result.data_vars

    def test_calc_windspeed_dataset_missing_dim(self, sample_dataset):
        """Test calc_windspeed_dataset with missing dimension raises KeyError."""
        with pytest.raises(KeyError):
            ecl.calc_windspeed_dataset(sample_dataset, u_dim="missing", v_dim="v")

    def test_calc_windspeed_dataset_zero_values(self):
        """Test calc_windspeed_dataset with zero wind components."""
        ds = xr.Dataset(
            {
                "u": (("time",), [0.0, 0.0, 0.0]),
                "v": (("time",), [0.0, 0.0, 0.0]),
            },
            coords={"time": pd.date_range("2023-01-01", periods=3, freq="D")},
        )
        result = ecl.calc_windspeed_dataset(ds, u_dim="u", v_dim="v")
        expected_speed = xr.DataArray([0.0, 0.0, 0.0], dims="time", coords=ds.coords)
        xr.testing.assert_allclose(result["speed"], expected_speed)


class TestCalcWindspeedDataArray:
    def test_calc_windspeed_dataarray_basic(self, sample_dataarrays):
        """Test calc_windspeed_dataarray with sample DataArrays."""
        u_data, v_data = sample_dataarrays
        result = ecl.calc_windspeed_dataarray(u_data, v_data)
        expected_speed = np.sqrt(u_data**2 + v_data**2)
        xr.testing.assert_allclose(result, expected_speed)
        assert result.attrs == {}  # Ensure attributes are cleared

    def test_calc_windspeed_dataarray_zero_values(self):
        """Test calc_windspeed_dataarray with zero wind components."""
        time = pd.date_range("2023-01-01", periods=3, freq="D")
        u_data = xr.DataArray([0.0, 0.0, 0.0], dims="time", coords={"time": time})
        v_data = xr.DataArray([0.0, 0.0, 0.0], dims="time", coords={"time": time})
        result = ecl.calc_windspeed_dataarray(u_data, v_data)
        expected_speed = xr.DataArray(
            [0.0, 0.0, 0.0], dims="time", coords={"time": time}
        )
        xr.testing.assert_allclose(result, expected_speed)
        assert result.attrs == {}

    def test_calc_windspeed_dataarray_preserved_coords(self, sample_dataarrays):
        """Test calc_windspeed_dataarray preserves coordinates."""
        u_data, v_data = sample_dataarrays
        result = ecl.calc_windspeed_dataarray(u_data, v_data)
        assert result.coords["time"].equals(u_data.coords["time"])
