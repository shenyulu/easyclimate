import pytest
import xarray as xr
import numpy as np
import pandas as pd
from easyclimate.core.monthstat import (
    calc_monthly_mean,
    calc_monthly_sum,
    calc_monthly_std,
    calc_monthly_var,
    calc_monthly_max,
    calc_monthly_min,
)


class TestCalcMonthlyMean:
    """Test class for calc_monthly_mean function"""

    @pytest.fixture
    def sample_dataarray(self):
        """Create sample DataArray for testing"""
        time_index = pd.date_range("2020-01-01", "2020-03-31", freq="D")
        data = np.random.rand(len(time_index), 3, 3)
        return xr.DataArray(data, dims=["time", "x", "y"], coords={"time": time_index})

    def test_basic_functionality(self, sample_dataarray):
        """Test basic monthly mean calculation"""
        result = calc_monthly_mean(sample_dataarray)

        # Check result type
        assert isinstance(result, xr.DataArray)

        # Check dimensions
        assert "time" in result.dims

        # Check time frequency is monthly
        assert len(result.time) == 3  # Jan, Feb, Mar

        # Check values are reasonable (mean should be between min and max of original)
        original_mean = sample_dataarray.mean()
        result_mean = result.mean()
        assert abs(original_mean - result_mean) < 1.0

    def test_custom_dimension(self, sample_dataarray):
        """Test with custom dimension name"""
        # Rename time dimension
        data_renamed = sample_dataarray.rename({"time": "custom_time"})
        result = calc_monthly_mean(data_renamed, dim="custom_time")

        assert "custom_time" in result.dims
        assert len(result.custom_time) == 3

    def test_kwargs_passed(self, sample_dataarray):
        """Test that kwargs are properly passed to underlying function"""
        # This should not raise an error
        result = calc_monthly_mean(sample_dataarray, skipna=True)
        assert isinstance(result, xr.DataArray)


class TestCalcMonthlySum:
    """Test class for calc_monthly_sum function"""

    @pytest.fixture
    def sample_dataarray_ones(self):
        """Create sample DataArray with ones for sum testing"""
        time_index = pd.date_range("2020-01-01", "2020-01-31", freq="D")
        data = np.ones(len(time_index))
        return xr.DataArray(data, dims=["time"], coords={"time": time_index})

    def test_monthly_sum_calculation(self, sample_dataarray_ones):
        """Test monthly sum calculation with known values"""
        result = calc_monthly_sum(sample_dataarray_ones)

        assert isinstance(result, xr.DataArray)
        assert len(result.time) == 1  # Only January

        # January has 31 days, sum should be 31
        assert result.values[0] == 31

    def test_multiple_months_sum(self):
        """Test sum over multiple months"""
        time_index = pd.date_range("2020-01-01", "2020-02-15", freq="D")
        data = np.ones(len(time_index))
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_sum(da)

        assert len(result.time) == 2  # January and February
        assert result.sel(time="2020-01").values == 31  # Jan has 31 days
        assert result.sel(time="2020-02").values == 15  # Feb has 15 days in test data


class TestCalcMonthlyStd:
    """Test class for calc_monthly_std function"""

    @pytest.fixture
    def sample_dataarray_constant(self):
        """Create sample DataArray with constant values for std testing"""
        time_index = pd.date_range("2020-01-01", "2020-01-31", freq="D")
        data = np.full(len(time_index), 5.0)  # All values are 5.0
        return xr.DataArray(data, dims=["time"], coords={"time": time_index})

    def test_std_of_constant(self, sample_dataarray_constant):
        """Test std of constant values should be zero"""
        result = calc_monthly_std(sample_dataarray_constant)

        assert isinstance(result, xr.DataArray)
        # Standard deviation of constant values should be 0
        assert result.values[0] == pytest.approx(0.0, abs=1e-10)

    def test_std_calculation(self):
        """Test std calculation with known variance"""
        time_index = pd.date_range("2020-01-01", "2020-01-10", freq="D")
        data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_std(da)

        expected_std = np.std([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        assert result.values[0] == pytest.approx(expected_std, rel=1e-10)


class TestCalcMonthlyVar:
    """Test class for calc_monthly_var function"""

    def test_variance_calculation(self):
        """Test variance calculation with known values"""
        time_index = pd.date_range("2020-01-01", "2020-01-05", freq="D")
        data = np.array([1, 2, 3, 4, 5])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_var(da)

        expected_var = np.var([1, 2, 3, 4, 5])
        assert result.values[0] == pytest.approx(expected_var, rel=1e-10)

    def test_relationship_with_std(self):
        """Test that variance is square of standard deviation"""
        time_index = pd.date_range("2020-01-01", "2020-01-10", freq="D")
        data = np.random.rand(10)
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        std_result = calc_monthly_std(da)
        var_result = calc_monthly_var(da)

        # Variance should be square of standard deviation
        assert var_result.values[0] == pytest.approx(
            std_result.values[0] ** 2, rel=1e-10
        )


class TestCalcMonthlyMax:
    """Test class for calc_monthly_max function"""

    def test_max_calculation(self):
        """Test maximum calculation with known values"""
        time_index = pd.date_range("2020-01-01", "2020-01-05", freq="D")
        data = np.array([1, 5, 3, 2, 4])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_max(da)

        assert result.values[0] == 5  # Maximum value

    def test_multiple_months_max(self):
        """Test max over multiple months"""
        time_index = pd.date_range("2020-01-25", "2020-02-06", freq="D")
        jan_data = np.array([1, 2, 3, 4, 5, 6, 7])  # Jan 25-31
        feb_data = np.array([10, 9, 8, 7, 6, 5])  # Feb 1-5
        data = np.concatenate([jan_data, feb_data])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_max(da)

        jan_max = result.sel(time="2020-01").values
        feb_max = result.sel(time="2020-02").values

        assert jan_max == 7
        assert feb_max == 10


class TestCalcMonthlyMin:
    """Test class for calc_monthly_min function"""

    def test_min_calculation(self):
        """Test minimum calculation with known values"""
        time_index = pd.date_range("2020-01-01", "2020-01-05", freq="D")
        data = np.array([5, 2, 3, 1, 4])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_min(da)

        assert result.values[0] == 1  # Minimum value

    def test_edge_cases(self):
        """Test edge cases for minimum calculation"""
        # Single day month
        time_index = pd.date_range("2020-01-01", "2020-01-01", freq="D")
        data = np.array([42])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        result = calc_monthly_min(da)

        assert result.values[0] == 42


class TestIntegration:
    """Integration tests for all functions"""

    @pytest.fixture
    def comprehensive_sample_data(self):
        """Create comprehensive sample data for integration testing"""
        time_index = pd.date_range("2020-01-01", "2020-12-31", freq="D")
        # Create data with some pattern
        data = np.sin(np.arange(len(time_index)) * 2 * np.pi / 365) + 10
        return xr.DataArray(
            data,
            dims=["time"],
            coords={"time": time_index},
            attrs={"description": "Test data"},
        )

    def test_all_functions_with_comprehensive_data(self, comprehensive_sample_data):
        """Test all functions with a full year of data"""
        # Test that all functions work without errors
        functions = [
            calc_monthly_mean,
            calc_monthly_sum,
            calc_monthly_std,
            calc_monthly_var,
            calc_monthly_max,
            calc_monthly_min,
        ]

        for func in functions:
            result = func(comprehensive_sample_data)
            assert isinstance(result, xr.DataArray)
            assert len(result.time) == 12  # 12 months
            # Check that attributes are preserved (if applicable)
            assert hasattr(result, "attrs")

    def test_result_consistency(self):
        """Test consistency between different statistical functions"""
        time_index = pd.date_range("2020-01-01", "2020-01-10", freq="D")
        data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        da = xr.DataArray(data, dims=["time"], coords={"time": time_index})

        min_val = calc_monthly_min(da).values[0]
        max_val = calc_monthly_max(da).values[0]
        mean_val = calc_monthly_mean(da).values[0]

        # Basic consistency checks
        assert min_val <= mean_val <= max_val
        assert min_val == 1
        assert max_val == 10
        assert mean_val == pytest.approx(5.5)
