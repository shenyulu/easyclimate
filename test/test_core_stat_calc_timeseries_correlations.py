import pytest
import xarray as xr
import numpy as np
import pandas as pd
import easyclimate as ecl


class TestCalcTimeseriesCorrelations:
    """Test suite for the calc_timeseries_correlations function."""

    @pytest.fixture
    def sample_data(self):
        """Create sample DataArrays for testing."""
        time = pd.date_range("2020-01-01", "2020-01-10", freq="D")
        data1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=float)
        data2 = data1 * 0.5 + np.random.randn(len(time)) * 0.1
        data3 = np.random.randn(len(time))

        da1 = xr.DataArray(data1, dims="time", coords={"time": time}, name="series1")
        da2 = xr.DataArray(data2, dims="time", coords={"time": time}, name="series2")
        da3 = xr.DataArray(data3, dims="time", coords={"time": time}, name="series3")

        return {
            "dict_input": {"series1": da1, "series2": da2, "series3": da3},
            "list_input": [da1, da2, da3],
        }

    def test_dict_input(self, sample_data):
        """Test function with dictionary input."""
        result = ecl.calc_timeseries_correlations(sample_data["dict_input"], dim="time")

        # Check output type and structure
        assert isinstance(result, xr.DataArray)
        assert result.dims == ("var1", "var2")
        assert list(result.coords["var1"].values) == ["series1", "series2", "series3"]
        assert list(result.coords["var2"].values) == ["series1", "series2", "series3"]
        assert result.name == "correlation"

        # Check matrix properties
        assert result.shape == (3, 3)
        assert np.allclose(result.values.diagonal(), 1.0)  # Diagonal should be 1
        assert np.allclose(result.values, result.values.T)  # Matrix should be symmetric

    def test_list_input(self, sample_data):
        """Test function with list input."""
        result = ecl.calc_timeseries_correlations(sample_data["list_input"], dim="time")

        # Check output type and structure
        assert isinstance(result, xr.DataArray)
        assert result.dims == ("var1", "var2")
        assert list(result.coords["var1"].values) == ["var_0", "var_1", "var_2"]
        assert list(result.coords["var2"].values) == ["var_0", "var_1", "var_2"]
        assert result.name == "correlation"

        # Check matrix properties
        assert result.shape == (3, 3)
        assert np.allclose(result.values.diagonal(), 1.0)
        assert np.allclose(result.values, result.values.T)

    def test_empty_input(self):
        """Test that empty input raises ValueError."""
        with pytest.raises(ValueError, match="data_arrays cannot be empty"):
            ecl.calc_timeseries_correlations({})
        with pytest.raises(ValueError, match="data_arrays cannot be empty"):
            ecl.calc_timeseries_correlations([])

    def test_invalid_input_type(self, sample_data):
        """Test that non-DataArray inputs raise TypeError."""
        invalid_input = sample_data["dict_input"].copy()
        invalid_input["series1"] = [1, 2, 3]  # Not a DataArray
        with pytest.raises(TypeError, match="All inputs must be xarray.DataArray"):
            ecl.calc_timeseries_correlations(invalid_input)

        invalid_list = sample_data["list_input"].copy()
        invalid_list[0] = [1, 2, 3]
        with pytest.raises(TypeError, match="All inputs must be xarray.DataArray"):
            ecl.calc_timeseries_correlations(invalid_list)

    def test_missing_dimension(self, sample_data):
        """Test that missing dimension raises ValueError."""
        # Create DataArray without 'time' dimension
        da = xr.DataArray([1, 2, 3], dims="x", coords={"x": [0, 1, 2]})
        invalid_input = sample_data["dict_input"].copy()
        invalid_input["series1"] = da
        with pytest.raises(
            ValueError, match="All DataArrays must contain the 'time' dimension"
        ):
            ecl.calc_timeseries_correlations(invalid_input)

    def test_correlation_values(self, sample_data):
        """Test correctness of correlation values."""
        result = ecl.calc_timeseries_correlations(sample_data["dict_input"], dim="time")

        # Verify specific correlation values
        # series1 and series2 should have high positive correlation
        corr_12 = result.sel(var1="series1", var2="series2").values
        assert (
            0.5 < corr_12 < 1.0
        )  # Expect positive correlation due to data construction

        # series1 and series3 should have near-zero correlation
        corr_13 = result.sel(var1="series1", var2="series3").values
        assert -0.5 < corr_13 < 0.5  # Expect near-zero correlation with random data

    def test_nan_handling(self):
        """Test handling of NaN values in input DataArrays."""
        time = pd.date_range("2020-01-01", "2020-01-10", freq="D")
        data1 = np.array([1, 2, np.nan, 4, 5, 6, 7, 8, 9, 10], dtype=float)
        data2 = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20], dtype=float)

        da1 = xr.DataArray(data1, dims="time", coords={"time": time}, name="series1")
        da2 = xr.DataArray(data2, dims="time", coords={"time": time}, name="series2")

        result = ecl.calc_timeseries_correlations(
            {"series1": da1, "series2": da2}, dim="time"
        )

        # Check that correlation is computed correctly despite NaNs
        assert not np.isnan(result.sel(var1="series1", var2="series2").values)
        assert np.allclose(result.values.diagonal(), 1.0)
