"""
pytest for field/equatorial_wave/wk_spectra
"""

import pytest
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import easyclimate as ecl
from pathlib import Path
from .const_define import DOCS_DATA_PATH

data = (
    xr.open_dataset(str(Path(DOCS_DATA_PATH, "olr_smooth_data.nc")))
    .sortby("lat")
    .olr.sel(lat=slice(-15, 15))
)

# var
spd = 1
nDayWin = 96
nDaySkip = -71

data_dt = ecl.field.equatorial_wave.remove_dominant_signals(
    data, spd, nDayWin, nDaySkip
)
data_as = ecl.field.equatorial_wave.decompose_symasym(data_dt)
psum = ecl.field.equatorial_wave.calc_spectral_coefficients(
    data_as, spd, nDayWin, nDaySkip
)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_remove_dominant_signals():
    fig, ax = plt.subplots()
    data_dt.isel(time=0).plot.contourf(levels=21)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_decompose_symasym():
    fig, ax = plt.subplots()
    data_as.isel(time=0).plot.contourf(levels=21)
    return fig


def test_decompose_symasym_value_error():
    """Test that ValueError is raised when input DataArray lacks the specified lat_dim."""

    # Create test data without latitude dimension
    test_data = np.random.rand(10, 5)  # 10x5 array with no explicit latitude
    da_no_lat = xr.DataArray(test_data, dims=["time", "lon"])

    # Test with default lat_dim="lat"
    with pytest.raises(ValueError) as excinfo:
        ecl.field.equatorial_wave.decompose_symasym(da_no_lat)
    assert "Input DataArray must have lat dimension" in str(excinfo.value)

    # Test with custom lat_dim name
    da_wrong_dim = xr.DataArray(test_data, dims=["time", "latitude"])
    with pytest.raises(ValueError) as excinfo:
        ecl.field.equatorial_wave.decompose_symasym(da_wrong_dim, lat_dim="lat")
    assert "Input DataArray must have lat dimension" in str(excinfo.value)

    # Test case where dimension exists but with different name
    da_diff_lat_name = xr.DataArray(test_data, dims=["latitude", "lon"])
    with pytest.raises(ValueError) as excinfo:
        ecl.field.equatorial_wave.decompose_symasym(da_diff_lat_name)
    assert "Input DataArray must have lat dimension" in str(excinfo.value)

    # Test case where we specify the correct dimension name
    # This should NOT raise an error
    da_with_lat = xr.DataArray(test_data, dims=["lat", "lon"])
    try:
        ecl.field.equatorial_wave.decompose_symasym(da_with_lat)
    except ValueError:
        pytest.fail("Unexpected ValueError for DataArray with correct lat dimension")

    # Test case with alternative dimension name when specified
    da_with_latitude = xr.DataArray(test_data, dims=["latitude", "lon"])
    try:
        ecl.field.equatorial_wave.decompose_symasym(
            da_with_latitude, lat_dim="latitude"
        )
    except ValueError:
        pytest.fail(
            "Unexpected ValueError when correct alternative lat_dim is specified"
        )


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_spectral_coefficients_1():
    fig, ax = plt.subplots()
    psum.psumanti_r.plot.contourf(ax=ax, levels=21, cmap="YlGnBu")
    ecl.field.equatorial_wave.draw_wk_anti_analysis()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=20)
def test_calc_spectral_coefficients_2():
    fig, ax = plt.subplots()
    psum.psumsym_r.plot.contourf(ax=ax, levels=21, cmap="YlGnBu")
    ecl.field.equatorial_wave.draw_wk_sym_analysis()
    return fig
