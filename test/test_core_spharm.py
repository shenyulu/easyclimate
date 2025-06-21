"""
pytest for spharm.py
"""

import pytest

import numpy as np
import xarray as xr
from easyclimate.core.spharm import (
    calc_gaussian_latitudes,
    calc_geodesic_points,
    calc_spherical_harmonic_coefficients,
    calc_legendre_functions,
    transfer_grid2spectral_transform,
    transfer_spectral_transform2grid,
)


# Test data setup
@pytest.fixture
def sample_grid_data():
    """Create sample grid data for testing transforms"""
    lats = np.linspace(-90, 90, 64)
    lons = np.linspace(0, 360, 128, endpoint=False)
    data = np.random.rand(len(lats), len(lons))
    return xr.DataArray(data, dims=("lat", "lon"), coords={"lat": lats, "lon": lons})


@pytest.fixture
def sample_spectral_data():
    """Create sample spectral data for testing inverse transforms"""
    nspec = 1056  # For T31 truncation (32*33/2)
    data = np.random.rand(nspec) + 1j * np.random.rand(nspec)
    return xr.DataArray(data, dims=("spec_dim",), coords={"spec_dim": np.arange(nspec)})


# Test cases
def test_calc_gaussian_latitudes():
    """Test Gaussian latitudes calculation"""
    nlat = 32
    result = calc_gaussian_latitudes(nlat)

    # Check return type
    assert isinstance(result, xr.Dataset)

    # Check dimensions
    assert "lats" in result
    assert "wts" in result
    assert len(result["lats"]) == nlat
    assert len(result["wts"]) == nlat

    # Check values
    assert np.all(result["lats"] <= 90)
    assert np.all(result["lats"] >= -90)
    assert np.all(result["wts"] > 0)

    # Check metadata
    assert result.attrs["nlat"] == nlat
    assert result["lats"].attrs["name"] == "Gaussian latitudes (degrees)"
    assert result["wts"].attrs["name"] == "Quadrature weights"


def test_calc_geodesic_points():
    """Test geodesic points calculation"""
    m = 10
    result = calc_geodesic_points(m)

    # Check return type
    assert isinstance(result, xr.Dataset)

    # Check dimensions
    assert "lats" in result
    assert "lons" in result
    npoints = 10 * (m - 1) ** 2 + 2
    assert len(result["lats"]) == npoints
    assert len(result["lons"]) == npoints

    # Check values
    assert np.all(result["lats"] <= 90)
    assert np.all(result["lats"] >= -90)
    assert np.all(result["lons"] >= 0)
    assert np.all(result["lons"] <= 360)

    # Check metadata
    assert result["lats"].attrs["m"] == m
    assert result["lons"].attrs["m"] == m


def test_calc_spherical_harmonic_coefficients():
    """Test spherical harmonic coefficients indices calculation"""
    ntrunc = 42
    result = calc_spherical_harmonic_coefficients(ntrunc)

    # Check return type
    assert isinstance(result, xr.Dataset)

    # Check dimensions
    assert "indxm" in result
    assert "indxn" in result
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    assert len(result["indxm"]) == nspec
    assert len(result["indxn"]) == nspec

    # Check values
    assert np.all(result["indxm"] >= 0)
    assert np.all(result["indxn"] >= 0)
    assert np.all(result["indxm"] <= result["indxn"])

    # Check metadata
    assert result["indxm"].attrs["ntrunc"] == ntrunc
    assert result["indxn"].attrs["ntrunc"] == ntrunc


def test_calc_legendre_functions():
    """Test Legendre functions calculation"""
    lat = 45.0
    ntrunc = 10
    result = calc_legendre_functions(lat, ntrunc)

    # Check return type
    assert isinstance(result, xr.DataArray)
    assert result.name == "pnm"

    # Check dimensions
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2
    assert len(result) == nspec

    # Check values are real numbers
    assert np.all(np.isreal(result))


def test_transfer_grid2spectral_transform(sample_grid_data):
    """Test grid to spectral transform"""
    # Test with regular grid
    result = transfer_grid2spectral_transform(
        sample_grid_data, grid_data_type="regular", ntrunc=31
    )

    # Check return type
    assert isinstance(result, xr.DataArray)

    # Check dimensions
    assert "spec_dim" in result.dims
    nspec = (31 + 1) * (31 + 2) // 2
    assert len(result["spec_dim"]) == nspec

    # Check values are complex
    assert np.any(np.iscomplex(result))

    # Check metadata
    assert "Spherical harmonic coefficients" in result.attrs["Description"]

    # Test with gaussian grid
    result_gauss = transfer_grid2spectral_transform(
        sample_grid_data, grid_data_type="gaussian", ntrunc=31
    )
    assert isinstance(result_gauss, xr.DataArray)


def test_transfer_spectral_transform2grid(sample_spectral_data):
    """Test spectral to grid transform"""
    # Test with regular grid
    result = transfer_spectral_transform2grid(
        sample_spectral_data, nlon=128, nlat=64, grid_data_type="regular"
    )

    # Check return type
    assert isinstance(result, xr.DataArray)

    # Check dimensions
    assert "lat" in result.dims
    assert "lon" in result.dims
    assert len(result["lat"]) == 64
    assert len(result["lon"]) == 128

    # Check values are real
    assert np.all(np.isreal(result))

    # Check coordinate values
    assert np.all(result["lon"] >= 0)
    assert np.all(result["lon"] < 360)
    assert np.all(result["lat"] >= -90)
    assert np.all(result["lat"] <= 90)

    # Test with gaussian grid
    result_gauss = transfer_spectral_transform2grid(
        sample_spectral_data, nlon=128, nlat=64, grid_data_type="gaussian"
    )
    assert isinstance(result_gauss, xr.DataArray)


def test_transform_roundtrip(sample_grid_data):
    """Test that grid->spectral->grid transform returns similar data"""
    # Transform to spectral space and back
    spec = transfer_grid2spectral_transform(
        sample_grid_data, grid_data_type="regular", ntrunc=31
    )
    grid = transfer_spectral_transform2grid(
        spec,
        nlon=len(sample_grid_data["lon"]),
        nlat=len(sample_grid_data["lat"]),
        grid_data_type="regular",
    )

    # Check that we got back similar values
    # Note: We can't expect exact equality due to truncation
    assert grid.shape == sample_grid_data.shape
    assert np.allclose(grid.mean(), sample_grid_data.mean(), rtol=0.1)


def test_invalid_inputs():
    """Test error handling for invalid inputs"""
    # Invalid grid type in transform
    with pytest.raises(ValueError):
        transfer_spectral_transform2grid(
            xr.DataArray(np.random.rand(10)), nlon=10, nlat=10, grid_data_type="invalid"
        )
