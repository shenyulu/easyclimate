"""
Easy climate top interface for the Pyspharm

This is the top layer of packaging for the Pyspharm package.
"""

from ..backend import spharm
import xarray as xr
import numpy as np
from typing import Literal


__all__ = [
    "calc_gaussian_latitudes",
    "calc_geodesic_points",
    "calc_spherical_harmonic_coefficients",
    "calc_legendre_functions",
    "transfer_grid2spectral_transform",
    "transfer_spectral_transform2grid",
]


def calc_gaussian_latitudes(nlat: int):
    """
    Calculate the gaussian latitudes (in degrees) and quadrature weights.

    Parameters
    ----------
    nlat: :py:class:`int <int>`.
        Number of gaussian latitudes desired.

    Returns
    -------
    The gaussian latitudes (in degrees north) and gaussian quadrature weights (:py:class:`xarray.Dataset<xarray.Dataset>`).

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/gaus.shtml
    """
    lats, wts = spharm.gaussian_lats_wts(nlat)

    lats_array = xr.DataArray(lats, name="lats")
    lats_array.attrs["name"] = "Gaussian latitudes (degrees)"
    lats_array.attrs["nlat"] = nlat
    wts_array = xr.DataArray(wts, name="wts")
    wts_array.attrs["name"] = "Quadrature weights"
    wts_array.attrs["nlat"] = nlat

    result = xr.Dataset({"lats": lats_array, "wts": wts_array})
    result.attrs["nlat"] = nlat
    return result


def calc_geodesic_points(m: int):
    """
    Calculate the lat/lon values of the points on the surface of the sphere
    corresponding to a twenty-sided (icosahedral) geodesic.

    Parameters
    ----------
    m: :py:class:`int <int>`.
        The number of points on the edge of a single geodesic triangle.
        There are :math:`10(m-1)^2+2` total geodesic points, including the poles.

    Returns
    -------
    The latitudes and longitudes of the geodesic points (in degrees).
    These points are nearly evenly distributed on the surface of the sphere. (:py:class:`xarray.Dataset<xarray.Dataset>`).
    """
    lats, lons = spharm.getgeodesicpts(m)

    lats_array = xr.DataArray(lats, name="lats")
    lats_array.attrs["name"] = "Geodesic point latitudes (degrees)"
    lats_array.attrs["m"] = m
    lons_array = xr.DataArray(lons, name="lons")
    lons_array.attrs["name"] = "Geodesic point longitudes (degrees)"
    lons_array.attrs["m"] = m

    result = xr.Dataset({"lats": lats_array, "lons": lons_array})
    return result


def calc_spherical_harmonic_coefficients(ntrunc: int):
    """
    Calculate indices of zonal wavenumber (indxm) and degree (indxn) for complex spherical harmonic coefficients.

    Parameters
    ----------
    ntrunc: :py:class:`int <int>`.
        The spherical harmonic triangular truncation limit, i.e, truncation wavenumber (e.g., ``T42``).

    Returns
    -------
    The latitudes and longitudes of the geodesic points (in degrees).
    These points are nearly evenly distributed on the surface of the sphere. (:py:class:`xarray.Dataset<xarray.Dataset>`).
    """
    indxm, indxn = spharm.getspecindx(ntrunc)

    indxm_array = xr.DataArray(indxm, name="indxm")
    indxm_array.attrs["name"] = "Zonal wavenumber of spherical harmonic coefficients"
    indxm_array.attrs["ntrunc"] = ntrunc
    indxn_array = xr.DataArray(indxn, name="indxn")
    indxn_array.attrs["name"] = "Zonal degree of spherical harmonic coefficients"
    indxn_array.attrs["ntrunc"] = ntrunc

    result = xr.Dataset({"indxm": indxm_array, "indxn": indxn_array})
    return result


def calc_legendre_functions(lat: float, ntrunc: int) -> xr.DataArray:
    """
    Calculate associated legendre functions for triangular truncation T(ntrunc), at a given latitude.

    Parameters
    ----------
    lat: :py:class:`float <float>`.
        The latitude (in degrees) to compute the associate legendre functions.
    ntrunc: :py:class:`int <int>`.
        The spherical harmonic triangular truncation limit, i.e, truncation wavenumber (e.g., ``T42``).

    Returns
    -------
    :math:`(\\mathrm{ntrunc} + 1) (\\mathrm{ntrunc} + 2) /2` associated legendre functions at latitude ``lat``.
    """
    pnm = spharm.legendre(lat, ntrunc)
    pnm_array = xr.DataArray(pnm, name="pnm")

    return pnm_array


def transfer_grid2spectral_transform(
    grid_data: xr.DataArray,
    grid_data_type: Literal["regular", "gaussian"],
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    ntrunc: int = None,
) -> xr.DataArray:
    """
    Transform grid data to spectral representation (spherical harmonic analysis).

    Parameters
    ------------
    grid_data: :py:class:`xarray.DataArray <xarray.DataArray>`.
        Input grid data, must contain longitude and latitude dimensions
    grid_data_type:
        Type of grid ('regular' or 'gaussian')
    lon_dim: :py:class:`str <str>`.
        Name of longitude dimension, default is ``'lon'``.
    lat_dim: :py:class:`str <str>`.
        Name of latitude dimension, default is ``'lat'``.
    ntrunc: :py:class:`int <int>`.
        Spectral truncation wavenumber, defaults to ``nlat-1``.

    Returns
    ------------

    :py:class:`xarray.DataArray <xarray.DataArray>` containing complex spherical harmonic coefficients with triangular spectral dimension
    """
    # Get number of latitude and longitude points
    nlat = len(grid_data[lat_dim])
    nlon = len(grid_data[lon_dim])

    # Set default truncation wavenumber if not specified
    if ntrunc is None:
        ntrunc = nlat - 1  # Default truncation wavenumber

    # Calculate number of triangular spectral coefficients
    nspec = (ntrunc + 1) * (ntrunc + 2) // 2

    # Initialize spherical harmonic transform object
    sph = spharm.Spharmt(nlon, nlat, gridtype=grid_data_type)

    def _grdtospec(data):
        """Perform actual grid-to-spectral space transformation"""
        return sph.grdtospec(data, ntrunc=ntrunc)

    # Apply the transformation using xarray's universal function
    result = xr.apply_ufunc(
        _grdtospec,
        grid_data,
        input_core_dims=[[lat_dim, lon_dim]],  # Input dimensions to transform
        output_core_dims=[["spec_dim"]],  # Output spectral dimension
        exclude_dims={lat_dim, lon_dim},  # Dimensions to exclude from output
        vectorize=True,  # Handle multiple arrays if needed
        dask="parallelized",  # Enable parallel processing for dask arrays
    )

    # Add coordinates for the triangular spectral coefficients
    result = result.assign_coords({"spec_dim": np.arange(nspec)})

    # Add metadata attributes
    result.attrs.update(
        {
            "Description": "Spherical harmonic coefficients",
            "Truncation": f"Triangular truncation at wavenumber {ntrunc}",
            "Order": "Complex coefficients in triangular order",
        }
    )

    return result


def transfer_spectral_transform2grid(
    spec_data: xr.DataArray,
    nlon: int,
    nlat: int,
    grid_data_type: Literal["regular", "gaussian"],
    spec_dim: str = "spec_dim",
    lon_dim: str = "lon",
    lat_dim: str = "lat",
) -> xr.DataArray:
    """
    Transform spectral data back to grid space representation (spherical harmonic synthesis)

    Parameters
    ------------
    spec_data: :py:class:`xarray.DataArray <xarray.DataArray>`.
        Input spectral coefficient data, must contain the spectral dimension.
    nlon: :py:class:`int <int>`.
        Number of longitude points in output grid.
    nlat: :py:class:`int <int>`.
        Number of latitude points in output grid.
    grid_data_type:
        Type of output grid (``'regular'`` or ``'gaussian'``).
    spec_dim:
        Name of spectral dimension, default is ``'spec_dim'``.
    lon_dim: :py:class:`str <str>`.
        Name for output longitude dimension, default is ``'lon'``.
    lat_dim: :py:class:`str <str>`.
        Name for output latitude dimension, default is ``'lat'``.

    Returns
    ------------

    :py:class:`xarray.DataArray <xarray.DataArray>` in grid space representation with ``(lat, lon)`` dimensions
    """
    # Initialize spherical harmonic transform object with specified grid parameters
    sph = spharm.Spharmt(nlon, nlat, gridtype=grid_data_type)

    def _spectogrd(data):
        """Perform actual spectral-to-grid space transformation"""
        return sph.spectogrd(data)

    # Apply the inverse transformation using xarray's universal function
    result = xr.apply_ufunc(
        _spectogrd,
        spec_data,
        input_core_dims=[[spec_dim]],  # Input spectral dimension
        output_core_dims=[[lat_dim, lon_dim]],  # Output grid dimensions
        exclude_dims={spec_dim},  # Dimension to exclude from output
        vectorize=True,  # Handle multiple arrays if needed
        dask="parallelized",  # Enable parallel processing for dask arrays
    )

    if grid_data_type == "gaussian":
        # For Gaussian grids, calculate and assign the proper latitude coordinates
        lat_array = calc_gaussian_latitudes(nlat)["lats"].data
        result[lat_dim] = lat_array
    elif grid_data_type == "regular":
        # For regular grids, create equally spaced latitudes from 90°N to 90°S
        result[lat_dim] = np.linspace(90, -90, nlat)
    else:
        # Raise error for unsupported grid types
        raise ValueError(
            f"Unsupported grid type: {grid_data_type}. "
            "Must be either 'regular' or 'gaussian'"
        )

    result[lon_dim] = np.linspace(0, 360, nlon, endpoint=False)
    return result
