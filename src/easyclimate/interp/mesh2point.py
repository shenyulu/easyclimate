"""
Interpolating mesh/grid data to point locations

This module provides functionality to interpolate values from a regular grid/mesh
(e.g., climate model output, remote sensing data) to specific point locations
(e.g., weather stations, observation points) using various interpolation methods.
"""

from __future__ import annotations

import xarray as xr
import pandas as pd
import numpy as np
from scipy.interpolate import RegularGridInterpolator

__all__ = ["interp_mesh2point"]

from typing import Literal


def interp_mesh2point(
    data_input: xr.DataArray,
    df: pd.DataFrame,
    lon_dim_mesh: str = "lon",
    lat_dim_mesh: str = "lat",
    lon_dim_df: str = "lon",
    lat_dim_df: str = "lat",
    method: Literal[
        "linear", "nearest", "slinear", "cubic", "quintic", "pchip"
    ] = "linear",
):
    """
    Interpolate values from a regular grid/mesh to specific point locations.

    Parameters
    ----------
    data_input : xr.DataArray
        Input grid data with latitude and longitude dimensions.
    df : pd.DataFrame
        DataFrame containing point locations with latitude and longitude columns.
    lon_dim_mesh : str, optional
        Name of the longitude dimension in the input grid, by default ``"lon"``.
    lat_dim_mesh : str, optional
        Name of the latitude dimension in the input grid, by default ``"lat"``.
    lon_dim_df : str, optional
        Name of the longitude column in the DataFrame, by default ``"lon"``.
    lat_dim_df : str, optional
        Name of the latitude column in the DataFrame, by default ``"lat"``.
    method : str, optional
        Interpolation method to use. Options are:

        - ``"linear"``: bilinear interpolation (default)
        - ``"nearest"``: nearest neighbor
        - ``"slinear"``: spline linear
        - ``"cubic"``: spline cubic
        - ``"quintic"``: spline quintic
        - ``"pchip"``: piecewise cubic Hermite interpolating polynomial

    Returns
    -------
    pd.DataFrame
        Original DataFrame with an additional column "interpolated_value" containing
        the interpolated values. Points outside the grid range will have NaN values.

    Raises
    ------
    ValueError
        If no points in the DataFrame fall within the grid's spatial extent.

    Examples
    --------
    >>> import xarray as xr
    >>> import pandas as pd
    >>> # Create sample grid data
    >>> lats = np.linspace(-90, 90, 181)
    >>> lons = np.linspace(-180, 180, 361)
    >>> data = xr.DataArray(np.random.rand(181, 361), dims=['lat', 'lon'],
    ...                     coords={'lat': lats, 'lon': lons})
    >>> # Create sample points
    >>> points = pd.DataFrame({'lat': [45.5, 30.2], 'lon': [-120.3, 150.7]})
    >>> # Interpolate
    >>> result = interp_mesh2point(data, points)

    .. seealso::

        - https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
        - https://www.ncl.ucar.edu/Applications/station.shtml
    """
    interp_func = RegularGridInterpolator(
        (data_input[lat_dim_mesh].values, data_input[lon_dim_mesh].values),
        data_input.transpose(lat_dim_mesh, lon_dim_mesh).values,
        method=method,
    )

    # Get grid bounds
    lat_min, lat_max = (
        data_input[lat_dim_mesh].min().data,
        data_input[lat_dim_mesh].max().data,
    )
    lon_min, lon_max = (
        data_input[lon_dim_mesh].min().data,
        data_input[lon_dim_mesh].max().data,
    )

    # Filter stations
    valid_mask = (
        (df[lat_dim_df] >= lat_min)
        & (df[lat_dim_df] <= lat_max)
        & (df[lon_dim_df] >= lon_min)
        & (df[lon_dim_df] <= lon_max)
    )
    valid_stations = df[valid_mask]

    # Check if any valid stations exist
    if len(valid_stations) == 0:
        raise ValueError("No stations fall within the grid data extent")

    # Perform interpolation
    station_values = interp_func(
        (valid_stations[lat_dim_df].values, valid_stations[lon_dim_df].values)
    )

    # Add results back to original DataFrame
    df = df.copy()
    df["interpolated_value"] = np.nan  # Create column filled with NaN
    df.loc[valid_mask, "interpolated_value"] = station_values
    return df
