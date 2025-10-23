easyclimate.interp.mesh2point
=============================

.. py:module:: easyclimate.interp.mesh2point

.. autoapi-nested-parse::

   Interpolating mesh/grid data to point locations

   This module provides functionality to interpolate values from a regular grid/mesh
   (e.g., climate model output, remote sensing data) to specific point locations
   (e.g., weather stations, observation points) using various interpolation methods.



Functions
---------

.. autoapisummary::

   easyclimate.interp.mesh2point.interp_mesh2point
   easyclimate.interp.mesh2point.interp_mesh2point_withtime


Module Contents
---------------

.. py:function:: interp_mesh2point(data_input: xarray.DataArray, df: pandas.DataFrame, lon_dim_mesh: str = 'lon', lat_dim_mesh: str = 'lat', lon_dim_df: str = 'lon', lat_dim_df: str = 'lat', method: Literal['linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'] = 'linear')

   Interpolate values from a regular grid/mesh to specific point locations.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input grid data with latitude and longitude dimensions.
   df : :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       DataFrame containing point locations with latitude and longitude columns.
   lon_dim_mesh : :py:class:`str <str>`, optional
       Name of the longitude dimension in the input grid, by default ``"lon"``.
   lat_dim_mesh : :py:class:`str <str>`, optional
       Name of the latitude dimension in the input grid, by default ``"lat"``.
   lon_dim_df : :py:class:`str <str>`, optional
       Name of the longitude column in the DataFrame, by default ``"lon"``.
   lat_dim_df : :py:class:`str <str>`, optional
       Name of the latitude column in the DataFrame, by default ``"lat"``.
   method : :py:class:`str <str>`, optional
       Interpolation method to use. Options are:

       - ``"linear"``: bilinear interpolation (default)
       - ``"nearest"``: nearest neighbor
       - ``"slinear"``: spline linear
       - ``"cubic"``: spline cubic
       - ``"quintic"``: spline quintic
       - ``"pchip"``: piecewise cubic Hermite interpolating polynomial

   Returns
   -------
   :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       Original DataFrame with an additional column ``"interpolated_value"`` containing
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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_interp_mesh2point.py


.. py:function:: interp_mesh2point_withtime(data_input: xarray.DataArray, stations_df: pandas.DataFrame, lon_dim_df: str, lat_dim_df: str, station_dim_df: str, lon_dim_mesh: str = 'lon', lat_dim_mesh: str = 'lat', time_dim_mesh: str = 'time', method: Literal['linear', 'nearest', 'slinear', 'cubic', 'quintic', 'pchip'] = 'linear') -> xarray.DataArray

   Interpolate gridded data to specific point locations with time series preservation.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input grid data with time, latitude and longitude dimensions.
   stations_df : :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       DataFrame containing station locations with ID, latitude and longitude.
   lon_dim_df : :py:class:`str <str>`.
       Name of the longitude column in the DataFrame.
   lat_dim_df : :py:class:`str <str>`.
       Name of the latitude column in the DataFrame.
   station_dim_df : :py:class:`str <str>`.
       Name of the station ID column in the DataFrame.
   lon_dim_mesh : :py:class:`str <str>`, default: "lon"
       Name of the longitude dimension in the input grid.
   lat_dim_mesh : :py:class:`str <str>`, default: "lat"
       Name of the latitude dimension in the input grid.
   time_dim_mesh : :py:class:`str <str>`, default: "time"
       Name of the time dimension in the input grid.
   method : :py:class:`str <str>`, optional
       Interpolation method:

       - "linear": bilinear interpolation (default)
       - "nearest": nearest neighbor
       - "slinear": spline linear
       - "cubic": spline cubic
       - "quintic": spline quintic
       - "pchip": piecewise cubic Hermite interpolating polynomial

   Returns
   -------
   :py:class:`xarray.Dataset <xarray.Dataset>`.
       :py:class:`xarray.Dataset <xarray.Dataset>` with dimensions ``(station, time)`` containing:

       - Interpolated values
       - Station coordinates


   Examples
   --------
   >>> import xarray as xr
   >>> import pandas as pd
   >>> import numpy as np
   >>> # Create sample grid
   >>> times = pd.date_range("2020-01-01", periods=5)
   >>> lats = np.linspace(-45, -10, 100)
   >>> lons = np.linspace(110, 156, 120)
   >>> data = xr.DataArray(
   ...     np.random.rand(5, 100, 120),
   ...     dims=["time", "lat", "lon"],
   ...     coords={"time": times, "lat": lats, "lon": lons}
   ... )
   >>> # Create stations
   >>> stations = pd.DataFrame({
   ...     "station_id_col": [1001, 1002],
   ...     "lat_col": [-15.5, -20.3],
   ...     "lon_col": [125.5, 130.2]
   ... })
   >>> # Interpolate
   >>> result = interp_mesh2point_withtime(
   ...     data, stations,
   ...     lon_dim_df="lon_col", lat_dim_df="lat_col", station_dim_df="station_id_col"
   ... )

   .. seealso::

       :py:class:`scipy.interpolate.RegularGridInterpolator <scipy.interpolate.RegularGridInterpolator>`.


