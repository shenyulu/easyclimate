"""
Regridding
"""

from __future__ import annotations

import xarray_regrid
import xarray as xr
import warnings
from ..core.utility import transfer_xarray_lon_from180TO360, generate_dataset_dispatcher

__all__ = ["interp_mesh2mesh"]


@generate_dataset_dispatcher
def interp_mesh2mesh(
    data_input: xr.DataArray | xr.Dataset,
    target_grid: xr.DataArray | xr.Dataset,
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    method: str = "linear",
):
    """
    Regridding regular or lat-lon grid data.

    Parameters
    ----------
    data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
    target_grid: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        Target grid to be regridding.

        :py:class:`xarray.DataArray<xarray.DataArray>` version sample

        .. code:: python

            target_grid = xr.DataArray(
                dims=('lat', 'lon'),
                coords={'lat': np.arange(-89, 89, 3) + 1 / 1.0, 'lon': np.arange(-180, 180, 3) + 1 / 1.0}
            )

        :py:class:`xarray.Dataset<xarray.Dataset>` version sample

        .. code:: python

            target_grid = xr.Dataset()
            target_grid['lat'] = np.arange(-89, 89, 3) + 1 / 1.0
            target_grid['lon'] = np.arange(-180, 180, 3) + 1 / 1.0

    lon_dim: :py:class:`str <str>`, default: `lon`.
        Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
    lat_dim: :py:class:`str <str>`, default: `lat`.
        Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
    method: :py:class:`str <str>`, default: `linear`.
        The methods of regridding.

        - `linear`: linear, bilinear, or higher dimensional linear interpolation.
        - `nearest`: nearest-neighbor regridding.
        - `cubic`: cubic spline regridding.
        - `conservative`: conservative regridding.

    Reference
    --------------
    https://github.com/EXCITED-CO2/xarray-regrid
    """
    target_grid_dims_len = len(target_grid.dims)
    if target_grid_dims_len != 2:
        raise ValueError(
            "The dimension should be 2, rather than %s." % target_grid_dims_len
        )

    # For the convenience of data processing
    target_grid = target_grid.transpose(lat_dim, lon_dim)

    for dims_name in target_grid.dims:
        try:
            data_input[dims_name]
        except Exception as r:
            print(
                "Latitude or Lontitude name should be same between `data_input` and `target_grid`, but here find unknown dimension name: %s."
                % r
            )

    lon_array_data_input = data_input[lon_dim].data
    lon_array_target_grid = target_grid[lon_dim].data

    if (lon_array_data_input < 0).any():
        warnings.warn(
            "It seems that the input data longitude range is from -180° to 180°. Currently automatically converted to it from 0° to 360°."
        )
        data_input = transfer_xarray_lon_from180TO360(data_input)
    if (lon_array_target_grid < 0).any():
        warnings.warn(
            "It seems that the input data longitude range is from -180° to 180°. Currently automatically converted to it from 0° to 360°."
        )
        target_grid = transfer_xarray_lon_from180TO360(target_grid)

    match method:
        case "linear":
            return data_input.regrid.linear(target_grid)
        case "nearest":
            return data_input.regrid.nearest(target_grid)
        case "cubic":
            return data_input.regrid.cubic(target_grid)
        case "conservative":
            return data_input.regrid.conservative(target_grid, latitude_coord=lat_dim)
