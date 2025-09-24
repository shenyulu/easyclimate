"""
Atmospheric Relative Angular Momentum
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from ...core.diff import calc_gradient
from ...core.utility import transfer_data_multiple_units
from typing import Literal

__all__ = ["calc_relative_angular_momentum"]


def calc_relative_angular_momentum(
    zonal_wind_speed_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
    lon_dim: str = "lon",
    lat_dim: str = "lat",
    weights=None,
):
    """
    Calculate atmospheric relative angular momentum.

    Parameters
    -----------
    zonal_wind_speed_data : :py:class:`xarray.DataArray <xarray.DataArray>` ( :math:`\\mathrm{m/s}` )
        Zonal wind component with the least similar dimensions ``(vertical_dim, lon_dim, lat_dim)``
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are ``hPa``, ``Pa``, ``mbar``.
    lon_dim: :py:class:`str <str>`, default: ``lon``.
        Longitude coordinate dimension name. By default extracting is applied over the ``lon`` dimension.
    lat_dim: :py:class:`str <str>`, default: ``lat``.
        Latitude coordinate dimension name. By default extracting is applied over the ``lat`` dimension.
    weights : :py:class:`xarray.DataArray <xarray.DataArray>`, optional
        Weights for each latitude, same dimension as lat.
        If None, computed as :math:`\\cos(lat)*\\mathrm{d}lat` with the values of ``latitude`` spacing.

    Returns
    --------
    aam : :py:class:`xarray.DataArray <xarray.DataArray>` ( :math:`\\mathrm{kg} \\cdot \\mathrm{m^2/s}` )
        Atmospheric angular momentum.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Built-in/angmom_atm.shtml
    """
    # Constants
    pi = np.pi
    rad = pi / 180.0
    twopi = 2.0 * pi
    re = 6.37122e6  # Earth's radius (m)
    re3 = re**3
    g = 9.81  # Gravity (m/s^2)

    def cos_weight(lat):
        """
        Compute latitude weights as cos(lat)*dlat for constant latitude spacing.

        Parameters:
        -----------
        lat : xr.DataArray
            Latitudes (degrees) with dimension (lat)

        Returns:
        --------
        wgt : xr.DataArray
            Weights for each latitude
        """
        rad = np.pi / 180.0
        # Assume constant latitude spacing
        dlatr = np.abs(lat[1] - lat[0]) * rad
        wgt = xr.where(np.abs(lat) == 90.0, 0.0, np.cos(lat * rad) * dlatr)
        return wgt

    # Pressure thickness
    uwnd_level_data = zonal_wind_speed_data[vertical_dim]
    uwnd_level_data = transfer_data_multiple_units(
        uwnd_level_data, vertical_dim_units, "Pa"
    )
    dp = calc_gradient(uwnd_level_data, dim=vertical_dim) * (-1)

    # Latitude
    lat = zonal_wind_speed_data[lat_dim].sortby(lat_dim, ascending=False)

    # Compute weights if not provided
    if weights is None:
        weights = cos_weight(lat)
    else:
        weights = xr.DataArray(weights)

    # Compute u*dp and sum over longitudes and levels
    udp = (zonal_wind_speed_data * dp).sum(dim=[lon_dim, vertical_dim], skipna=True)

    # Zonal average
    udp = udp / zonal_wind_speed_data.sizes[lon_dim]

    # Latitude contribution to aam
    aam = (udp * np.cos(lat * rad) * weights * twopi).sum(dim=lat_dim, skipna=True)

    # Final scaling
    aam = (re3 / g) * aam

    # clean other attrs
    aam.attrs = dict()
    aam.attrs["units"] = "kg m^2/s"
    aam.name = "angular_momentum"
    return aam
