"""
Coriolis Force
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from ...core.utility import transfer_deg2rad


__all__ = ["get_coriolis_parameter"]


def get_coriolis_parameter(
    lat_data: xr.DataArray | np.array, omega: float = 7.292e-5
) -> xr.DataArray | np.array:
    """
    Calculate the Coriolis parameter at each point.

    .. math::
        f = 2 \\Omega \\sin(\\phi)

    Parameters
    ----------
    lat_data: :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.
        Latitude at each point.
    omega: :py:class:`float <float>`, default: `7.292e-5` ( :math:`\\mathrm{rad/s}` ).
        The angular speed of the earth.

    Returns
    -------
    Corresponding Coriolis force at each point ( :math:`\\mathrm{s^{-1}}` ).
        :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.

    Reference
    --------------
    - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

    .. seealso::
        - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
        - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_co_coeff.py
    """
    lat_data = lat_data.astype("float64")
    return_data = 2 * omega * np.sin(transfer_deg2rad(lat_data))

    if isinstance(return_data, xr.DataArray):
        return_data.attrs["units"] = "s^-1"
        return_data.name = "coriolis_parameter"
    elif isinstance(return_data, np.ndarray):
        warnings.warn(
            "The unit for the output in `get_coriolis_parameter` is " + "s^-1"
        )
    else:
        raise TypeError("`lat_data` shuold be `xr.DataArray` or `np.array`.")
    return return_data
