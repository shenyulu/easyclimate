"""
Atmospheric Stability
"""

from __future__ import annotations
import numpy as np
import xarray as xr
import warnings
from typing import Literal
from ...core.diff import calc_gradient, calc_p_gradient
from ..temperature.potential_temperature import calc_potential_temperature_vertical

__all__ = ["calc_brunt_vaisala_frequency_atm", "calc_static_stability"]


def calc_brunt_vaisala_frequency_atm(
    potential_temperature_data: xr.DataArray,
    z_data: xr.DataArray,
    vertical_dim: str,
    g: float = 9.8,
) -> xr.DataArray:
    """
    Calculation of the Brunt-väisälä frequency for the vertical atmosphere.

    .. math::

        N = \\left( \\frac{g}{\\theta} \\frac{\\mathrm{d}\\theta}{\\mathrm{d}z} \\right)^\\frac{1}{2}

    Parameters
    ----------
    potential_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric potential temperature.
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{m}` ).
        Vertical atmospheric geopotential height.

    .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    g: :py:class:`float <float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    Brunt-väisälä frequency, units according to ``potential_temperature_data`` :math:`^{1/2}`.
        :py:class:`xarray.DataArray<xarray.DataArray>`

    Reference
    --------------
    - `Brunt-väisälä frequency - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Brunt-v%C3%A4is%C3%A4l%C3%A4_frequency>`__

    .. seealso::
        - `brunt_vaisala_frequency — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.brunt_vaisala_frequency.html>`__
        - `brunt_vaisala_atm - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/brunt_vaisala_atm.shtml>`__

    """
    dp = 1.0
    dtheta_dp = calc_gradient(potential_temperature_data, dim=vertical_dim) / dp
    dz_dp = calc_gradient(z_data, dim=vertical_dim) / dp
    N = np.sqrt((g / potential_temperature_data) * (dtheta_dp / dz_dp))

    # clean other attrs
    N.attrs = dict()
    try:
        z_data_units = z_data.attrs["units"]
        N.attrs["units"] = f"{z_data_units}^(1/2)"
    except:
        N.attrs["units"] = f"[L]^(1/2)"
        warnings.warn(
            "The variable of `z_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    N.name = "brunt_vaisala_frequency_atm"
    return N


def calc_static_stability(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate the static stability within a vertical profile.

    .. math::

        \\sigma = - T \\frac{\\partial \\ln \\theta}{\\partial p}

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    Static stability, units according to ``temper_data_units^2 vertical_dim_units^-1``.
        :py:class:`xarray.DataArray<xarray.DataArray>`.

    Reference
    --------------
    - Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1

    .. seealso::
        - `static_stability - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/static_stability.shtml>`__
        - `static_stability — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.static_stability.html>`__
        - `Static stability parameters · Issue #2535 · Unidata/MetPy <https://github.com/Unidata/MetPy/issues/2535>`__
    """
    theta = calc_potential_temperature_vertical(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    ln_theta = np.log(theta)
    part = calc_p_gradient(
        ln_theta, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    return_data = -temper_data * part

    # clean other attrs
    return_data.attrs = dict()
    try:
        temper_data_units = temper_data.attrs["units"]
        return_data.attrs["units"] = f"{temper_data_units}^2 {vertical_dim_units}^-1"
    except:
        return_data.attrs["units"] = f"[T]^2 {vertical_dim_units}^-1"
        warnings.warn(
            "The variable of `temper_data` do not have `units` attribution, so the attribution of units in it is assigned to the symbol for dimension!"
        )
    return_data.name = "static_stability"
    return return_data
