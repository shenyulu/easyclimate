"""
Functions for Weather and climate variable diagnosis.
"""
from __future__ import annotations
import xarray as xr
import numpy as np
from .diff import calc_gradient, calc_p_gradient
from .utility import transfer_deg2rad, transfer_units_coeff


def calc_brunt_vaisala_frequency_atm(
    potential_temperature_data: xr.DataArray, 
    z_data: xr.DataArray, 
    vertical_dim: str, 
    g: float = 9.8
) -> xr.DataArray:
    """
    Calculation of the Brunt-väisälä frequency for the vertical atmosphere.

    .. math::
        N = \\left( \\frac{g}{\\theta} \\frac{\\mathrm{d}\\theta}{\\mathrm{d}z} \\right)^\\frac{1}{2}

    Parameters
    ----------
    potential_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric potential temperature.
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric geopotential height.

    .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    vertical_dim: :py:class:`str<python.str>`.
        Vertical coordinate dimension name.
    g: :py:class:`float<python.float>`, default: `9.8`.
        The acceleration of gravity.

    Returns
    -------
    Brunt-väisälä frequency (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Brunt-väisälä frequency - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Brunt-v%C3%A4is%C3%A4l%C3%A4_frequency>`__

    .. seealso::
        - `brunt_vaisala_frequency — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.brunt_vaisala_frequency.html>`__
        - `brunt_vaisala_atm - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/brunt_vaisala_atm.shtml>`__
        
    """
    dp = 1.0
    dtheta_dp = calc_gradient(potential_temperature_data, dim = vertical_dim) /dp
    dz_dp = calc_gradient(z_data, dim = vertical_dim) /dp
    N = np.sqrt((g /potential_temperature_data) *(dtheta_dp/ dz_dp) )
    return N

def get_coriolis_parameter(
    lat_data: xr.DataArray | np.array, 
    omega: float = 7.292e-5
) -> xr.DataArray | np.array:
    """
    Calculate the Coriolis parameter at each point.

    .. math::
        f = 2 \\Omega \\sin(\\phi)

    Parameters
    ----------
    lat_data: :py:class:`xarray.DataArray <xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`.
        Latitude at each point.
    omega: :py:class:`float<python.float>`, default: `7.292e-5`.
        The angular speed of the earth.

    Returns
    -------
    Corresponding Coriolis force at each point (:py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`numpy.array <numpy.array>`).

    Reference
    --------------
    - `Coriolis parameter - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Coriolis_parameter>`__

    .. seealso::
        - `coriolis_parameter — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.coriolis_parameter.html>`__
        - `coriolis_param - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/coriolis_param.shtml>`__
    """
    lat_data = lat_data.astype('float64')
    return 2 *omega *np.sin(transfer_deg2rad(lat_data))


def get_potential_temperature(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    kappa: float = 287/1005.7
) -> xr.DataArray:
    """
    Calculate the potential temperature.

    Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

    .. math::
        \\theta = T \\left( \\frac{p_0}{p} \\right) ^\\kappa

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str<python.str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str<python.str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    kappa: :py:class:`float<python.float>`, default: `287/1005.7`.
        Poisson constant :math:`\\kappa`.

        .. note::
            `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

    Returns
    -------
    Potential temperature corresponding to the temperature and pressure (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - `Potential temperature - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Potential_temperature>`__
    - `Potential-temperature.pdf <http://weatherclimatelab.mit.edu/wp-content/uploads/2018/02/Potential-temperature.pdf>`__
    - `大气位温、相当位温、饱和相当位温、静力稳定度 <https://renqlsysu.github.io/2019/10/23/potential_temperature/>`__

    .. seealso::
        - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
        - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__
    """
    p_base = transfer_units_coeff(vertical_dim_units, 'Pa')
    if p_base == 1.:
        P_0 = 1000e2
    elif p_base == 100.0:
        P_0 = 1000
    else:
        raise ValueError('`vertical_dim_units` be `Pa`, `hPa`, `mbar`.')
    
    return temper_data *(P_0/ temper_data[vertical_dim]) **(kappa)

def calc_static_stability(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str
) -> xr.DataArray:
    """
    Calculate the static stability within a vertical profile.
    
    .. math::
        \\sigma = - T \\frac{\\partial \\ln \\theta}{\\partial p}

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    vertical_dim: :py:class:`str<python.str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str<python.str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.

    Returns
    -------
    Static stability (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1

    .. seealso::
        - `static_stability - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/static_stability.shtml>`__
        - `static_stability — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.static_stability.html>`__
        - `Static stability parameters · Issue #2535 · Unidata/MetPy <https://github.com/Unidata/MetPy/issues/2535>`__
    """
    theta = get_potential_temperature(temper_data, vertical_dim = vertical_dim, vertical_dim_units = vertical_dim_units)
    ln_theta = np.log(theta)
    part = calc_p_gradient(ln_theta, vertical_dim = vertical_dim, vertical_dim_units = vertical_dim_units)
    return -temper_data *part