"""
Functions for Weather and climate variable diagnosis.
"""

from __future__ import annotations
import xarray as xr
import numpy as np
from .diff import calc_gradient, calc_p_gradient
from .utility import transfer_deg2rad, transfer_units_coeff, transfer_data_units

__all__ = [
    "calc_brunt_vaisala_frequency_atm",
    "get_coriolis_parameter",
    "calc_potential_temperature",
    "calc_virtual_temperature_Hobbs2006",
    "calc_virtual_temperature",
    "calc_static_stability",
]


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
    z_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Vertical atmospheric geopotential height.

    .. attention:: The unit of `z_data` should be **meters**, NOT :math:`\\mathrm{m^2 \\cdot s^2}` which is the unit used in the representation of potential energy.

    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    g: :py:class:`float <float>`, default: `9.8`.
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
    dtheta_dp = calc_gradient(potential_temperature_data, dim=vertical_dim) / dp
    dz_dp = calc_gradient(z_data, dim=vertical_dim) / dp
    N = np.sqrt((g / potential_temperature_data) * (dtheta_dp / dz_dp))
    return N


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
    omega: :py:class:`float <float>`, default: `7.292e-5`.
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
    lat_data = lat_data.astype("float64")
    return 2 * omega * np.sin(transfer_deg2rad(lat_data))


def calc_potential_temperature(
    temper_data: xr.DataArray,
    vertical_dim: str,
    vertical_dim_units: str,
    kappa: float = 287 / 1005.7,
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
    vertical_dim: :py:class:`str <str>`.
        Vertical coordinate dimension name.
    vertical_dim_units: :py:class:`str <str>`.
        The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
    kappa: :py:class:`float <float>`, default: `287/1005.7`.
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
    p_base = transfer_units_coeff(vertical_dim_units, "Pa")
    if p_base == 1.0:
        P_0 = 1000e2
    elif p_base == 100.0:
        P_0 = 1000
    else:
        raise ValueError("`vertical_dim_units` be `Pa`, `hPa`, `mbar`.")

    return temper_data * (P_0 / temper_data[vertical_dim]) ** (kappa)


def calc_virtual_temperature(
    temper_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    specific_humidity_units: str,
    epsilon=0.608,
) -> xr.DataArray:
    """
    Calculate virtual temperature.

    The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
    In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
    The virtual temperature is always greater than  the absolute air temperature.

    .. math::
        T_v = T(1+ \\epsilon q)

    where :math:`\\epsilon = 0.608` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\\mathrm{g \ g^{-1}}`.

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    epsilon: :py:class:`float <float>`.
        A constant.

    Reference
    --------------
    - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
    - https://en.wikipedia.org/wiki/Virtual_temperature
    - https://glossary.ametsoc.org/wiki/Virtual_temperature

    .. seealso::
        - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
        - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__
    """
    specific_humidity_data = transfer_data_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )

    T_v = (1 + epsilon * specific_humidity_data) * temper_data
    T_v.name = "virtual_temperature"
    return T_v


def calc_virtual_temperature_Hobbs2006(
    temper_data: xr.DataArray,
    specific_humidity_data: xr.DataArray,
    specific_humidity_units: str,
    epsilon=0.6219569100577033,
) -> xr.DataArray:
    """
    Calculate virtual temperature.

    The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
    In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
    The virtual temperature is always greater than  the absolute air temperature.

    This calculation must be given an air parcel's temperature and mixing ratio. The implementation uses the formula outlined in [Hobbs2006] pg.67 & 80.

    .. math::
        T_v = T \\frac{\\text{q} + \\epsilon}{\\epsilon\\,(1 + \\text{q})}

    where :math:`\\epsilon \\approx 0.622` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\\mathrm{g \ g^{-1}}`.

    Parameters
    ----------
    temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Air temperature.
    specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The absolute humidity data.
    specific_humidity_units: :py:class:`str <str>`.
        The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/kg` and so on.
    epsilon: :py:class:`float <float>`.
        The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air. Defaults to the ratio for water vapor to dry air. (:math:`\\epsilon \\approx 0.622`)

    Reference
    --------------
    - Hobbs, P. V., and J. M. Wallace, 2006: Atmospheric Science: An Introductory Survey. 2nd ed. Academic Press, 504 pp. https://www.sciencedirect.com/book/9780127329512/atmospheric-science
    - Doswell , C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2.
    - https://en.wikipedia.org/wiki/Virtual_temperature
    - https://glossary.ametsoc.org/wiki/Virtual_temperature

    .. seealso::
        - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
        - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__
    """
    specific_humidity_data = transfer_data_units(
        specific_humidity_data, specific_humidity_units, "g/g"
    )

    T_v = temper_data * (
        (specific_humidity_data + epsilon) / (epsilon * (1 + specific_humidity_data))
    )
    T_v.name = "virtual_temperature"
    return T_v


def calc_static_stability(
    temper_data: xr.DataArray, vertical_dim: str, vertical_dim_units: str
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
    Static stability (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    Howard B. Bluestein. (1992). Synoptic-Dynamic Meteorology in Midlatitudes: Principles of Kinematics and Dynamics, Vol. 1

    .. seealso::
        - `static_stability - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/static_stability.shtml>`__
        - `static_stability — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.static_stability.html>`__
        - `Static stability parameters · Issue #2535 · Unidata/MetPy <https://github.com/Unidata/MetPy/issues/2535>`__
    """
    theta = calc_potential_temperature(
        temper_data, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    ln_theta = np.log(theta)
    part = calc_p_gradient(
        ln_theta, vertical_dim=vertical_dim, vertical_dim_units=vertical_dim_units
    )
    return -temper_data * part
