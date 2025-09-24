"""
Latent Heat Flux for Water
"""

from __future__ import annotations


__all__ = ["calc_latent_heat_water"]


import numpy as np
import xarray as xr
from typing import Literal
from ...core.utility import transfer_data_temperature_units


def calc_latent_heat_water(
    temperature_data,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    latent_heat_type: Literal[
        "evaporation_condensation", "melting_freezing", "sublimation_deposition"
    ],
) -> xr.DataArray:
    """
    Estimate latent heat flux for water: evaporization (condensation), melting (freezing) or sublimation (deposition).

    .. tip::

        This function returns the latent heat of

        - evaporation/condensation
        - melting/freezing
        - sublimation/deposition

        for water. The latent heatis a function of temperature t. The formulas are polynomial approximations
        to the values in Table 92, p. 343 of the Smithsonian Meteorological Tables, Sixth Revised Edition,
        1963 by Roland List. The approximations were developed by Eric Smith at Colorado State University.

        - Source: Thomas W. Schlatter and Donald V. Baker: PROFS Program Office, NOAA Environmental Research Laboratories, Boulder, Colorado.

    Parameters
    -----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``celsius``, ``kelvin``, ``fahrenheit``.
    latent_heat_type: :py:class:`str <str>`.
        The type of latent heat to estimate. Optional values are ``evaporation_condensation``, ``melting_freezing``, ``sublimation_deposition``.

    Returns
    --------
    Latent heat flux for water ( :math:`\\mathrm{J/kg}` ).
        :py:class:`xarray.DataArray <xarray.DataArray>`.

    .. seealso::

        - https://www.ncl.ucar.edu/Document/Functions/Contributed/latent_heat_water.shtml
    """
    # Polynomial Coefficients
    dtype = np.float64

    if latent_heat_type == "evaporation_condensation":
        a = np.array([3337118.5, -3642.8583, 2.1263947], dtype=dtype)
        lname = "Latent Heat: evaporation/condensation"
    elif latent_heat_type == "melting_freezing":
        a = np.array([-1161004.0, 9002.2648, -12.931292], dtype=dtype)
        lname = "Latent Heat: melting/freezing"
    elif latent_heat_type == "sublimation_deposition":
        a = np.array([2632536.8, 1726.9659, -3.6248111], dtype=dtype)
        lname = "Latent Heat: sublimation/deposition"
    else:
        raise ValueError(f"Illegal key value: latent_heat_type={latent_heat_type}")

    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )  # K

    # Calculate latent heat
    heatl = a[0] + a[1] * temperature_data + a[2] * temperature_data**2  # J/kg

    # clean other attrs
    heatl.attrs = dict()
    heatl.name = "latent_heat_water"
    heatl.attrs["long_name"] = lname
    heatl.attrs["units"] = "J/kg"
    return heatl
