"""
Human index
"""

import xarray as xr
import numpy as np
from typing import Literal
from ...core.utility import (
    transfer_data_temperature_units,
    transfer_data_multiple_units,
    transfer_data_difference_units,
)
from ...backend import human_index_mod, human_index_mod_old

__all__ = [
    "calc_apparent_temperature",
    "calc_simplified_human_discomfort_index",
    "calc_simplified_human_discomfort_index_stull",
    "calc_swamp_cooler_temperatures",
    "calc_heat_thic_thip",
    "calc_simplified_wbgt_index",
    "calc_human_feels_temperature",
]


def calc_apparent_temperature(
    temperature_data: xr.DataArray,
    vapor_pressure_data: xr.DataArray,
    wind_10m_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    vapor_pressure_data_units: Literal["hPa", "Pa", "mbar"],
    wind_10m_data_units: Literal["m/s"] = "m/s",
) -> xr.DataArray:
    """
    Calculate apparent temperature.

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vapor pressure.
    wind_10m_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The 10-meter winds.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    vapor_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.
    wind_10m_data_units: :py:class:`str <str>`, default ``m/s``.
        The unit corresponding to `wind_10m_data` value. default value is  ``m/s``.

    Returns
    ---------
    The apparent temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_apptemp.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - http://www.bom.gov.au/info/thermal_stress/#atapproximation
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    - Steadman, R.G., 1994: Norms of apparent temperature in Australia, Aust. Met. Mag., 43, 1-16.
    """
    # conversion unit
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    vapor_pressure_data = transfer_data_multiple_units(
        vapor_pressure_data, vapor_pressure_data_units, "Pa"
    )
    wind_10m_data = transfer_data_multiple_units(
        wind_10m_data, wind_10m_data_units, "m/s"
    )

    # core function
    def _apptemp(t, vp, w10):
        val = human_index_mod.apptemp(t, vp, w10)
        return np.array([val])

    # apply_ufunc
    result = xr.apply_ufunc(
        _apptemp,
        temperature_data,
        vapor_pressure_data,
        wind_10m_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )
    result.attrs["long_name"] = "apparent temperature"
    result.attrs["units"] = "degC"

    return result


def calc_simplified_human_discomfort_index(
    temperature_data: xr.DataArray,
    vapor_pressure_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    vapor_pressure_data_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate a simplified human discomfort index.

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vapor pressure.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    vapor_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

    Returns
    ---------
    The simplified human discomfort index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_discoi.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    - Steadman, R.G., 1994: Norms of apparent temperature in Australia, Aust. Met. Mag., 43, 1-16.
    """
    # conversion unit
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    vapor_pressure_data = transfer_data_multiple_units(
        vapor_pressure_data, vapor_pressure_data_units, "Pa"
    )

    # core function
    def _dis_coi(t, vp):
        val = human_index_mod.dis_coi(t, vp)
        return np.array([val])

    # apply_ufunc
    result = xr.apply_ufunc(
        _dis_coi,
        temperature_data,
        vapor_pressure_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )
    result.attrs["long_name"] = "discomfort index"
    result.attrs["units"] = "degC"

    return result


def calc_simplified_human_discomfort_index_stull(
    temperature_2m_data: xr.DataArray,
    stull_wet_bulb_temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_2m_data_units: Literal["degC", "degF", "degK"],
    stull_wet_bulb_temperature_data_units: Literal["degC", "degF", "degK"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate the human discomfort index due to excessive heat and humidity using the Stull wet bulb temperature.

    Parameters
    ------------
    temperature_2m_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The 2m temperature(s).
    stull_wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The Stull wet bulb temperature (`wetbulb_stull <https://www.ncl.ucar.edu/Document/Functions/Contributed/wetbulb_stull.shtml>`__).
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_2m_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_2m_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    stull_wet_bulb_temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `stull_wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    ---------
    The human discomfort index due to excessive heat and humidity using the Stull wet bulb temperature (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_discoi_stull.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    - Epstein, Y., and D.S. Moran (2006) Thermal comfort and the heat stress indices, Ind. Health, 44, 388-398 doi:https://doi.org/10.2486/indhealth.44.388.
    """
    # conversion unit
    temperature_2m_data = transfer_data_temperature_units(
        temperature_2m_data, temperature_2m_data_units, "degC"
    )
    stull_wet_bulb_temperature_data = transfer_data_temperature_units(
        stull_wet_bulb_temperature_data, stull_wet_bulb_temperature_data_units, "degC"
    )
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    # core function
    def _dis_cois(t, twb, rh):
        val = human_index_mod.dis_cois(t, twb, rh)
        return np.array([val])

    # apply_ufunc
    result = xr.apply_ufunc(
        _dis_cois,
        temperature_2m_data,
        stull_wet_bulb_temperature_data,
        relative_humidity_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )
    result.attrs["long_name"] = "discomfort index with Stull wet bulb temperature"
    result.attrs["units"] = "degC"

    return result


def calc_swamp_cooler_temperatures(
    temperature_data: xr.DataArray,
    wet_bulb_temperature_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    wet_bulb_temperature_data_units: Literal["degC", "degF", "degK"],
) -> xr.DataArray:
    """
    Calculate the swamp cooler temperatures at 65% amd 80% efficiency.

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The wet bulb temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    wet_bulb_temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.

    Returns
    ---------
    The swamp cooler temperatures at 65% amd 80% efficiency (:py:class:`xarray.Dataset<xarray.Dataset>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_swamp_cooleff.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    """
    # conversion unit
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    wet_bulb_temperature_data = transfer_data_temperature_units(
        wet_bulb_temperature_data, wet_bulb_temperature_data_units, "degC"
    )

    # core function
    def _swampcooleff(t, twb):
        val1, val2 = human_index_mod.swampcooleff(t, twb)
        return np.array([val1, val2])

    # apply_ufunc
    tmp_result = xr.apply_ufunc(
        _swampcooleff,
        temperature_data,
        wet_bulb_temperature_data,
        output_core_dims=[["parameter"]],
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 1}, "allow_rechunk": True},
    )
    result1 = tmp_result[..., 0]
    result2 = tmp_result[..., 1]

    result1.attrs["long_name"] = "Swamp Cooler: 80% Efficient"
    result1.attrs["units"] = "degC"
    result2.attrs["long_name"] = "Swamp Cooler: 65% Efficient"
    result2.attrs["units"] = "degC"

    result = xr.Dataset({"swp80": result1, "swp65": result2})
    return result


def calc_heat_thic_thip(
    temperature_data: xr.DataArray,
    wet_bulb_temperature_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    wet_bulb_temperature_data_units: Literal["degC", "degF", "degK"],
) -> xr.DataArray:
    """
    Calculate the thermal humidity comfort index (thic) and the thermal humidity physiology index (thip).

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    wet_bulb_temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The wet bulb temperature.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    wet_bulb_temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `wet_bulb_temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.

    Returns
    ---------
    The thermal humidity comfort index (thic) and the thermal humidity physiology index (thip).

    Quantified estimates for Comfort (THIC) and Physiology (THIP)

    +------------+------------------------+
    |    THIC    |    Description         |
    +============+========================+
    |   75-78    |    alert               |
    +------------+------------------------+
    |   79-83    |    dangerous           |
    +------------+------------------------+
    |    84+     |    very dangerous      |
    +------------+------------------------+

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_thic_thip.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    """
    # conversion unit
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    wet_bulb_temperature_data = transfer_data_temperature_units(
        wet_bulb_temperature_data, wet_bulb_temperature_data_units, "degC"
    )

    # core function
    def _thindex(t, twb):
        val1, val2 = human_index_mod.thindex(t, twb)
        return np.array([val1, val2])

    # apply_ufunc
    tmp_result = xr.apply_ufunc(
        _thindex,
        temperature_data,
        wet_bulb_temperature_data,
        output_core_dims=[["parameter"]],
        vectorize=True,
        dask_gufunc_kwargs={"output_sizes": {"parameter": 1}, "allow_rechunk": True},
    )
    result1 = tmp_result[..., 0]
    result2 = tmp_result[..., 1]

    result1.attrs["long_name"] = "Temperature Humidity Index Comfort: Livestock"
    result1.attrs["units"] = "dimensionless"
    result2.attrs["long_name"] = "Temperature Humidity Index Physiology"
    result2.attrs["units"] = "dimensionless"

    result = xr.Dataset({"thic": result1, "thip": result2})
    return result


def calc_simplified_wbgt_index(
    temperature_data: xr.DataArray,
    vapor_pressure_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    vapor_pressure_data_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate Simplified WBGT index.

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vapor pressure.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    vapor_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

    Returns
    ---------
    The simplified WBGT index (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_wbgt_simplified.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Willett, K.M. and Sherwood, S. (2012), Exceedance of heat index thresholds for 15 regions under a warming climate using the wet-bulb globe temperature. Int. J. Climatol., 32: 161-177. https://doi.org/10.1002/joc.2257
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    """
    # conversion unit
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    vapor_pressure_data = transfer_data_multiple_units(
        vapor_pressure_data, vapor_pressure_data_units, "Pa"
    )

    # core function
    def _swbgt(t, vp):
        val = human_index_mod.swbgt(t, vp)
        return np.array([val])

    # apply_ufunc
    result = xr.apply_ufunc(
        _swbgt,
        temperature_data,
        vapor_pressure_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )
    result.attrs["long_name"] = "simplified wet-bulb index"
    result.attrs["units"] = "dimensionless"

    return result


def calc_human_feels_temperature(
    temperature_data: xr.DataArray,
    vapor_pressure_data: xr.DataArray,
    temperature_data_units: Literal["degC", "degF", "degK"],
    vapor_pressure_data_units: Literal["hPa", "Pa", "mbar"],
) -> xr.DataArray:
    """
    Calculate the 'feels-like' temperature for humans.

    Parameters
    ------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The temperature(s).
    vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The vapor pressure.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are ``degC``, ``degF``, ``degK`` and so on.
    vapor_pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``Pa``, ``hPa``, ``mbar`` and so on.

    Returns
    ---------
    The 'feels-like' temperature for humans (:py:class:`xarray.DataArray<xarray.DataArray>`).

    Reference
    --------------
    - https://www.ncl.ucar.edu/Document/Functions/Heat_stress/heat_humidex.shtml
    - https://github.com/jrbuzan/HumanIndexMod_2020
    - Masterson, J., and F. Richardson, 1979: Humidex, a method of quantifying human discomfort due to excessive heat and humidity CLI 1-79, Environment Canada, Atmosheric Environment Servic website: https://publications.gc.ca/site/eng/9.865813/publication.html, https://publications.gc.ca/collections/collection_2018/eccc/En57-23-1-79-eng.pdf.
    - Buzan, J. R., Oleson, K., and Huber, M.: Implementation and comparison of a suite of heat stress metrics within the Community Land Model version 4.5, Geosci. Model Dev., 8, 151–170, https://doi.org/10.5194/gmd-8-151-2015, 2015.
    """
    # conversion unit
    temperature_data_c = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    vapor_pressure_data = transfer_data_multiple_units(
        vapor_pressure_data, vapor_pressure_data_units, "Pa"
    )

    # core function
    def _hmdex(t_c, vp):
        val = human_index_mod.hmdex(t_c, vp)
        return np.array([val])

    # apply_ufunc
    result = xr.apply_ufunc(
        _hmdex,
        temperature_data_c,
        vapor_pressure_data,
        vectorize=True,
        dask="parallelized",
        dask_gufunc_kwargs={
            "allow_rechunk": True,
        },
    )
    result.attrs["long_name"] = "Humidex"
    result.attrs["description"] = "human feels-like temperature"
    result.attrs["units"] = "dimensionless"

    return result
