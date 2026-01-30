"""
Wet-bulb Temperature
"""

from __future__ import annotations
import numpy as np
import xarray as xr
from typing import Literal

from ...core.units import (
    transfer_data_multiple_units,
    transfer_data_temperature_units,
)
from .vapor_pressure import calc_saturation_vapor_pressure
from ..temperature.equivalent_potential_temperature import (
    calc_equivalent_potential_temperature,
)

__all__ = [
    "calc_wet_bulb_temperature_iteration",
    "calc_wet_bulb_potential_temperature_iteration",
    "calc_wet_bulb_potential_temperature_davies_jones2008",
    "calc_wet_bulb_temperature_stull2011",
    "calc_wet_bulb_temperature_sadeghi2013",
]


def calc_wet_bulb_temperature_iteration(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    A: float = 0.662 * 10 ** (-3),
    tolerance: float = 0.01,
    max_iter: int = 100,
    method: Literal["easyclimate-backend", "easyclimate-rust"] = "easyclimate-rust",
) -> xr.DataArray:
    """
    Calculate wet-bulb potential temperature using iteration.

    The iterative formula

    .. math::

        e = e_{tw} - AP(t-t_{w})

    - :math:`e` is the water vapor pressure
    - :math:`e_{tw}` is the saturation water vapor pressure over a pure flat ice surface at wet-bulb temperature :math:`t_w` (when the wet-bulb thermometer is frozen, this becomes the saturation vapor pressure over a pure flat ice surface)
    - :math:`A` is the psychrometer constant
    - :math:`P` is the sea-level pressure
    - :math:`t` is the dry-bulb temperature
    - :math:`t_w` is the wet-bulb temperature

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    A: :py:class:`float <float>`.
        Psychrometer coefficients.

        +-----------------------------------------+---------------------------------+-------------------------------+
        | Psychrometer Type and Ventilation Rate  | Wet Bulb Unfrozen (10^-3/°C^-1) | Wet Bulb Frozen (10^-3/°C^-1) |
        +=========================================+=================================+===============================+
        | Ventilated Psychrometer (2.5 m/s)       | 0.662                           | 0.584                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Spherical Psychrometer (0.4 m/s)        | 0.857                           | 0.756                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Cylindrical Psychrometer (0.4 m/s)      | 0.815                           | 0.719                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Chinese Spherical Psychrometer (0.8 m/s)| 0.7949                          | 0.7949                        |
        +-----------------------------------------+---------------------------------+-------------------------------+

    tolerance: :py:class:`float <float>`.
        Minimum acceptable deviation of the iterated value from the true value.
    max_iter: :py:class:`int <float>`.
        Maximum number of iterations.
    method : {"easyclimate-backend","easyclimate-rust"}
        Backend implementation.


    Returns
    ---------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{degC}` )
        Wet-bulb temperature

    Examples
    ----------------

    .. code:: python

        >>> import xarray as xr
        >>> import numpy as np

        # Create sample data
        >>> temp = xr.DataArray(np.array([20, 25, 30]), dims=['point'])
        >>> rh = xr.DataArray(np.array([50, 60, 70]), dims=['point'])
        >>> pressure = xr.DataArray(np.array([1000, 950, 900]), dims=['point'])

        # Calculate wet-bulb potential temperature
        >>> theta_w = calc_wet_bulb_potential_temperature_iteration(
        ...     temperature_data=temp,
        ...     relative_humidity_data=rh,
        ...     pressure_data=pressure,
        ...     temperature_data_units="celsius",
        ...     relative_humidity_data_units="%",
        ...     pressure_data_units="hPa"
        ... )

        # Example with 2D data
        >>> temp_2d = xr.DataArray(np.random.rand(10, 10) * 30, dims=['lat', 'lon'])
        >>> rh_2d = xr.DataArray(np.random.rand(10, 10) * 100, dims=['lat', 'lon'])
        >>> pres_2d = xr.DataArray(np.random.rand(10, 10) * 200 + 800, dims=['lat', 'lon'])
        >>> theta_w_2d = calc_wet_bulb_potential_temperature_iteration(
        ...     temp_2d, rh_2d, pres_2d, "celsius", "%", "hPa"
        ... )

    .. seealso::
        - Fan, J. (1987). Determination of the Psychrometer Coefficient A of the WMO Reference Psychrometer by Comparison with a Standard Gravimetric Hygrometer. Journal of Atmospheric and Oceanic Technology, 4(1), 239-244. https://journals.ametsoc.org/view/journals/atot/4/1/1520-0426_1987_004_0239_dotpco_2_0_co_2.xml
        - Wang Haijun. (2011). Two Wet-Bulb Temperature Estimation Methods and Error Analysis. Meteorological Monthly (Chinese), 37(4): 497-502. website: http://qxqk.nmc.cn/html/2011/4/20110415.html
        - Cheng Zhi, Wu Biwen, Zhu Baolin, et al, (2011). Wet-Bulb Temperature Looping Iterative Scheme and Its Application. Meteorological Monthly (Chinese), 37(1): 112-115. website: http://qxqk.nmc.cn/html/2011/1/20110115.html

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wet_bulk.py
    """
    from ...backend import _wet_bulb_temperature, calc_wet_bulb_temperature_rs

    # ----------------------------
    # 0) basic validation
    # ----------------------------
    if max_iter <= 0:
        raise ValueError(f"max_iter must be positive, got {max_iter}.")
    if tolerance <= 0:
        raise ValueError(f"tolerance must be positive, got {tolerance}.")
    if A <= 0:
        raise ValueError(f"A must be positive, got {A}.")

    if method not in ("easyclimate-backend", "easyclimate-rust"):
        raise ValueError(
            f"Unknown method={method!r}. "
            "Expected 'easyclimate-backend' or 'easyclimate-rust'."
        )

    # ----------------------------
    # 1) unit conversion to working units
    # ----------------------------
    pressure_hpa = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temp_degc = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    rh_pct = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    # If RH came as dimensionless but transfer function is permissive,
    # enforce expected scale. (If your transfer already handles this, this is harmless.)
    if relative_humidity_data_units == "dimensionless":
        # common convention is [0,1]; convert to %
        rh_pct = rh_pct * 100.0

    # Clip RH to physical range; helps avoid edge-case divergence
    rh_pct = rh_pct.clip(min=0.0, max=100.0)

    # Force float64 everywhere to avoid silent float32 propagation
    temp_degc = temp_degc.astype("float64")
    rh_pct = rh_pct.astype("float64")
    pressure_hpa = pressure_hpa.astype("float64")

    # ----------------------------
    # 2) define scalar kernel (returns scalar, not length-1 array)
    # ----------------------------
    if method == "easyclimate-backend":

        def _kernel(t_dry: np.float64, rh: np.float64, p: np.float64) -> np.float64:
            # NaN propagation
            if np.isnan(t_dry) or np.isnan(rh) or np.isnan(p):
                return np.float64(np.nan)
            return np.float64(
                _wet_bulb_temperature.wet_bulb_temperature(
                    t_dry, rh, p, tolerance, A, max_iter
                )
            )

    else:

        def _kernel(t_dry: np.float64, rh: np.float64, p: np.float64) -> np.float64:
            if np.isnan(t_dry) or np.isnan(rh) or np.isnan(p):
                return np.float64(np.nan)
            return np.float64(
                calc_wet_bulb_temperature_rs(t_dry, rh, p, tolerance, A, max_iter)
            )

    # ----------------------------
    # 3) vectorize over xarray with a true scalar output
    # ----------------------------
    tw = xr.apply_ufunc(
        _kernel,
        temp_degc,
        rh_pct,
        pressure_hpa,
        input_core_dims=[[], [], []],
        output_core_dims=[[]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[np.float64],
        keep_attrs=False,
        dask_gufunc_kwargs={"allow_rechunk": True},
    )

    # ----------------------------
    # 4) metadata
    # ----------------------------
    tw = tw.rename("tw")
    tw.attrs = {
        "standard_name": "wet_bulb_temperature",
        "long_name": "Wet-bulb temperature",
        "units": "degC",
        "method": method,
        "tolerance": float(tolerance),
        "max_iter": int(max_iter),
        "A": float(A),
    }

    return tw


def calc_wet_bulb_potential_temperature_iteration(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    pressure_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit", "degC", "degK"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    A: float = 0.662e-3,
    tolerance: float = 0.01,
    max_iter: int = 100,
    method: Literal["easyclimate-backend", "easyclimate-rust"] = "easyclimate-rust",
) -> xr.DataArray:
    """
    Calculate wet-bulb potential temperature (:math:`\\theta_w`).

    The iterative formula for wet-bulb temperature

    .. math::

        e = e_{tw} - AP(t-t_{w})

    - :math:`e` is the water vapor pressure
    - :math:`e_{tw}` is the saturation water vapor pressure over a pure flat ice surface at wet-bulb temperature :math:`t_w` (when the wet-bulb thermometer is frozen, this becomes the saturation vapor pressure over a pure flat ice surface)
    - :math:`A` is the psychrometer constant
    - :math:`P` is the sea-level pressure
    - :math:`t` is the dry-bulb temperature
    - :math:`t_w` is the wet-bulb temperature

    Wet-bulb potential temperature (:math:`\\theta_w`) is defined as the temperature that an air parcel would have
    if it were first brought to saturation at its ambient pressure (i.e., cooled to the wet-bulb temperature, :math:`T_w`),
    and then brought dry-adiabatically to a reference pressure, conventionally (:math:`p_0 = 1000 \\mathrm{hPa}`).

    This quantity is therefore obtained from two steps:

    - Compute the wet-bulb temperature (:math:`T_w`) at the parcel’s pressure (:math:`p`);
    - Apply the dry-adiabatic (Poisson) transformation from (:math:`p`) to (:math:`p_0`).

    Under this definition, once (:math:`T_w`) is known, (:math:`\\theta_w`) follows directly as

    .. math::

        \\theta_w = (T_w + 273.15) ( \\frac{p_0}{p}) ^{\\kappa} - 273.15

    where :math:`\\kappa = \\frac{R_d}{c_p} \\approx 0.2854`, and :math:`p_0 = 1000 \\mathrm{hPa}`.

    This formulation makes clear that the iterative/nonlinear part of the calculation is confined to determining (:math:`T_w`);
    the mapping from (:math:`T_w`) to (:math:`\\theta_w`) is purely algebraic via the dry-adiabatic relation.

    Parameters
    ----------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    A: :py:class:`float <float>`.
        Psychrometer coefficients.

        +-----------------------------------------+---------------------------------+-------------------------------+
        | Psychrometer Type and Ventilation Rate  | Wet Bulb Unfrozen (10^-3/°C^-1) | Wet Bulb Frozen (10^-3/°C^-1) |
        +=========================================+=================================+===============================+
        | Ventilated Psychrometer (2.5 m/s)       | 0.662                           | 0.584                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Spherical Psychrometer (0.4 m/s)        | 0.857                           | 0.756                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Cylindrical Psychrometer (0.4 m/s)      | 0.815                           | 0.719                         |
        +-----------------------------------------+---------------------------------+-------------------------------+
        | Chinese Spherical Psychrometer (0.8 m/s)| 0.7949                          | 0.7949                        |
        +-----------------------------------------+---------------------------------+-------------------------------+

    tolerance: :py:class:`float <float>`.
        Minimum acceptable deviation of the iterated value from the true value.
    max_iter: :py:class:`int <float>`.
        Maximum number of iterations.
    method : {"easyclimate-backend","easyclimate-rust"}
        Backend implementation.

    Notes
    -----
    :math:`\\theta_w` is obtained by first computing the wet-bulb temperature (Tw)
    and then reducing it dry-adiabatically to 1000 hPa.
    """

    # 1) First compute wet-bulb temperature (Tw)
    tw = calc_wet_bulb_temperature_iteration(
        temperature_data=temperature_data,
        relative_humidity_data=relative_humidity_data,
        pressure_data=pressure_data,
        temperature_data_units=temperature_data_units,
        relative_humidity_data_units=relative_humidity_data_units,
        pressure_data_units=pressure_data_units,
        A=A,
        tolerance=tolerance,
        max_iter=max_iter,
        method=method,
    )

    # 2) Convert pressure to hPa (for consistency)
    p_hpa = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    ).astype("float64")

    # 3) Dry-adiabatic reduction to 1000 hPa
    kappa = 0.2854
    theta_w = (tw + 273.15) * (1000.0 / p_hpa) ** kappa - 273.15

    # 4) metadata
    theta_w = theta_w.rename("theta_w")
    theta_w.attrs = {
        "standard_name": "wet_bulb_potential_temperature",
        "long_name": "Wet-bulb potential temperature",
        "units": "degC",
        "reference_pressure": "1000 hPa",
        "method": method,
    }

    return theta_w


def calc_wet_bulb_potential_temperature_davies_jones2008(
    pressure_data: xr.DataArray,
    temperature_data: xr.DataArray,
    dewpoint_data: xr.DataArray,
    pressure_data_units: Literal["hPa", "Pa", "mbar"],
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    dewpoint_data_units: Literal["celsius", "kelvin", "fahrenheit"],
) -> xr.DataArray:
    """
    Calculate wet-bulb potential temperature using Robert Davies-Jones (2008) approximation.

    Parameters
    ----------
    pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The pressure data set.
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The dewpoint temperature.
    pressure_data_units: :py:class:`str <str>`.
        The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`, `mbar`.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    dewpoint_data_units: :py:class:`str <str>`.
        The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

    Returns
    -------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K}` )
        Wet-bulb temperature

    .. seealso::

        - Davies-Jones, R. (2008). An Efficient and Accurate Method for Computing the Wet-Bulb Temperature along Pseudoadiabats. Monthly Weather Review, 136(7), 2764-2785. https://doi.org/10.1175/2007MWR2224.1
        - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wet_bulk.py
    """
    pressure_data = transfer_data_multiple_units(
        pressure_data, pressure_data_units, "hPa"
    )
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "K"
    )
    dewpoint_data = transfer_data_temperature_units(
        dewpoint_data, dewpoint_data_units, "K"
    )

    # Calculate the equivalent potential temperature
    theta_e = calc_equivalent_potential_temperature(
        pressure_data=pressure_data,
        temperature_data=temperature_data,
        dewpoint_data=dewpoint_data,
        pressure_data_units="hPa",
        temperature_data_units="K",
        dewpoint_data_units="K",
    )  # units: K

    # Create the resulting DataArray, keeping the same coordinates and properties
    theta_w = xr.full_like(theta_e, np.nan)

    # Calculate x value
    x = theta_e / 273.15

    # see (3.8) in Davies-Jones, R. (2008)
    # Create mask condition
    mask = theta_e >= 173.15

    # The points that satisfy the conditions are calculated
    x_masked = x.where(mask)
    x2 = x_masked * x_masked
    x3 = x2 * x_masked
    x4 = x3 * x_masked

    a = 7.101574 - 20.68208 * x_masked + 16.11182 * x2 + 2.574631 * x3 - 5.205688 * x4
    b = 1 - 3.552497 * x_masked + 3.781782 * x2 - 0.6899655 * x3 - 0.5929340 * x4

    # Assign the result to theta_w
    theta_w = xr.where(mask, theta_e - np.exp(a / b), theta_e)

    # clean other attrs
    theta_w.attrs = dict()
    theta_w.name = "tw"
    theta_w.attrs["standard_name"] = "wet_bulb_temperature"
    theta_w.attrs["units"] = "K"
    return theta_w


def calc_wet_bulb_temperature_stull2011(
    temperature_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate wet-bulb temperature using Stull (2011) empirical formula.

    .. math::
        T_{w} =T\\operatorname{atan}[0.151977(\\mathrm{RH} \\% +8.313659)^{1/2}]+\\operatorname{atan}(T+\\mathrm{RH}\\%)-\\operatorname{atan}(\\mathrm{RH} \\% -1.676331)
        +0.00391838(\\mathrm{RH}\\%)^{3/2}\\operatorname{atan}(0.023101\\mathrm{RH}\\%)-4.686035.

    .. tip::

        This methodology was not valid for ambient conditions with low values of :math:`T_a` (dry-bulb temperature; i.e., <10°C),
        and/or with low values of RH  (5% < RH < 10%).
        The Stull methodology was also only valid at sea level.

    Parameters
    ----------------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    ----------------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{K}` )
        Wet-bulb temperature

    .. seealso::
        - Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature. Journal of Applied Meteorology and Climatology, 50(11), 2267-2269. https://doi.org/10.1175/JAMC-D-11-0143.1
        - Stull, R. (2011): Meteorology for Scientists and Engineers. 3rd ed. Discount Textbooks, 924 pp. [Available online at https://www.eoas.ubc.ca/books/Practical_Meteorology/, https://www.eoas.ubc.ca/courses/atsc201/MSE3.html]
        - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wet_bulk.py
    """
    # Convert to Celsius
    t_c = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    # Convert to `%`
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "%"
    )

    tw_c = (
        t_c * np.arctan(0.151977 * np.sqrt(relative_humidity_data + 8.313659))
        + np.arctan(t_c + relative_humidity_data)
        - np.arctan(relative_humidity_data - 1.676331)
        + 0.00391838
        * (relative_humidity_data**1.5)
        * np.arctan(0.023101 * relative_humidity_data)
        - 4.686035
    )

    # Convert back to Kelvin
    tw = transfer_data_temperature_units(tw_c, "degC", "K")

    # clean other attrs
    tw.attrs = dict()
    tw.name = "tw"
    tw.attrs["standard_name"] = "wet_bulb_temperature"
    tw.attrs["units"] = "K"
    return tw


def calc_wet_bulb_temperature_sadeghi2013(
    temperature_data: xr.DataArray,
    height_data: xr.DataArray,
    relative_humidity_data: xr.DataArray,
    temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
    height_data_units: Literal["m", "km"],
    relative_humidity_data_units: Literal["%", "dimensionless"],
) -> xr.DataArray:
    """
    Calculate wet-bulb temperature using Sadeghi et. al (2011) empirical formula.

    Parameters
    ----------------------
    temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        Atmospheric temperature.
    height_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The elevation.
    relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
        The relative humidity.
    temperature_data_units: :py:class:`str <str>`.
        The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
    height_data_units: :py:class:`str <str>`.
        The unit corresponding to `height_data` value. Optional values are `m`, `km`.
    relative_humidity_data_units: :py:class:`str <str>`.
        The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

    Returns
    ----------------------
    tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\\mathrm{degC}` )
        Wet-bulb temperature

    .. seealso::
        - Sadeghi, S., Peters, T. R., Cobos, D. R., Loescher, H. W., & Campbell, C. S. (2013). Direct Calculation of Thermodynamic Wet-Bulb Temperature as a Function of Pressure and Elevation. Journal of Atmospheric and Oceanic Technology, 30(8), 1757-1765. https://doi.org/10.1175/JTECH-D-12-00191.1

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_wet_bulk.py
    """
    gamma = 0.4 * 0.001  # units: g/(kg K) => g/(g K)
    a = 0.611

    # T_a: dry-bulb temperature (K)
    temperature_data = transfer_data_temperature_units(
        temperature_data, temperature_data_units, "degC"
    )
    # H: elevation (m).
    height_data = transfer_data_multiple_units(height_data, height_data_units, "m")
    relative_humidity_data = transfer_data_multiple_units(
        relative_humidity_data, relative_humidity_data_units, "dimensionless"
    )

    # saturation vapor pressure
    e_s_Ta = calc_saturation_vapor_pressure(
        temperature_data=temperature_data, temperature_data_units="degC"
    )  # units: hPa
    # vapor_pressure
    e_a_Ta = e_s_Ta * relative_humidity_data
    e_a = transfer_data_multiple_units(e_a_Ta, "hPa", "kPa")

    # see (11)
    lambda_value = 0.0014 * np.exp(0.027 * temperature_data)
    # see (12)
    zeta = (
        (-3 * 10 ** (-7) * temperature_data**3)
        - (10 ** (-5)) * temperature_data**2
        + (2 * 10 ** (-5)) * temperature_data
        + 4.44 * 10 ** (-2)
    )
    # see (3c)
    p_a = 101.3 * np.exp(-height_data / 8200.0)  # meter for `height_data`
    # see (9c)
    phi = zeta + gamma * p_a
    # see (9a)
    psi = a - gamma * p_a * temperature_data - e_a

    # see (8b)
    delta = phi**2 - 4 * lambda_value * psi

    # see (8a)
    tw = (-phi + np.sqrt(delta)) / (2 * lambda_value)

    # clean other attrs
    tw.attrs = dict()
    tw.name = "tw"
    tw.attrs["standard_name"] = "wet_bulb_temperature"
    tw.attrs["units"] = "degC"
    return tw
