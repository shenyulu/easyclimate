easyclimate.physics.moisture.wet_bulb
=====================================

.. py:module:: easyclimate.physics.moisture.wet_bulb

.. autoapi-nested-parse::

   Wet-bulb Temperature



Functions
---------

.. autoapisummary::

   easyclimate.physics.moisture.wet_bulb.calc_wet_bulb_temperature_iteration
   easyclimate.physics.moisture.wet_bulb.calc_wet_bulb_potential_temperature_iteration
   easyclimate.physics.moisture.wet_bulb.calc_wet_bulb_potential_temperature_davies_jones2008
   easyclimate.physics.moisture.wet_bulb.calc_wet_bulb_temperature_stull2011
   easyclimate.physics.moisture.wet_bulb.calc_wet_bulb_temperature_sadeghi2013


Module Contents
---------------

.. py:function:: calc_wet_bulb_temperature_iteration(temperature_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, pressure_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], relative_humidity_data_units: Literal['%', 'dimensionless'], pressure_data_units: Literal['hPa', 'Pa', 'mbar'], A: float = 0.662 * 10**(-3), tolerance: float = 0.01, max_iter: int = 100, method: Literal['easyclimate-backend', 'easyclimate-rust'] = 'easyclimate-rust') -> xarray.DataArray

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
   tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\mathrm{degC}` )
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


.. py:function:: calc_wet_bulb_potential_temperature_iteration(temperature_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, pressure_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit', 'degC', 'degK'], relative_humidity_data_units: Literal['%', 'dimensionless'], pressure_data_units: Literal['hPa', 'Pa', 'mbar'], A: float = 0.000662, tolerance: float = 0.01, max_iter: int = 100, method: Literal['easyclimate-backend', 'easyclimate-rust'] = 'easyclimate-rust') -> xarray.DataArray

   Calculate wet-bulb potential temperature (:math:`\theta_w`).

   The iterative formula for wet-bulb temperature

   .. math::

       e = e_{tw} - AP(t-t_{w})

   - :math:`e` is the water vapor pressure
   - :math:`e_{tw}` is the saturation water vapor pressure over a pure flat ice surface at wet-bulb temperature :math:`t_w` (when the wet-bulb thermometer is frozen, this becomes the saturation vapor pressure over a pure flat ice surface)
   - :math:`A` is the psychrometer constant
   - :math:`P` is the sea-level pressure
   - :math:`t` is the dry-bulb temperature
   - :math:`t_w` is the wet-bulb temperature

   Wet-bulb potential temperature (:math:`\theta_w`) is defined as the temperature that an air parcel would have
   if it were first brought to saturation at its ambient pressure (i.e., cooled to the wet-bulb temperature, :math:`T_w`),
   and then brought dry-adiabatically to a reference pressure, conventionally (:math:`p_0 = 1000 \mathrm{hPa}`).

   This quantity is therefore obtained from two steps:

   - Compute the wet-bulb temperature (:math:`T_w`) at the parcel’s pressure (:math:`p`);
   - Apply the dry-adiabatic (Poisson) transformation from (:math:`p`) to (:math:`p_0`).

   Under this definition, once (:math:`T_w`) is known, (:math:`\theta_w`) follows directly as

   .. math::

       \theta_w = (T_w + 273.15) ( \frac{p_0}{p}) ^{\kappa} - 273.15

   where :math:`\kappa = \frac{R_d}{c_p} \approx 0.2854`, and :math:`p_0 = 1000 \mathrm{hPa}`.

   This formulation makes clear that the iterative/nonlinear part of the calculation is confined to determining (:math:`T_w`);
   the mapping from (:math:`T_w`) to (:math:`\theta_w`) is purely algebraic via the dry-adiabatic relation.

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
   :math:`\theta_w` is obtained by first computing the wet-bulb temperature (Tw)
   and then reducing it dry-adiabatically to 1000 hPa.


.. py:function:: calc_wet_bulb_potential_temperature_davies_jones2008(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

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
   tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\mathrm{K}` )
       Wet-bulb temperature

   .. seealso::

       - Davies-Jones, R. (2008). An Efficient and Accurate Method for Computing the Wet-Bulb Temperature along Pseudoadiabats. Monthly Weather Review, 136(7), 2764-2785. https://doi.org/10.1175/2007MWR2224.1
       - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wet_bulk.py


.. py:function:: calc_wet_bulb_temperature_stull2011(temperature_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], relative_humidity_data_units: Literal['%', 'dimensionless']) -> xarray.DataArray

   Calculate wet-bulb temperature using Stull (2011) empirical formula.

   .. math::
       T_{w} =T\operatorname{atan}[0.151977(\mathrm{RH} \% +8.313659)^{1/2}]+\operatorname{atan}(T+\mathrm{RH}\%)-\operatorname{atan}(\mathrm{RH} \% -1.676331)
       +0.00391838(\mathrm{RH}\%)^{3/2}\operatorname{atan}(0.023101\mathrm{RH}\%)-4.686035.

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
   tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\mathrm{K}` )
       Wet-bulb temperature

   .. seealso::
       - Stull, R. (2011). Wet-Bulb Temperature from Relative Humidity and Air Temperature. Journal of Applied Meteorology and Climatology, 50(11), 2267-2269. https://doi.org/10.1175/JAMC-D-11-0143.1
       - Stull, R. (2011): Meteorology for Scientists and Engineers. 3rd ed. Discount Textbooks, 924 pp. [Available online at https://www.eoas.ubc.ca/books/Practical_Meteorology/, https://www.eoas.ubc.ca/courses/atsc201/MSE3.html]
       - Knox, J. A., Nevius, D. S., & Knox, P. N. (2017). Two Simple and Accurate Approximations for Wet-Bulb Temperature in Moist Conditions, with Forecasting Applications. Bulletin of the American Meteorological Society, 98(9), 1897-1906. https://doi.org/10.1175/BAMS-D-16-0246.1

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wet_bulk.py


.. py:function:: calc_wet_bulb_temperature_sadeghi2013(temperature_data: xarray.DataArray, height_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], height_data_units: Literal['m', 'km'], relative_humidity_data_units: Literal['%', 'dimensionless']) -> xarray.DataArray

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
   tw: :py:class:`xarray.DataArray<xarray.DataArray>` ( :math:`\mathrm{degC}` )
       Wet-bulb temperature

   .. seealso::
       - Sadeghi, S., Peters, T. R., Cobos, D. R., Loescher, H. W., & Campbell, C. S. (2013). Direct Calculation of Thermodynamic Wet-Bulb Temperature as a Function of Pressure and Elevation. Journal of Atmospheric and Oceanic Technology, 30(8), 1757-1765. https://doi.org/10.1175/JTECH-D-12-00191.1

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_wet_bulk.py


