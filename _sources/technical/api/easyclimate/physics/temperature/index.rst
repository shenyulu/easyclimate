easyclimate.physics.temperature
===============================

.. py:module:: easyclimate.physics.temperature


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/physics/temperature/equivalent_potential_temperature/index
   /technical/api/easyclimate/physics/temperature/potential_temperature/index
   /technical/api/easyclimate/physics/temperature/virtual_temperature/index


Functions
---------

.. autoapisummary::

   easyclimate.physics.temperature.calc_equivalent_potential_temperature
   easyclimate.physics.temperature.calc_equivalent_potential_temperature_davies_jones2009
   easyclimate.physics.temperature.calc_potential_temperature
   easyclimate.physics.temperature.calc_potential_temperature_vertical
   easyclimate.physics.temperature.calc_virtual_temperature
   easyclimate.physics.temperature.calc_virtual_temperature_Hobbs2006


Package Contents
----------------

.. py:function:: calc_equivalent_potential_temperature(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate equivalent potential temperature using Bolton (1980) approximation.


   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dew point temperature.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   Equivalent potential temperature ( :math:`\mathrm{K}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml


.. py:function:: calc_equivalent_potential_temperature_davies_jones2009(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.DataArray

   Calculate equivalent potential temperature using Robert Davies-Jones (2009) approximation.


   Parameters
   ----------
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   dewpoint_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The dew point temperature.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.
   dewpoint_data_units: :py:class:`str <str>`.
       The unit corresponding to `dewpoint_data` value. Optional values are `celsius`, `kelvin`, `fahrenheit`.

   Returns
   -------
   Equivalent potential temperature ( :math:`\mathrm{K}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Davies-Jones, R. (2009). On Formulas for Equivalent Potential Temperature. Monthly Weather Review, 137(9), 3137-3148. https://doi.org/10.1175/2009MWR2774.1


.. py:function:: calc_potential_temperature(temper_data: xarray.DataArray, pressure_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], kappa: float = 287 / 1005.7) -> xarray.DataArray

   Calculate the potential temperature for **dry air**.

   Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

   .. math::
       \theta = T \left( \frac{p_0}{p} \right) ^\kappa

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The pressure data set.
   pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `pressure_data` value. Optional values are `hPa`, `Pa`.
   kappa: :py:class:`float <float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   Returns
   -------
   Potential temperature, units according to ``temper_data``.
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   - `Potential temperature - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Potential_temperature>`__
   - `Potential-temperature.pdf <http://weatherclimatelab.mit.edu/wp-content/uploads/2018/02/Potential-temperature.pdf>`__
   - `大气位温、相当位温、饱和相当位温、静力稳定度 <https://renqlsysu.github.io/2019/10/23/potential_temperature/>`__
   - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml

   .. seealso::
       - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
       - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__


.. py:function:: calc_potential_temperature_vertical(temper_data: xarray.DataArray, vertical_dim: str, vertical_dim_units: Literal['hPa', 'Pa', 'mbar'], kappa: float = 287 / 1005.7) -> xarray.DataArray

   Calculate the potential temperature for vertical variables.

   Uses the Poisson equation to calculation the potential temperature given pressure and temperature.

   .. math::
       \theta = T \left( \frac{p_0}{p} \right) ^\kappa

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   vertical_dim: :py:class:`str <str>`.
       Vertical coordinate dimension name.
   vertical_dim_units: :py:class:`str <str>`.
       The unit corresponding to the vertical p-coordinate value. Optional values are `hPa`, `Pa`, `mbar`.
   kappa: :py:class:`float <float>`, default: `287/1005.7`.
       Poisson constant :math:`\kappa`.

       .. note::
           `Poisson constant - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Poisson_constant>`__

   Returns
   -------
   Potential temperature, units according to ``temper_data``.
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   - `Potential temperature - Glossary of Meteorology <https://glossary.ametsoc.org/wiki/Potential_temperature>`__
   - `Potential-temperature.pdf <http://weatherclimatelab.mit.edu/wp-content/uploads/2018/02/Potential-temperature.pdf>`__
   - `大气位温、相当位温、饱和相当位温、静力稳定度 <https://renqlsysu.github.io/2019/10/23/potential_temperature/>`__

   .. seealso::
       - `potential_temperature — MetPy 1.5 <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.potential_temperature.html>`__
       - `pot_temp - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml>`__


.. py:function:: calc_virtual_temperature(temper_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], epsilon: float = 0.608) -> xarray.DataArray

   Calculate virtual temperature.

   The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
   In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
   The virtual temperature is always greater than  the absolute air temperature.

   .. math::
       T_v = T(1+ \epsilon q)

   where :math:`\epsilon = 0.608` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\mathrm{g \cdot g^{-1}}`.

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
   epsilon: :py:class:`float <float>`.
       A constant.

   Returns
   -------
   The virtual temperature, units according to ``temper_data``.
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   - Doswell, C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://journals.ametsoc.org/view/journals/wefo/9/4/1520-0434_1994_009_0625_teontv_2_0_co_2.xml
   - https://en.wikipedia.org/wiki/Virtual_temperature
   - https://glossary.ametsoc.org/wiki/Virtual_temperature

   .. seealso::
       - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
       - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__


.. py:function:: calc_virtual_temperature_Hobbs2006(temper_data: xarray.DataArray, specific_humidity_data: xarray.DataArray, specific_humidity_data_units: Literal['kg/kg', 'g/g', 'g/kg'], epsilon: float = 0.6219569100577033) -> xarray.DataArray

   Calculate virtual temperature.

   The virtual temperature (:math:`T_v`) is the temperature at which dry air would have the same density as the moist air, at a given pressure.
   In other words, two air samples with the same virtual temperature have the same density, regardless of their actual temperature or relative humidity.
   The virtual temperature is always greater than  the absolute air temperature.

   This calculation must be given an air parcel's temperature and mixing ratio. The implementation uses the formula outlined in [Hobbs2006] pg.67 & 80.

   .. math::
       T_v = T \frac{\text{q} + \epsilon}{\epsilon\,(1 + \text{q})}

   where :math:`\epsilon \approx 0.622` when the mixing ratio (specific humidity) :math:`q` is expressed in :math:`\mathrm{g \ g^{-1}}`.

   Parameters
   ----------
   temper_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Air temperature.
   specific_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The absolute humidity data.
   specific_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `specific_humidity` value. Optional values are `kg/kg`, `g/g`, `g/kg` and so on.
   epsilon: :py:class:`float <float>`.
       The molecular weight ratio, which is molecular weight of the constituent gas to that assumed for air. Defaults to the ratio for water vapor to dry air. (:math:`\epsilon \approx 0.622`)

   Returns
   -------
   The virtual temperature, units according to ``temper_data``.
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   Reference
   --------------
   - Hobbs, P. V., and J. M. Wallace, 2006: Atmospheric Science: An Introductory Survey. 2nd ed. Academic Press, 504 pp. https://www.sciencedirect.com/book/9780127329512/atmospheric-science
   - Doswell, C. A., and E. N. Rasmussen, 1994: The Effect of Neglecting the Virtual Temperature Correction on CAPE Calculations. Wea. Forecasting, 9, 625–629, https://journals.ametsoc.org/view/journals/wefo/9/4/1520-0434_1994_009_0625_teontv_2_0_co_2.xml
   - https://en.wikipedia.org/wiki/Virtual_temperature
   - https://glossary.ametsoc.org/wiki/Virtual_temperature

   .. seealso::
       - `virtual_temperature — MetPy <https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.virtual_temperature.html>`__
       - `temp_virtual - NCL <https://www.ncl.ucar.edu/Document/Functions/Contributed/temp_virtual.shtml>`__


