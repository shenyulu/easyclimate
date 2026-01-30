easyclimate.physics.condensation.lcl
====================================

.. py:module:: easyclimate.physics.condensation.lcl

.. autoapi-nested-parse::

   Lifting Condensation Level (LCL)



Functions
---------

.. autoapisummary::

   easyclimate.physics.condensation.lcl.calc_lifting_condensation_level_bolton1980
   easyclimate.physics.condensation.lcl.calc_lifting_condensation_level_Bohren_Albrecht2023


Module Contents
---------------

.. py:function:: calc_lifting_condensation_level_bolton1980(temperature_data: xarray.DataArray, relative_humidity_data: xarray.DataArray, temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], relative_humidity_data_units: Literal['%', 'dimensionless']) -> xarray.DataArray

   Calculate lifting condensation level using Bolton (1980) approximation.

   .. math::

       T_L = \frac{1}{\frac{1}{T_K - 55} - \frac{\ln (U/100)}{2840}} + 55


   Parameters
   ----------
   temperature_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Atmospheric temperature.
   relative_humidity_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The relative humidity.
   temperature_data_units: :py:class:`str <str>`.
       The unit corresponding to `temperature_data` value. Optional values are ``celsius``, ``kelvin``, ``fahrenheit``.
   relative_humidity_data_units: :py:class:`str <str>`.
       The unit corresponding to `vapor_pressure_data` value. Optional values are ``%``, ``dimensionless``.

   Returns
   -------
   The lifting condensation level ( :math:`\mathrm{K}` ).
       :py:class:`xarray.DataArray<xarray.DataArray>`.

   .. seealso::
       - Bolton, D. (1980). The Computation of Equivalent Potential Temperature. Monthly Weather Review, 108(7), 1046-1053. https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml


.. py:function:: calc_lifting_condensation_level_Bohren_Albrecht2023(pressure_data: xarray.DataArray, temperature_data: xarray.DataArray, dewpoint_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], temperature_data_units: Literal['celsius', 'kelvin', 'fahrenheit'], dewpoint_data_units: Literal['celsius', 'kelvin', 'fahrenheit']) -> xarray.Dataset

   Calculate lifting condensation level using Bohren & Albrecht (2023) approximation.

   According to formulation (6.32) in Bohren & Albrecht (2023),

   .. math::

       T_{LCL}=\frac{1-AT_d}{1/T_d + B\ln(T/T_d)-A}

   Where

   .. math::

       A=-\left(\frac{c_{pv}-c_{pw}}{l_r}-\frac{c_{pd}}{\epsilon l_v}\right),\quad B=\frac{c_{pd}}{\epsilon l_v}

   and

   .. math::

       l_v=l_{vr} - (c_{pv}-c_{pw})T_r, \quad l_v(T_d)=l_r+(c_{pd} - c_{pw}) T_d


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
   The lifting condensation level (:py:class:`xarray.Dataset<xarray.Dataset>`)

   - p_lcl: lifting condensation level pressure ( :math:`\mathrm{hPa}` ).
   - t_lcl: lifting condensation level temperature ( :math:`\mathrm{K}` ).

   .. seealso::
       Bohren, C. F., and B. A. Albrecht, 2023: Atmospheric Thermodynamics Second Edition. Oxford University Press, 579 pp. Website: http://gen.lib.rus.ec/book/index.php?md5=AA3B25841BE3AEBA2628EF9961F58C52


