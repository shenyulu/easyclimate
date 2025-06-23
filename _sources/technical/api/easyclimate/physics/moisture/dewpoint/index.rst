easyclimate.physics.moisture.dewpoint
=====================================

.. py:module:: easyclimate.physics.moisture.dewpoint

.. autoapi-nested-parse::

   Dewpoint



Functions
---------

.. autoapisummary::

   easyclimate.physics.moisture.dewpoint.calc_dewpoint


Module Contents
---------------

.. py:function:: calc_dewpoint(vapor_pressure_data: xarray.DataArray, vapor_pressure_data_units: Literal['hPa', 'Pa', 'mbar']) -> xarray.DataArray

   Calculate the ambient dewpoint given the vapor pressure.

   Parameters
   ----------
   vapor_pressure_data: :py:class:`xarray.DataArray<xarray.DataArray>`.
       Water vapor partial pressure.
   total_pressure_data_units: :py:class:`str <str>`.
       The unit corresponding to `total_pressure_data` value. Optional values are `hPa`, `Pa`.

   Returns
   -------
   The dew point (:py:class:`xarray.DataArray<xarray.DataArray>`), degrees Celsius.

   .. seealso::
       - https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.dewpoint.html


