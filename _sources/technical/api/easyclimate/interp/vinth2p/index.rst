easyclimate.interp.vinth2p
==========================

.. py:module:: easyclimate.interp.vinth2p

.. autoapi-nested-parse::

   Interpolates Community Atmosphere Model (CAM) or Community Earth System Model (CESM) hybrid coordinates to pressure coordinates

   .. seealso::

       - https://ncar.github.io/CAM/doc/build/html/
       - https://www.cesm.ucar.edu/models/cam
       - https://www.cesm.ucar.edu/



Functions
---------

.. autoapisummary::

   easyclimate.interp.vinth2p.interp_vinth2p_dp
   easyclimate.interp.vinth2p.interp_vinth2p_ecmwf
   easyclimate.interp.vinth2p.interp_vintp2p_ecmwf


Module Contents
---------------

.. py:function:: interp_vinth2p_dp(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], hybrid_A_coefficients: xarray.DataArray, hybrid_B_coefficients: xarray.DataArray, vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', p0_hPa=1000.0) -> xarray.DataArray

   Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure coordinates to pressure levels.

   This function performs vertical interpolation of atmospheric data (typically temperature)
   from hybrid sigma-pressure coordinates to specified pressure levels using the NCL's
   ``vinth2p`` algorithm implemented in Fortran.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Input 3D temperature field on hybrid levels with dimensions, e.g., (time, lev, lat, lon).
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the surface pressure data.
   hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Hybrid A coefficients (pressure term) for the model levels.
   hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`.
       Hybrid B coefficients (sigma term) for the model levels.
   vertical_output_level : list[:py:class:`int <in>`t | :py:class:`float <float>`]
       List of target pressure levels for interpolation.
   vertical_input_dim : :py:class:`str <str>`.
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`.
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`.
       Units for the output pressure levels (must be convertible to hPa).
   interp_method : Literal["linear", "log", "loglog"], optional
       Interpolation method:

       - ``"linear"``: Linear interpolation
       - ``"log"``: Logarithmic interpolation
       - ``"loglog"``: Log-log interpolation

       Default is ``"linear"``.
   extrapolation : :py:class:`bool <bool>`., optional
       Whether to extrapolate below the lowest model level when needed.
       Default is ``False``.
   lon_dim : :py:class:`str <str>`., optional
       Name of the longitude dimension. Default is ``"lon"``.
   lat_dim : :py:class:`str <str>`., optional
       Name of the latitude dimension. Default is ``"lat"``.
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is ``"time"``.
   p0_hPa : :py:class:`float <float>`., optional
       Reference pressure in hPa for hybrid level calculation. Default is ``1000.0`` hPa.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to vertical_output_level.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels
   >>> interp_vinth2p_dp(
   ...     data_input=temp_data,
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     hybrid_A_coefficients=hyam,
   ...     hybrid_B_coefficients=hybm,
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     interp_method="log"
   ... )


.. py:function:: interp_vinth2p_ecmwf(data_input: xarray.DataArray, surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], hybrid_A_coefficients: xarray.DataArray, hybrid_B_coefficients: xarray.DataArray, vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, variable_flag: Literal['T', 'Z', 'other'], temperature_bottom_data: Optional[xarray.DataArray] = None, surface_geopotential_data: Optional[xarray.DataArray] = None, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = True, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time', p0_hPa: float = 1000.0) -> xarray.DataArray

   Interpolate atmospheric data from Community Atmosphere Model (CAM) hybrid sigma-pressure
   coordinates to pressure levels using ECMWF extrapolation methods.

   This function performs vertical interpolation of atmospheric data from hybrid sigma-pressure
   coordinates to specified pressure levels using the NCL's ``vinth2p_ecmwf`` algorithm implemented
   in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
   height ('Z') below the lowest hybrid level.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input 3D field (e.g., temperature or geopotential) on hybrid levels with dimensions,
       e.g., (time, lev, lat, lon).
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : ``Literal["hPa", "Pa", "mbar"]``
       Units of the surface pressure data.
   hybrid_A_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
       Hybrid A coefficients (pressure term) for the model levels.
   hybrid_B_coefficients : :py:class:`xarray.DataArray <xarray.DataArray>`
       Hybrid B coefficients (sigma term) for the model levels.
   vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
       List of target pressure levels for interpolation (in specified units).
   vertical_input_dim : :py:class:`str <str>`
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`
       Units for the output pressure levels (must be convertible to hPa).
   variable_flag : ``Literal["T", "Z", "other"]``
       Indicates the type of variable being interpolated:
       - "T": Temperature (uses ECMWF extrapolation if enabled)
       - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
       - "other": Any other variable (uses lowest level value for extrapolation)
   temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Temperature at the lowest model level (required for 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Surface geopotential (required for 'T' or 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   interp_method : ``Literal["linear", "log", "loglog"]``, optional
       Interpolation method:

       - ``"linear"``: Linear interpolation
       - ``"log"``: Logarithmic interpolation
       - ``"loglog"``: Log-log interpolation

       Default is ``"linear"``.
   extrapolation : :py:class:`bool <bool>`, optional
       Whether to extrapolate below the lowest model level when needed.
       Default is False.
   lon_dim : :py:class:`str <str>`, optional
       Name of the longitude dimension. Default is "lon".
   lat_dim : :py:class:`str <str>`, optional
       Name of the latitude dimension. Default is "lat".
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is "time".
   p0_hPa : :py:class:`float <float>`, optional
       Reference pressure in **hPa** for hybrid level calculation. Default is ``1000.0 hPa``.

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to ``vertical_output_level``.

   Notes
   -----
   - The hybrid level pressure is calculated as: :math:`P = A \cdot p_0 + B \cdot \mathrm{psfc}`
   - Output pressure levels are converted to hPa internally for calculations.
   - Missing values are converted to NaN in the output.
   - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
   - For 'T' or 'Z', tbot and phis must be provided when ``extrapolation=True``.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vinth2p_ecmwf.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
   >>> interp_vinth2p_ecmwf(
   ...     data_input=temp_data,
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     hybrid_A_coefficients=hyam,
   ...     hybrid_B_coefficients=hybm,
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     variable_flag="T",
   ...     temperature_bottom_data=tbot_data,
   ...     surface_geopotential_data=phis_data,
   ...     interp_method="log",
   ...     extrapolation=True
   ... )


.. py:function:: interp_vintp2p_ecmwf(data_input: xarray.DataArray, pressure_data: xarray.DataArray, pressure_data_units: Literal['hPa', 'Pa', 'mbar'], surface_pressure_data: xarray.DataArray, surface_pressure_data_units: Literal['hPa', 'Pa', 'mbar'], vertical_output_level: list[int | float], vertical_input_dim: str, vertical_output_dim: str, vertical_output_dim_units: str, variable_flag: Literal['T', 'Z', 'other'], temperature_bottom_data: Optional[xarray.DataArray] = None, surface_geopotential_data: Optional[xarray.DataArray] = None, interp_method: Literal['linear', 'log', 'loglog'] = 'linear', extrapolation: bool = False, lon_dim: str = 'lon', lat_dim: str = 'lat', time_dim: str = 'time') -> xarray.DataArray

   Interpolates data at multidimensional pressure levels to constant pressure coordinates and uses an ECMWF formulation to extrapolate values below ground.

   This function performs vertical interpolation of atmospheric data from input pressure levels
   to specified output pressure levels using the NCL's `vintp2p_ecmwf` algorithm implemented
   in Fortran. It supports ECMWF-specific extrapolation for temperature ('T') and geopotential
   height ('Z') below the lowest pressure level.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray <xarray.DataArray>`
       Input 3D field (e.g., temperature or geopotential) on pressure levels with dimensions,
       e.g., (time, lev, lat, lon).
   pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       3D pressure field corresponding to the input data levels, with dimensions,
       e.g., (time, lev, lat, lon).
   pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the pressure data.
   surface_pressure_data : :py:class:`xarray.DataArray <xarray.DataArray>`
       Surface pressure field with dimensions, e.g., (time, lat, lon).
   surface_pressure_data_units : Literal["hPa", "Pa", "mbar"]
       Units of the surface pressure data.
   vertical_output_level : list[:py:class:`int <int>` | :py:class:`float <float>`]
       List of target pressure levels for interpolation (in specified units).
   vertical_input_dim : :py:class:`str <str>`
       Name of the vertical dimension in the input data.
   vertical_output_dim : :py:class:`str <str>`
       Name to use for the vertical dimension in the output data.
   vertical_output_dim_units : :py:class:`str <str>`
       Units for the output pressure levels (must be convertible to hPa).
   variable_flag : ``Literal["T", "Z", "other"]``
       Indicates the type of variable being interpolated:
       - "T": Temperature (uses ECMWF extrapolation if enabled)
       - "Z": Geopotential height (uses ECMWF extrapolation if enabled)
       - "other": Any other variable (uses lowest level value for extrapolation)
   temperature_bottom_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Temperature at the lowest model level (required for 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   surface_geopotential_data : Optional[:py:class:`xarray.DataArray <xarray.DataArray>`], optional
       Surface geopotential (required for 'T' or 'Z' extrapolation).
       Dimensions, e.g., (time, lat, lon). Default is None.
   interp_method : ``Literal["linear", "log", "loglog"]``, optional
       Interpolation method:

       - "linear": Linear interpolation
       - "log": Logarithmic interpolation
       - "loglog": Log-log interpolation

       Default is "linear".
   extrapolation : :py:class:`bool <bool>`, optional
       Whether to extrapolate below the lowest pressure level when needed.
       Default is False.
   lon_dim : :py:class:`str <str>`, optional
       Name of the longitude dimension. Default is "lon".
   lat_dim : :py:class:`str <str>`, optional
       Name of the latitude dimension. Default is "lat".
   time_dim: :py:class:`str <str>`, default: `time`.
       The time coordinate dimension name. Default is "time".

   Returns
   -------
   :py:class:`xarray.DataArray <xarray.DataArray>`
       Interpolated data on pressure levels with dimensions, e.g., (time, plev, lat, lon),
       where plev corresponds to vertical_output_level.

   Notes
   -----
   - The input pressure levels are provided directly via ``pressure_data``.
   - Output pressure levels and surface pressure are converted to hPa internally for calculations.
   - Missing values are converted to NaN in the output.
   - ECMWF extrapolation is applied only when ``extrapolation=True`` and variable_flag is 'T' or 'Z'.
   - For 'T' or 'Z', ``temperature_bottom_data`` and ``surface_geopotential_data`` must be provided when ``extrapolation=True``.

   .. seealso::

       - https://www.ncl.ucar.edu/Document/Functions/Built-in/vintp2p_ecmwf.shtml

   Examples
   --------
   >>> # Interpolate temperature to pressure levels with ECMWF extrapolation
   >>> interp_vintp2p_ecmwf(
   ...     data=temp_data,
   ...     pressure_data=pres_data,
   ...     pressure_data_units="Pa",
   ...     surface_pressure_data=psfc_data,
   ...     surface_pressure_data_units="Pa",
   ...     vertical_output_level=[1000, 850, 700, 500, 300],
   ...     vertical_input_dim="lev",
   ...     vertical_output_dim="plev",
   ...     vertical_output_dim_units="hPa",
   ...     variable_flag="T",
   ...     temperature_bottom_data=tbot_data,
   ...     surface_geopotential_data=phis_data,
   ...     interp_method="log",
   ...     extrapolation=True
   ... )


