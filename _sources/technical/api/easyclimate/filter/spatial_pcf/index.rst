easyclimate.filter.spatial_pcf
==============================

.. py:module:: easyclimate.filter.spatial_pcf

.. autoapi-nested-parse::

   2D spatial parabolic cylinder function



Functions
---------

.. autoapisummary::

   easyclimate.filter.spatial_pcf.filter_2D_spatial_parabolic_cylinder_function


Module Contents
---------------

.. py:function:: filter_2D_spatial_parabolic_cylinder_function(zonal_wind_speed_data: xarray.DataArray, meridional_wind_speed_data: xarray.DataArray, z_data: xarray.DataArray, period_min=3.0, period_max=30.0, wavenumber_min=2, wavenumber_max=20, trapping_scale_deg=6.0, lon_dim='lon', lat_dim='lat', time_dim: str = 'time', complex_dtype=np.complex128, real_dtype=np.float64)

   Perform space-time spectral analysis to extract equatorial wave components.

   This function filters atmospheric fields to isolate equatorial wave modes, including Kelvin waves,
   westward-moving mixed Rossby-gravity (WMRG) waves, and equatorial Rossby waves of the first and
   second kind (R1 and R2), within user-specified period and wavenumber ranges. It processes input
   data using a combination of detrending, windowing, Fourier transforms, and projection onto
   parabolic cylinder functions to decompose the fields into these wave types.

   Parameters
   ----------

   zonal_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
       The zonal (east-west) wind speed component.
   meridional_wind_speed_data: :class:`xarray.DataArray <xarray.DataArray>`
       The meridional (north-south) wind speed component.
   z_data: :class:`xarray.DataArray <xarray.DataArray>`
       The geopotential height.
   period_min: :class:`float`, optional (default=3.0)
       The minimum period (in days) of the waves to be extracted.
   period_max: :class:`float`, optional (default=30.0)
       The maximum period (in days) of the waves to be extracted.
   wavenumber_min: :class:`int`, optional (default=2)
       The minimum zonal wavenumber of the waves to be extracted.
   wavenumber_max: :class:`int`, optional (default=20)
       The maximum zonal wavenumber of the waves to be extracted.
   trapping_scale_deg: :class:`float`, optional (default=6.0)
       The meridional trapping scale in degrees, defining the width of the equatorial waveguide
       for the parabolic cylinder functions.
   lon_dim: :class:`str`, optional (default="lon")
       The name of the longitude dimension in the input data arrays.
   lat_dim: :class:`str`, optional (default="lat")
       The name of the latitude dimension in the input data arrays.
   time_dim: :class:`str`, optional (default="time")
       The name of the time dimension in the input data arrays.

   Returns
   -------

   :class:`xarray.Dataset <xarray.Dataset>`

   A dataset containing the filtered fields for each wave type:
       - ``'u'``: Zonal wind component (m/s)
       - ``'v'``: Meridional wind component (m/s)
       - ``'z'``: Geopotential height (m)
   Each variable includes an additional ``'wave_type'`` dimension with values:
   ``['kelvin', 'wmrg', 'r1', 'r2']``.

   Mathematical Explanation
   ------------------------

   The function implements a space-time spectral analysis method to decompose atmospheric fields into
   equatorial wave modes, based on the normal modes of the shallow water equations on an equatorial
   beta-plane. The meridional structure of these waves is represented by parabolic cylinder functions,
   which arise as solutions to the quantum harmonic oscillator problem adapted to the equatorial
   waveguide.

   Key Steps
   ~~~~~~~~~

   1. Detrending and Windowing:
       - The input fields are spatially detrended to remove large-scale trends.
       - A Tukey window (alpha=0.1) is applied along the time dimension to minimize spectral leakage.

   2. Variable Transformation:
       Following Yang et al. (2003), new variables :math:`q` and :math:`r` are defined:

       .. math::
           q = \frac{g}{c_e} z + u, \quad r = \frac{g}{c_e} z - u

       where:
           - :math:`g = 9.8 \text{m/s}^2` is gravitational acceleration.
           - :math:`c_e = 2 \beta L^2` is a characteristic speed.
           - :math:`\beta = 2.3 \times 10^{-11}, \text{m}^{-1}\text{s}^{-1}` is the beta parameter.
           - :math:`L = \frac{2\pi R_e \cdot \text{trapping_scale_deg}}{360}` is the trapping scale in meters (:math:`R_e = 6.371 \times 10^6, \text{m}` is Earth's radius).
           - :math:`z` is geopotential height, and :math:`u` is zonal wind speed.

   3. Fourier Transform:
       A 2D Fast Fourier Transform (FFT) is applied along the time and longitude dimensions to convert
       the fields into frequency-wavenumber space.

   4. Projection onto Parabolic Cylinder Functions:
       - The spectral coefficients are projected onto parabolic cylinder functions :math:`D_n(y)`, where
           :math:`y = \text{latitude} / \text{trapping_scale_deg}` is the scaled latitude.
       - The first four modes (:math:`n = 0, 1, 2, 3`) are computed with normalization factors:

       .. math::

           D_0(y) = e^{-y^2/4}, \quad D_1(y) = y e^{-y^2/4}, \quad D_2(y) = (y^2 - 1) e^{-y^2/4}, \quad D_3(y) = y (y^2 - 3) e^{-y^2/4}

       - These functions define the meridional structure of the waves.

   5. Wave Type Identification and Filtering:
       - Kelvin Wave: Uses :math:`n = 0` mode, eastward propagating within specified frequency (:math:`1/\text{period_max}` to :math:`1/\text{period_min}`) and wavenumber ranges.
       - WMRG Wave: Combines :math:`n = 1` for :math:`q` and :math:`n = 0` for meridional wind :math:`v`, westward propagating.
       - R1 Wave: Uses :math:`n = 2` for :math:`q`, :math:`n = 0` for :math:`r`, and :math:`n = 1` for :math:`v`, westward propagating.
       - R2 Wave: Uses :math:`n = 3` for :math:`q`, :math:`n = 1` for :math:`r`, and :math:`n = 2` for :math:`v`, westward propagating.
       - Frequencies and wavenumbers are selected based on ``period_min``, ``period_max``, ``wavenumber_min``, and ``wavenumber_max``.

   6. Reconstruction:
       - The filtered spectral coefficients are recombined with the parabolic cylinder functions and transformed back to physical space using an inverse 2D FFT.
       - The physical fields are reconstructed as:
           - :math:`u = (q + r) / 2`.
           - :math:`v` directly from its coefficients.
           - :math:`z = (q - r) \cdot c_e / (2g)`.

   References
   ----------
   - Gill, A.E. (1980), Some simple solutions for heat-induced tropical circulation. Q.J.R. Meteorol. Soc., 106: 447-462. https://doi.org/10.1002/qj.49710644905
   - Li, X.f., Cho, HR. Development and propagation of equatorial waves. Adv. Atmos. Sci. 14, 323–338 (1997). https://doi.org/10.1007/s00376-997-0053-6
   - Yang, G.-Y., Hoskins, B., & Slingo, J. (2003). Convectively coupled equatorial waves: A new methodology for identifying wave structures in observational data. *Journal of the Atmospheric Sciences*, 60(14), 1637-1654.
   - Knippertz, P., Gehne, M., Kiladis, G.N., Kikuchi, K., Rasheeda Satheesh, A., Roundy, P.E., et al. (2022) The intricacies of identifying equatorial waves. Quarterly Journal of the Royal Meteorological Society, 148(747), 2814–2852. Available from: https://doi.org/10.1002/qj.4338

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_spatial_pcf.py


