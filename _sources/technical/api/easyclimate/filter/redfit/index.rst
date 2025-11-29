easyclimate.filter.redfit
=========================

.. py:module:: easyclimate.filter.redfit

.. autoapi-nested-parse::

   Red-noise spectra estimating



Functions
---------

.. autoapisummary::

   easyclimate.filter.redfit.calc_redfit
   easyclimate.filter.redfit.calc_redfit_cross


Module Contents
---------------

.. py:function:: calc_redfit(data: xarray.DataArray, timearray: numpy.array = None, nsim: int = 1000, mctest: bool = False, rhopre: float = -99.0, ofac: float = 1.0, hifac: float = 1.0, n50: int = 1, iwin: Literal['rectangular', 'welch', 'hanning', 'triangular', 'blackmanharris'] = 'rectangular')

   Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

   Parameters
   ----------
   data: :py:class:`xarray.DataArray<xarray.DataArray>`
       Input time series data
   timearray: :py:class:`numpy.array<numpy.array>`
       Time series data array
   nsim: :py:class:`int<int>`
       Number of Monte-Carlo simulations (1000-2000 should be o.k. in most cases)
   mctest: :py:class:`bool<bool>`
       Toggle calculation of false-alarm levels based on Monte-Carlo simulation,
       if set to `True` : perform Monte-Carlo test,
       if set to `False` : skip Monte-Carlo test (default).
   rhopre: :py:class:`float<float>`
       Prescibed value for :math:`\rho`; unused if < 0 (default = -99.0)
   ofac: :py:class:`float<float>`
       Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0)
   hifac: :py:class:`float<float>`
       Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0)
   n50: :py:class:`int<int>`
       Number of WOSA segments (with 50 % overlap)
   iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
       Window-type identifier used to suppress sidelobes in spectral analysis:
       ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional)

   .. caution::
       Parameters `ofac`, `hifac`, `n50` and window type are identical to the SPECTRUM program
       (see Schulz and Stattegger, 1997 for further details).
       Except mctest, hifac and rhopre all parameters must be specified.

   Returns
   -------
   The red-noise spectra (:py:class:`xarray.Dataset<xarray.Dataset>`).

   .. seealso::
       - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
       - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_redfit.py


.. py:function:: calc_redfit_cross(data_x: xarray.DataArray, data_y: xarray.DataArray, timearray_x: numpy.array = None, timearray_y: numpy.array = None, x_sign: bool = False, y_sign: bool = False, nsim: int = 1000, mctest: bool = True, mctest_phi: bool = True, rhopre_1: float = -999.0, rhopre_2: float = -999.0, ofac: float = 1.0, hifac: float = 1.0, n50: int = 1, alpha: float = 0.05, iwin: Literal['rectangular', 'welch', 'hanning', 'triangular', 'blackmanharris'] = 'rectangular')

   Estimating red-noise spectra directly from unevenly spaced paleoclimatic time series.

   Parameters
   ----------
   data_x::py:class:`xarray.DataArray<xarray.DataArray>`
       First input time series data
   data_y: :py:class:`xarray.DataArray<xarray.DataArray>`
       Second input time series data
   timearray_x: :py:class:`numpy.array<numpy.array>`
       First time series data array
   timearray_y: :py:class:`numpy.array<numpy.array>`
       Second time series data array
   x_sign: :py:class:`bool<bool>`
       Change the sign of the first time series:
       if `True`: The sign of the data is changed
       if `False`: The sign of the data is not changed (default)
   y_sign: :py:class:`bool<bool>`
       Change the sign of the second time series:
       if `True`: The sign of the data is changed
       if `False`: The sign of the data is not changed (default)
   nsim: :py:class:`int<int>`
       Number of Monte Carlo simulations (1000-2000 is recommended)
   mctest: :py:class:`bool<bool>`
       Estimate the significance of auto and coherency spectrum with Monte Carlo simulations
       if `True`: perform Monte Carlo simulations
       if `False`: do not perform Monte Carlo simulations
   mctest_phi: :py:class:`bool<bool>`
       Estimate Monte Carlo confidence interval for the phase spectrum
       if `True`: perform Monte Carlo simulations (mctest needs to be true as well)
       if `False`: do not perform Monte Carlo simulations
   rhopre_1: :py:class:`float<float>`
       Prescribed value for :math:`\rho` for the first time series, not used if :math:`\rho < 0` (default = -999.0).
   rhopre_2: :py:class:`float<float>`
       Prescribed value for :math:`\rho` for the second time series, not used if :math:`\rho< 0` (default = -999.0).
   ofac: :py:class:`float<float>`
       Oversampling factor for Lomb-Scargle Fourier transform (typical values: 2.0-4.0).
   hifac: :py:class:`float<float>`
       Max. frequency to analyze is set to hifac * <fNyq> (default = 1.0).
   n50: :py:class:`int<int>`
       Number of WOSA segments (with 50 % overlap)
   alpha: :py:class:`float<float>`
       Significance level (Note: only 0.01, 0.05 [default], or 0.1 are allowed).
   iwin: {"rectangular", "welch", "hanning", "triangular", "blackmanharris"}
       Window-type identifier used to suppress sidelobes in spectral analysis:
       ({"rectangular", "welch", "hanning", "triangular", "blackmanharris"}, optional).

   .. caution::
       Parameters ofac, hifac, n50 and window type are identical to the SPECTRUM program
       (see Schulz and Stattegger, 1997 for further details).
       Except mctest, hifac, rhopre(1) and rhopre(2) all parameters must be specified.

   .. seealso::
       - Schulz, M., & Mudelsee, M. (2002). REDFIT: estimating red-noise spectra directly from unevenly spaced paleoclimatic time series [Software]. Computers & Geosciences, 28(3), 421-426. https://doi.org/10.1016/S0098-3004(01)00044-9
       - https://www.marum.de/Prof.-Dr.-michael-schulz/Michael-Schulz-Software.html


