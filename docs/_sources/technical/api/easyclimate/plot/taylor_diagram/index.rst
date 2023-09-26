:py:mod:`easyclimate.plot.taylor_diagram`
=========================================

.. py:module:: easyclimate.plot.taylor_diagram

.. autoapi-nested-parse::

   Functions for Taylor Diagrams plotting.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.plot.taylor_diagram.calc_correlation_coefficient
   easyclimate.plot.taylor_diagram.calc_standard_deviation
   easyclimate.plot.taylor_diagram.calc_centeredRMS
   easyclimate.plot.taylor_diagram.calc_Taylor_skill_score
   easyclimate.plot.taylor_diagram.calc_TaylorDiagrams_values
   easyclimate.plot.taylor_diagram.calc_TaylorDiagrams_metadata
   easyclimate.plot.taylor_diagram.draw_TaylorDiagrams_base
   easyclimate.plot.taylor_diagram.draw_TaylorDiagrams_metadata



.. py:function:: calc_correlation_coefficient(f: xarray.DataArray, r: xarray.DataArray) -> xarray.DataArray

   Calculate the correlation coefficient.

   .. math::
       R = \frac{\frac{1}{N} \sum_{n=1}^N (f_n -\bar{f})(r_n - \bar{r})}{\sigma_f \sigma_r}

   Parameters
   ----------
   - f: (xarray DataArray, required)
       A spatial array of models to be compared.
   - r: (xarray DataArray, required)
       A spatial array for model reference comparisons (observations).

   .. attention::
       `f` and `r` must have the same two-dimensional spatial dimension.

   Returns
   -------
   - numpy array
       Pattern correlation coefficient array.



.. py:function:: calc_standard_deviation(f: xarray.DataArray, ddof=0) -> xarray.DataArray

   Calculate the standard deviation.

   .. math::
       STD = \frac{1}{N} \ \sum_{n=1}^{N} \left ( x_n - \bar{x} \ \right ) ^2

   Parameters
   ----------
   - f: (xarray DataArray, required)
       A spatial array of models to be compared.
   - ddof: (int, optional)
       "Delta Degrees of Freedom": the divisor used in the calculation is `N - ddof`, where `N` represents the number of elements.

   .. attention::
       `f` and `r` must have the same two-dimensional spatial dimension.

   Returns
   -------
   - numpy array
       Pattern standard deviation array.



.. py:function:: calc_centeredRMS(f: xarray.DataArray, r: xarray.DataArray) -> xarray.DataArray

   Calculate the center root mean square error.

   .. math::
       E'= \left \lbrace \frac{1}{N} \sum_{n=1}^N \left[ (f_n - \bar{f}) - (r_n - \bar{r}) \right] ^2  \right \rbrace ^{1/2}

   Parameters
   ----------
   - f: (xarray DataArray, required)
       A spatial array of models to be compared.
   - r: (xarray DataArray, required)
       A spatial array for model reference comparisons (observations).

   .. attention::
       `f` and `r` must have the same two-dimensional spatial dimension.
           
   Returns
   -------
   - centerRMS : numpy array
       Pattern center root mean square error array.



.. py:function:: calc_Taylor_skill_score(r: xarray.DataArray, sigma_f: xarray.DataArray, sigma_r: xarray.DataArray, r0=0.999) -> xarray.DataArray

   Calculate Taylor skill score (TSS).

   .. math::
       TSS = \frac{4 * (1+r)^4}{(SDR + \frac{1}{1+SDR})^2 (1+r_0)^4}

   Parameters
   ----------
   - r: (Float, required)
       The correlation coefficient value.
   - sigma_f: (Float, required)
       The standard deviation value of the model.
   - sigma_r: (Float, required)
       The standard deviation value of the observation.
   - r0: (Float, optional)
       Maximum correlation obtainable.


.. py:function:: calc_TaylorDiagrams_values(f: xarray.DataArray, r: xarray.DataArray, model_name: str, weighted=False, lat_dim='lat', normalized=True, r0=0.999) -> pandas.DataFrame

   Calculate the center root mean square error.

   where :math:`N` is the number of points in spatial pattern.

   Parameters
   ----------
   - f : (xarray DataArray, required)
       A spatial array of models to be compared.
   - r : (xarray DataArray, required)
       A spatial array for model reference comparisons (observations).
   - model_name: (str, required)
       The name of the model.
   - weighted: (bool, default `False`)
       Whether to weight the data by latitude or not? The default value is `False`.
   - lat_dim: (str, default `lat`)
       The name of `latitude` coordinate name.
   - normalized: (bool, default `True`, optional)
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
   - r0 : (Float, optional)
       Maximum correlation obtainable.

   .. attention::
       `f` and `r` must have the same two-dimensional spatial dimension.
           
   Returns
   -------
   pandas DataFrame.

   Reference
   --------------
   Taylor, K. E. (2001), Summarizing multiple aspects of model 
     performance in a single diagram, J. Geophys. Res., 106(D7),
     7183-7192, doi:`10.1029/2000JD900719 <https://doi.org/10.1029/2000JD900719>`__.



.. py:function:: calc_TaylorDiagrams_metadata(f: xarray.DataArray, r: xarray.DataArray, models_name=[], weighted=False, lat_dim='lat', normalized=True)

   Calculating Taylor diagram metadata

   Parameters
   ----------
   - f: (xarray DataArray, required)
       A spatial array of models to be compared.
   - r: (xarray DataArray, required)
       A spatial array for model reference comparisons (observations).
   - weighted: (bool, default `False`)
       Whether to weight the data by latitude or not? The default value is `False`.
   - lat_dim: (str, default `lat`)
       The name of `latitude` coordinate name.
   - model_name: (list str, required)
       The `list` of the models' name.
   - normalized: (bool, default `True`, optional)
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.

   Returns
   --------------
   pandas.DataFrame

   Examples
   ---------------

   .. code:: python

       >>> import xarray as xr
       >>> import pandas as pd
       >>> import numpy as np
       >>> import easyclimate as ecl
       >>> da_a = xr.DataArray(
       ...:     np.array([[1, 2, 3], [0.1, 0.2, 0.3], [3.2, 0.6, 1.8]]),
       ...:     dims = ("lat", "time"),
       ...:     coords= {'lat': np.array([-30, 0, 30]),
       ...:              'time': pd.date_range("2000-01-01", freq="D", periods=3)
       ...:              }
       ...:)
       >>> da_a
       <xarray.DataArray (lat: 3, time: 3)>
       array([[1. , 2. , 3. ],
           [0.1, 0.2, 0.3],
           [3.2, 0.6, 1.8]])
       Coordinates:
       * lat      (lat) int32 -30 0 30
       * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
       >>>  da_b = xr.DataArray(
       ...:     np.array([[0.2, 0.4, 0.6], [15, 10, 5], [3.2, 0.6, 1.8]]),
       ...:     dims = ("lat", "time"),
       ...:     coords= {'lat': np.array([-30, 0, 30]),
       ...:              'time': pd.date_range("2000-01-01", freq="D", periods=3)
       ...:              }
       ...:)
       >>>  da_b
       <xarray.DataArray (lat: 3, time: 3)>
       array([[ 0.2,  0.4,  0.6],
           [15. , 10. ,  5. ],
           [ 3.2,  0.6,  1.8]])
       Coordinates:
       * lat      (lat) int32 -30 0 30
       * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
       >>>  da_obs = (da_a + da_b) / 1.85
       >>>  da_obs
       <xarray.DataArray (lat: 3, time: 3)>
       array([[0.64864865, 1.2972973 , 1.94594595],
           [8.16216216, 5.51351351, 2.86486486],
           [3.45945946, 0.64864865, 1.94594595]])
       Coordinates:
       * lat      (lat) int32 -30 0 30
       * time     (time) datetime64[ns] 2000-01-01 2000-01-02 2000-01-03
       >>>  ecl.calc_TaylorDiagrams_metadata(
       ...:     f = [da_a, da_b],
       ...:     r = [da_obs, da_obs],
       ...:     models_name = ['f1', 'f2'],
       ...:     weighted = True,
       ...:     normalized = True,
       ...:)
       item       std                   cc  centeredRMS       TSS
       0  Obs  1.000000                  1.0     0.000000  1.002003
       1   f1  0.404621  -0.4293981636461462     1.229311  0.003210
       2   f2  2.056470    0.984086060161888     1.087006  0.600409



.. py:function:: draw_TaylorDiagrams_base(ax=None, half_circle=False, normalized=True, std_min=0, std_max=2, std_interval=0.25, arc_label='Correlation', arc_label_pad=0.2, arc_label_kwargs={'fontsize': 12}, arc_ticker_kwargs={'lw': 0.8, 'c': 'black'}, arc_tickerlabel_kwargs={'labelsize': 12, 'pad': 8}, arc_ticker_length=0.02, arc_minorticker_length=0.01, x_label='Std (Normalized)', x_label_pad=0.25, x_label_kwargs={'fontsize': 12}, x_ticker_length=0.02, x_tickerlabel_kwargs={'fontsize': 12}, x_ticker_kwargs={'lw': 0.8, 'c': 'black'}, y_ticker_kwargs={'lw': 0.8, 'c': 'black'}) -> matplotlib.collections.Collection

   Drawing Taylor Graphics Basic Framework

   Parameters
   ----------
   - ax: (matplotlib axes object, optional)
       Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
   - half_circle: (bool, default `False`, optional)
       Whether to draw the `'half-circle'` version of the Taylor diagram.
   - normalized: (bool, default `True`, optional)
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations. 
       This parameter mainly affects the label `x=1` on the `x` axis, if normalized to True, it is rewritten to `REF`.
   - std_min: (float, default `0.0`, optional)
       Minimum value of x-axis (standard deviation) on Taylor diagram.

       .. note:: The value of `std_min` shoud be 0 in the `'half-circle'` version of the Taylor diagram.

   - std_max: (float, default `2.0`, optional)
       Maximum value of x-axis (standard deviation) on Taylor diagram.
   - std_interval: (float, default `0.25`, optional)
       The interval between the ticker on the x-axis (standard deviation) between the minimum and maximum values on the Taylor diagram.
   - arc_label: (str, default `'Correlation'`, optional)
       Label on Taylor chart arc, default value is `'Correlation'`.
   - arc_label_pad: (float, default `0.2`, optional)
       The offset of the title from the top of the arc, based on x-axis based coordinate system.
   - arc_label_kwargs: (dict, default `{'fontsize': 12}`, optional)
       Additional keyword arguments passed on to labels on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   - arc_ticker_kwargs: (dict, default `{'lw': 0.8, 'c': 'black'}`, optional)
       Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   - arc_tickerlabel_kwargs: (dict, default `{'labelsize': 12, 'pad': 8}`, optional)
       Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.tick_params`.
   - arc_ticker_length: (float, default `0.02`, optional)
       Ticker length on arc.
   - arc_minorticker_length: (float, default `0.01`, optional)
       Minor ticker length on arc.
   - x_label: (str, default `'Std (Normalized)'`, optional)
       Label on Taylor chart x axis, default value is `'Std (Normalized)'`.
   - x_label_pad: (float, default `0.25`, optional)
       The offset of the title from the top of the x-axis, based on x-axis based coordinate system.
   - x_label_kwargs: (dict, default `{'fontsize': 12}`, optional)
       Additional keyword arguments passed on to labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   - x_ticker_length: (float, default `0.02`, optional)
       Ticker length on x-axis
   - x_tickerlabel_kwargs: (dict, default `{'fontsize': 12}`, optional)
       Additional keyword arguments passed on to tickers' labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   - x_ticker_kwargs: (dict, default `{'lw': 0.8, 'c': 'black'}`, optional)
       Additional keyword arguments passed on to tickers on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   - y_ticker_kwargs: (dict, default `{'lw': 0.8, 'c': 'black'}`, optional)
       Additional keyword arguments passed on to tickers on y-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.

   Returns
   -------
   `matplotlib.collections.Collection`


.. py:function:: draw_TaylorDiagrams_metadata(taylordiagrams_metadata, marker_list, color_list, label_list, legend_list, ax=None, normalized=True, cc='cc', std='std', point_label_xoffset=0, point_label_yoffset=0.05, point_kwargs={'alpha': 1, 'markersize': 6.5}, point_label_kwargs={'fontsize': 14}) -> matplotlib.collections.Collection

   Draw points to Taylor Graphics Basic Framework according to Taylor diagram metadata.

   Parameters
   ----------
   - taylordiagrams_metadata: (pandas.DataFrame, required)
       Taylor diagram metadata generated by the function `calc_TaylorDiagrams_metadata`.
   - marker_list: (list, required)
       The list of markers. The order of `marker` in `marker_list` is determined by the order in `taylordiagrams_metadata`. 
       See `matplotlib.markers` for full description of possible arguments.
   - color_list: (list, required)
       The list of colors. The order of `color` in `color_list` is determined by the order in `taylordiagrams_metadata`. 
   - label_list: (list, required)
       The list of data point labels (marked next to plotted points). 
       The order of label in `label_list` is determined by the order in `taylordiagrams_metadata`.
   - legend_list: (list, required)
       The list of legend label. 
       The order of label in `legend_list` is determined by the order in `taylordiagrams_metadata`.
   - ax: (matplotlib axes object, optional)
       Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
   - normalized: (bool, default `True`, optional)
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
   - cc: (str, default `'cc'`, optional)
       The name of correlation coefficient in `taylordiagrams_metadata`.
   - std: (str, default `'std'`, optional)
       The name of standard deviation in `taylordiagrams_metadata`.
   - point_label_xoffset: (float, optional)
       The offset of the labels from the points, based on x-axis based coordinate system.
   - point_label_yoffset: (float, optional)
       The offset of the labels from the points, based on y-axis based coordinate system.
   - point_kwargs: (dict, optional)
       Additional keyword arguments passed on to data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   - point_label_kwargs: (dict, optional)
       Additional keyword arguments passed on to the labels of data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.

   Returns
   -------
   `matplotlib.collections.Collection`


