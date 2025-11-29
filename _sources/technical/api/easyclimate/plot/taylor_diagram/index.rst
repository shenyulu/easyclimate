easyclimate.plot.taylor_diagram
===============================

.. py:module:: easyclimate.plot.taylor_diagram

.. autoapi-nested-parse::

   Functions for Taylor Diagrams plotting.



Functions
---------

.. autoapisummary::

   easyclimate.plot.taylor_diagram.calc_TaylorDiagrams_values
   easyclimate.plot.taylor_diagram.calc_TaylorDiagrams_metadata
   easyclimate.plot.taylor_diagram.draw_TaylorDiagrams_base
   easyclimate.plot.taylor_diagram.draw_TaylorDiagrams_metadata


Module Contents
---------------

.. py:function:: calc_TaylorDiagrams_values(f: xarray.DataArray, r: xarray.DataArray, model_name: str, weighted: bool = False, lat_dim: str = 'lat', normalized: bool = True, r0: float = 0.999) -> pandas.DataFrame

   Calculate the center root mean square error.

   where :math:`N` is the number of points in spatial pattern.

   Parameters
   ----------
   f : :py:class:`xarray.DataArray<xarray.DataArray>`, required.
       A spatial array of models to be compared.
   r : :py:class:`xarray.DataArray<xarray.DataArray>`, required.
       A spatial array for model reference comparisons (observations).
   model_name: :py:class:`str <str>`, required.
       The name of the model.
   weighted: :py:class:`bool <bool>`, default `False`.
       Whether to weight the data by latitude or not? The default value is `False`.
   lat_dim: :py:class:`str <str>`, default `lat`.
       The name of `latitude` coordinate name.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
   r0 : :py:class:`float <float>`, optional.
       Maximum correlation obtainable.

   .. attention::
       `f` and `r` must have the same two-dimensional spatial dimension.

   Returns
   -------
   :py:class:`pandas.DataFrame <pandas.DataFrame>`.

   Reference
   --------------
   - Taylor, K. E. (2001), Summarizing multiple aspects of model performance in a single diagram, J. Geophys. Res., 106(D7), 7183-7192, doi:`10.1029/2000JD900719 <https://doi.org/10.1029/2000JD900719>`__.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_taylor_diagram.py


.. py:function:: calc_TaylorDiagrams_metadata(f: xarray.DataArray, r: xarray.DataArray, models_name: list[str] = [], weighted: bool = False, lat_dim: str = 'lat', normalized: bool = True)

   Calculating Taylor diagram metadata

   Parameters
   ----------
   f: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
       A spatial array of models to be compared.
   r: :py:class:`xarray.DataArray<xarray.DataArray>`, required.
       A spatial array for model reference comparisons (observations).
   model_name: :py:class:`list[str]`, required.
       The `list` of the models' name.
   weighted: :py:class:`bool <bool>`, default `False`.
       Whether to weight the data by latitude or not? The default value is `False`.
   lat_dim: :py:class:`str <str>`, default `lat`.
       The name of `latitude` coordinate name.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.

   Returns
   --------------
   :py:class:`pandas.DataFrame <pandas.DataFrame>`.

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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_taylor_diagram.py


.. py:function:: draw_TaylorDiagrams_base(ax: matplotlib.axes.Axes = None, half_circle: bool = False, normalized: bool = True, std_min: float = 0, std_max: float = 2, std_interval: float = 0.25, arc_label: str = 'Correlation', arc_label_pad: float = 0.2, arc_label_kwargs: dict = {'fontsize': 12}, arc_ticker_kwargs: dict = {'lw': 0.8, 'c': 'black'}, arc_tickerlabel_kwargs: dict = {'labelsize': 12, 'pad': 8}, arc_ticker_length: float = 0.02, arc_minorticker_length: float = 0.01, x_label: str = 'Std (Normalized)', x_label_pad: float = 0.25, x_label_kwargs: dict = {'fontsize': 12}, x_ticker_length: float = 0.02, x_tickerlabel_kwargs: dict = {'fontsize': 12}, x_ticker_kwargs: dict = {'lw': 0.8, 'c': 'black'}, y_ticker_kwargs: dict = {'lw': 0.8, 'c': 'black'}) -> matplotlib.collections.Collection

   Drawing Taylor Graphics Basic Framework

   Parameters
   ----------
   ax: :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional.
       Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
   half_circle: :py:class:`bool <bool>`, default `False`, optional.
       Whether to draw the `'half-circle'` version of the Taylor diagram.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
       This parameter mainly affects the label `x=1` on the `x` axis, if normalized to True, it is rewritten to `REF`.
   std_min: :py:class:`float <float>`, default `0.0`, optional.
       Minimum value of x-axis (standard deviation) on Taylor diagram.

       .. note:: The value of `std_min` shoud be 0 in the `'half-circle'` version of the Taylor diagram.

   std_max: :py:class:`float <float>`, default `2.0`, optional.
       Maximum value of x-axis (standard deviation) on Taylor diagram.
   std_interval: :py:class:`float <float>`, default `0.25`, optional.
       The interval between the ticker on the x-axis (standard deviation) between the minimum and maximum values on the Taylor diagram.
   arc_label: :py:class:`str <str>`, default `'Correlation'`, optional.
       Label on Taylor chart arc, default value is `'Correlation'`.
   arc_label_pad: :py:class:`float <float>`, default `0.2`, optional.
       The offset of the title from the top of the arc, based on x-axis based coordinate system.
   arc_label_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
       Additional keyword arguments passed on to labels on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   arc_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
       Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   arc_tickerlabel_kwargs: :py:class:`dict <dict>`, default `{'labelsize': 12, 'pad': 8}`, optional.
       Additional keyword arguments passed on to tickers on arcs, according to other miscellaneous parameters in`matplotlib.axes.Axes.tick_params`.
   arc_ticker_length: :py:class:`float <float>`, default `0.02`, optional.
       Ticker length on arc.
   arc_minorticker_length: :py:class:`float <float>`, default `0.01`, optional.
       Minor ticker length on arc.
   x_label: :py:class:`str <str>`, default `'Std (Normalized)'`, optional.
       Label on Taylor chart x axis, default value is `'Std (Normalized)'`.
   x_label_pad: :py:class:`float <float>`, default `0.25`, optional.
       The offset of the title from the top of the x-axis, based on x-axis based coordinate system.
   x_label_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
       Additional keyword arguments passed on to labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   x_ticker_length: :py:class:`float <float>`, default `0.02`, optional.
       Ticker length on x-axis
   x_tickerlabel_kwargs: :py:class:`dict <dict>`, default `{'fontsize': 12}`, optional.
       Additional keyword arguments passed on to tickers' labels on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.
   x_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
       Additional keyword arguments passed on to tickers on x-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   y_ticker_kwargs: :py:class:`dict <dict>`, default `{'lw': 0.8, 'c': 'black'}`, optional.
       Additional keyword arguments passed on to tickers on y-axis, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.

   Returns
   -------
   :py:class:`matplotlib.collections.Collection <matplotlib.collections.Collection>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_taylor_diagram.py


.. py:function:: draw_TaylorDiagrams_metadata(taylordiagrams_metadata: pandas.DataFrame, marker_list: list, color_list: list, label_list: list, legend_list: list, ax: matplotlib.axes.Axes = None, normalized: bool = True, cc: str = 'cc', std: str = 'std', point_label_xoffset: float = 0, point_label_yoffset: float = 0.05, point_kwargs: dict = {'alpha': 1, 'markersize': 6.5}, point_label_kwargs: dict = {'fontsize': 14}) -> matplotlib.collections.Collection

   Draw points to Taylor Graphics Basic Framework according to Taylor diagram metadata.

   Parameters
   ----------
   taylordiagrams_metadata: :py:class:`pandas.DataFrame <pandas.DataFrame>`, required.
       Taylor diagram metadata generated by the function `calc_TaylorDiagrams_metadata`.
   marker_list: :py:class:`list <list>`, required.
       The list of markers. The order of `marker` in `marker_list` is determined by the order in `taylordiagrams_metadata`.
       See `matplotlib.markers` for full description of possible arguments.
   color_list: :py:class:`list <list>`, required.
       The list of colors. The order of `color` in `color_list` is determined by the order in `taylordiagrams_metadata`.
   label_list: :py:class:`list <list>`, required.
       The list of data point labels (marked next to plotted points).
       The order of label in `label_list` is determined by the order in `taylordiagrams_metadata`.
   legend_list: :py:class:`list <list>`, required.
       The list of legend label.
       The order of label in `legend_list` is determined by the order in `taylordiagrams_metadata`.
   ax: :py:class:`matplotlib.axes.Axes <matplotlib.axes.Axes>`, optional.
       Axes on which to plot. By default, use the current axes, i.e. `ax = plt.gca()`.
   normalized: :py:class:`bool <bool>`, default `True`, optional.
       Whether the standard deviations is normalized, that is standard deviation of the model divided by that of the observations.
   cc: :py:class:`str <str>`, default `'cc'`, optional.
       The name of correlation coefficient in `taylordiagrams_metadata`.
   std: :py:class:`str <str>`, default `'std'`, optional.
       The name of standard deviation in `taylordiagrams_metadata`.
   point_label_xoffset: :py:class:`float <float>`, optional.
       The offset of the labels from the points, based on x-axis based coordinate system.
   point_label_yoffset: :py:class:`float <float>`, optional.
       The offset of the labels from the points, based on y-axis based coordinate system.
   point_kwargs: :py:class:`dict <dict>`, optional.
       Additional keyword arguments passed on to data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.plot`.
   point_label_kwargs: :py:class:`dict <dict>`, optional.
       Additional keyword arguments passed on to the labels of data points, according to other miscellaneous parameters in`matplotlib.axes.Axes.text`.

   Returns
   -------
   :py:class:`matplotlib.collections.Collection <matplotlib.collections.Collection>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_taylor_diagram.py


