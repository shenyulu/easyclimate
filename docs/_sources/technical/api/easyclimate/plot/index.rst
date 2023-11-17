:py:mod:`easyclimate.plot`
==========================

.. py:module:: easyclimate.plot


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   axisticker/index.rst
   projection/index.rst
   significance_plot/index.rst
   taylor_diagram/index.rst
   wind/index.rst


Package Contents
----------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.plot.transfer_xarray_lon_from180TO360
   easyclimate.plot.assert_compared_version
   easyclimate.plot.draw_Circlemap_PolarStereo
   easyclimate.plot.add_lon_cyclic
   easyclimate.plot.draw_significant_area_contourf
   easyclimate.plot.get_significance_point
   easyclimate.plot.draw_significant_area_scatter
   easyclimate.plot.calc_correlation_coefficient
   easyclimate.plot.calc_standard_deviation
   easyclimate.plot.calc_centeredRMS
   easyclimate.plot.calc_Taylor_skill_score
   easyclimate.plot.calc_TaylorDiagrams_values
   easyclimate.plot.calc_TaylorDiagrams_metadata
   easyclimate.plot.draw_TaylorDiagrams_base
   easyclimate.plot.draw_TaylorDiagrams_metadata
   easyclimate.plot.set_lon_format_axis
   easyclimate.plot.set_lat_format_axis
   easyclimate.plot.set_p_format_axis



Attributes
~~~~~~~~~~

.. autoapisummary::

   easyclimate.plot.check_return


.. py:function:: transfer_xarray_lon_from180TO360(data_input: xr.DataArray | xr.Dataset, lon_dim: str = 'lon') -> xr.DataArray | xr.Dataset

   Longitude conversion -180-180 to 0-360.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`.

   .. seealso::
       :py:func:`transfer_xarray_lon_from360TO180 <transfer_xarray_lon_from360TO180>`


.. py:function:: assert_compared_version(ver1: float, ver2: float) -> int

   Compare python library versions.

   .. attention::
       - Only for incoming version numbers without alphabetic characters.
       - Based on this method, the version number comparison should result in the following `"10.12.2.6.5">"10.12.2.6"`.

   Parameters
   ----------
   - ver1: Version number 1
   - ver2: Version number 2

   Returns
   -------
   :py:class:`int<python.int>`.

   .. note::
       If `ver1<ver2`, return `-1`; If `ver1=ver2`, return `0`; If `ver1>ver2`, return `1`.

   Examples
   --------

   .. code:: python

       >>> import easyclimate as ecl
       >>> result = ecl.assert_compared_version("10.12.2.6.5", "10.12.2.6")
       >>> print(result)
       1


.. py:data:: check_return

   

.. py:function:: draw_Circlemap_PolarStereo(*, lat_range: tuple | list, add_gridlines: bool = True, lon_step: float = None, lat_step: float = None, ax: matplotlib.axes.Axes = None, draw_labels: bool = True, set_map_boundary_kwargs: dict = {}, gridlines_kwargs: dict = {})

   Utility function to set the boundary of ax to a path that surrounds a
   given region specified by latitude and longitude coordinates. This boundary
   is drawn in the projection coordinates and therefore follows any curves
   created by the projection. As of now, this works consistently for the
   North/South Polar Stereographic Projections.

   Parameters
   ----------
   lat_range : :py:class:`tuple`, :py:class:`list`.
       The two-tuple containing the start and end of the desired range of
       latitudes. The first entry must be smaller than the second entry.
       Both entries must be between [-90 , 90].
   add_gridlines: :py:class:`bool`.
       whether or not add gridlines and tick labels to a map.
   lon_step: :py:class:`float`.
       The step of grid lines in longitude.
   lat_step: :py:class:`float`.
       The step of grid lines in latitude.
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   draw_labels: :py:class:`bool`.
       Whether to draw labels. Defaults to `True`.
   **set_map_boundary_kwargs: :py:class:`dict`.
       Additional keyword arguments to wrapped :py:func:`geocat.viz.util.set_map_boundary <geocat.viz:geocat.viz.util.set_map_boundary>`.
   **gridlines_kwargs: :py:class:`dict`.
       Additional keyword arguments to wrapped :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy:cartopy.mpl.gridliner.Gridliner>`.
   .. seealso
       :py:func:`geocat.viz.util.set_map_boundary <geocat.viz:geocat.viz.util.set_map_boundary>`, :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy:cartopy.mpl.gridliner.Gridliner>`.


.. py:function:: add_lon_cyclic(data_input: xarray.DataArray, inter: float, lon_dim: str = 'lon')

   Add a cyclic point to an array and optionally a corresponding coordinate.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   inter: :py:class:`float<float>`
       Longitude interval (assuming longitude is arranged in a sequence of equal differences).
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   .. seealso
       :py:func:`xarray.DataArray.pad <xarray:xarray.DataArray.pad>`, :py:func:`cartopy.util.add_cyclic_point <cartopy:cartopy.util.add_cyclic_point>`


.. py:function:: draw_significant_area_contourf(p_value: xarray.DataArray, thresh: float = 0.05, lon_dim: str = 'lon', lat_dim: str = 'lat', ax: matplotlib.axes.Axes = None, hatches: str = '...', hatch_colors: str = 'k', reverse_level_plot: bool = False, **kwargs) -> matplotlib.contour.QuadContourSet

   Draw significant area by :py:func:`matplotlib.axes.Axes.contourf<matplotlib.axes.Axes.contourf>`.

   Parameters
   ----------
   p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The p value data.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   ax : :py:class:`matplotlib.axes.Axes`, optional.
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   hatches: `list[str]`, default: `...`
       A list of cross hatch patterns to use on the filled areas. If None, no hatching will be added to the contour. Hatching is supported in the PostScript, PDF, SVG and Agg backends only.
   hatch_colors, default: `k`.
       The colors of the hatches.
   reverse_level_plot: :py:class:`bool<python.bool>`, default: `False`.
       Whether to reverse the drawing area.
   **kwargs, optional:
       Additional keyword arguments to :py:func:`xarray.plot.contourf<xarray.plot.contourf>`.

       .. attention::
           You must specify `transform = ccrs.PlateCarree()` (`import cartopy.crs as ccrs`) in the cartopy `GeoAxes` or `GeoAxesSubplot`, otherwise projection errors may occur.

   Returns
   -------
   :py:class:`matplotlib.contour.QuadContourSet<matplotlib.contour.QuadContourSet>`.


.. py:function:: get_significance_point(p_value: xarray.DataArray, thresh: float = 0.05, lon_dim: str = 'lon', lat_dim: str = 'lat') -> (numpy.array, numpy.array)

   Obtain longitude and latitude array values that meet the conditions within the threshold from a two-dimensional array of p-values

   Parameters
   ----------
   p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The p value data.
   thresh: :py:class:`float<python.float>`.
       The threshold value.
   lon_dim: :py:class:`str<python.str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str<python.str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.    

   Returns
   -------
   :py:class:`matplotlib.contour.QuadContourSet<matplotlib.contour.QuadContourSet>`.


.. py:function:: draw_significant_area_scatter(point_lon: numpy.array, point_lat: numpy.array, ax: matplotlib.axes.Axes = None, **kwargs)

   Draw significant area by :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

   Parameters
   ----------
   point_lon: :py:class:`numpy.array<numpy.array>`.
       The longitude of significant points.
   point_lat: :py:class:`numpy.array<numpy.array>`.
       The latitude of significant points.
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   **kwargs, optional:
       Additional keyword arguments to :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

       .. attention::
           You must specify `transform = ccrs.PlateCarree()` (`import cartopy.crs as ccrs`) in the cartopy `GeoAxes` or `GeoAxesSubplot`, otherwise projection errors may occur.


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


.. py:function:: set_lon_format_axis(ax: matplotlib.axes.Axes, axis: str = 'x', **kwargs)

   Setting the axes in longitude format.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   axis: {'x', 'y'}, default: 'x'
       The axis to which the parameters are applied.
   **kwargs
       Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.


.. py:function:: set_lat_format_axis(ax: matplotlib.axes.Axes, axis: str = 'y', **kwargs)

   Setting the axes in latitude format.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   axis: {'x', 'y'}, default: 'y'
       The axis to which the parameters are applied.
   **kwargs
       Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.


.. py:function:: set_p_format_axis(ax: matplotlib.axes.Axes, axis: str = 'y', axis_limits: tuple = (1000, 100), ticker_step: float = 100)

   Setting the axes in logarithmic vertical barometric pressure format.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   axis: {'x', 'y'}, default: 'y'
       The axis to which the parameters are applied.
   axis_limits: :py:class:`tuple`, default `(1000, 100)`.
       Assuming that the distribution of coordinates exhibits an isotropic series distribution, 
       this item sets the maximum value (near surface air pressure) and the minimum value (near overhead air pressure).
   ticker_step: :py:class:`float`, default `100`.
       Assuming an isotropic series of coordinate distributions, the term sets the tolerance.


