easyclimate.plot
================

.. py:module:: easyclimate.plot


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/plot/axisticker/index
   /technical/api/easyclimate/plot/bar/index
   /technical/api/easyclimate/plot/curved_quiver_plot/index
   /technical/api/easyclimate/plot/line/index
   /technical/api/easyclimate/plot/modplot/index
   /technical/api/easyclimate/plot/projection/index
   /technical/api/easyclimate/plot/quick_draw/index
   /technical/api/easyclimate/plot/significance_plot/index
   /technical/api/easyclimate/plot/taylor_diagram/index


Functions
---------

.. autoapisummary::

   easyclimate.plot.draw_Circlemap_PolarStereo
   easyclimate.plot.add_lon_cyclic
   easyclimate.plot.add_lon_cyclic_lonarray
   easyclimate.plot.draw_significant_area_contourf
   easyclimate.plot.get_significance_point
   easyclimate.plot.draw_significant_area_scatter
   easyclimate.plot.calc_TaylorDiagrams_metadata
   easyclimate.plot.draw_TaylorDiagrams_base
   easyclimate.plot.draw_TaylorDiagrams_metadata
   easyclimate.plot.calc_TaylorDiagrams_values
   easyclimate.plot.set_lon_format_axis
   easyclimate.plot.set_lat_format_axis
   easyclimate.plot.set_p_format_axis
   easyclimate.plot.quick_draw_spatial_basemap
   easyclimate.plot.quick_draw_rectangular_box
   easyclimate.plot.curved_quiver
   easyclimate.plot.add_curved_quiverkey
   easyclimate.plot.bar_plot_with_threshold
   easyclimate.plot.line_plot_with_threshold


Package Contents
----------------

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

   .. seealso::

       :py:func:`geocat.viz.util.set_map_boundary <geocat.viz:geocat.viz.util.set_map_boundary>`, :py:class:`cartopy.mpl.gridliner.Gridliner <cartopy:cartopy.mpl.gridliner.Gridliner>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_formatting_coordinates.py
       ./dynamic_docs/plot_da_bbo.py


.. py:function:: add_lon_cyclic(data_input: xarray.DataArray, inter: float, lon_dim: str = 'lon')

   Add a cyclic point to an array and optionally a corresponding coordinate.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   inter: :py:class:`float<float>`
       Longitude interval (assuming longitude is arranged in a sequence of equal differences).
   lon_dim: :py:class:`str<str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   .. seealso
       :py:func:`xarray.DataArray.pad <xarray:xarray.DataArray.pad>`, :py:func:`cartopy.util.add_cyclic_point <cartopy:cartopy.util.add_cyclic_point>`


.. py:function:: add_lon_cyclic_lonarray(data_input: xarray.DataArray, lon_array: numpy.array, lon_dim: str = 'lon')

   Add a cyclic point to an array and optionally a corresponding coordinate.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   inter: :py:class:`float<float>`
       Longitude interval (assuming longitude is arranged in a sequence of equal differences).
   lon_dim: :py:class:`str<str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   .. seealso
       :py:func:`xarray.DataArray.pad <xarray:xarray.DataArray.pad>`, :py:func:`cartopy.util.add_cyclic_point <cartopy:cartopy.util.add_cyclic_point>`


.. py:function:: draw_significant_area_contourf(p_value: xarray.DataArray, thresh: float = 0.05, lon_dim: str = 'lon', lat_dim: str = 'lat', ax: matplotlib.axes.Axes = None, hatches: str = '...', hatch_colors: str = 'k', reverse_level_plot: bool = False, **kwargs) -> matplotlib.contour.QuadContourSet

   Draw significant area by :py:func:`matplotlib.axes.Axes.contourf<matplotlib.axes.Axes.contourf>`.

   Parameters
   ----------
   p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The p value data.
   thresh: :py:class:`float <float>`.
       The threshold value.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   ax : :py:class:`matplotlib.axes.Axes`, optional.
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   hatches: :py:class:`list[str]`, default: `...`
       A list of cross hatch patterns to use on the filled areas. If None, no hatching will be added to the contour. Hatching is supported in the PostScript, PDF, SVG and Agg backends only.
   hatch_colors: :py:class:`str <str>`, default: `k`.
       The colors of the hatches.

       .. warning::

           The parameter ``hatch_colors`` is not support to changed now.

   reverse_level_plot: :py:class:`bool<bool>`, default: `False`.
       Whether to reverse the drawing area.
   **kwargs, optional:
       Additional keyword arguments to :py:func:`xarray.plot.contourf<xarray.plot.contourf>`.

   Returns
   -------
   :py:class:`matplotlib.contour.QuadContourSet<matplotlib.contour.QuadContourSet>`.


.. py:function:: get_significance_point(p_value: xarray.DataArray, thresh: float = 0.05, lon_dim: str = 'lon', lat_dim: str = 'lat') -> pandas.DataFrame

   Obtain longitude and latitude array values that meet the conditions within the threshold from a two-dimensional array of p-values

   Parameters
   ----------
   p_value: :py:class:`xarray.DataArray<xarray.DataArray>`.
       The p value data.
   thresh: :py:class:`float <float>`.
       The threshold value.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.

   Returns
   -------
   :py:class:`pandas.DataFrame <pandas.DataFrame>`.


.. py:function:: draw_significant_area_scatter(significant_points_dataframe: pandas.DataFrame, lon_dim: str = 'lon', lat_dim: str = 'lat', ax: matplotlib.axes.Axes = None, **kwargs)

   Draw significant area by :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

   Parameters
   ----------
   significant_points_dataframe: :py:class:`pandas.DataFrame <pandas.DataFrame>`.
       The data contains the significant points, which is obtained by the :py:func:`easyclimate.plot.get_significance_point <easyclimate.plot.get_significance_point>`.
   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   ax : :py:class:`matplotlib.axes.Axes`, optional
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   **kwargs, optional:
       Additional keyword arguments to :py:func:`matplotlib.axes.Axes.scatter<matplotlib.axes.Axes.scatter>`.

       .. attention::
           You must specify `kwargs = {'transform': ccrs.PlateCarree()}` (`import cartopy.crs as ccrs`) in the cartopy `GeoAxes` or `GeoAxesSubplot`, otherwise projection errors may occur.


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


.. py:function:: set_lon_format_axis(ax: matplotlib.axes.Axes = None, axis: str = 'x', **kwargs)

   Setting the axes in longitude format.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   axis: {'x', 'y'}, default: 'x'
       The axis to which the parameters are applied.
   **kwargs
       Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_formatting_coordinates.py


.. py:function:: set_lat_format_axis(ax: matplotlib.axes.Axes = None, axis: str = 'y', **kwargs)

   Setting the axes in latitude format.

   Parameters
   ----------
   ax : :py:class:`matplotlib.axes.Axes`
       The axes to which the boundary will be applied.
   axis: {'x', 'y'}, default: 'y'
       The axis to which the parameters are applied.
   **kwargs
       Additional keyword arguments to wrapped :py:func:`matplotlib.axis.Axis.set_major_formatter <matplotlib:matplotlib.axis.Axis.set_major_formatter>`.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_formatting_coordinates.py


.. py:function:: set_p_format_axis(ax: matplotlib.axes.Axes = None, axis: str = 'y', axis_limits: tuple = (1000, 100), ticker_step: float = 100)

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

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_formatting_coordinates.py


.. py:function:: quick_draw_spatial_basemap(nrows: int = 1, ncols: int = 1, figsize=None, central_longitude: float = 0.0, draw_labels: str | bool | list | dict = ['bottom', 'left'], gridlines_color: str = 'grey', gridlines_alpha: float = 0.5, gridlines_linestyle: str = '--', coastlines_edgecolor: str = 'black', coastlines_kwargs: dict = {'lw': 0.5})

   Create geographical and spatial base map.

   Parameters
   ----------
   nrows, ncols :py:class:`int <int>`, default: 1
       Number of rows/columns of the subplot grid.
   figsize: (:py:class:`float <float>`, :py:class:`float <float>`)
       Width, height in inches.
   central_longitude: :py:class:`float <float>`, default: 0.
       The central longitude for `cartopy.crs.PlateCarree` projection.
   draw_labels: :py:class:`str <str>` | :py:class:`bool <bool>` | :py:class:`list <list>` | :py:class:`dict <dict>`, default: ["bottom", "left"].
       Toggle whether to draw labels. For finer control, attributes of Gridliner may be modified individually.

       - string: `"x"` or `"y"` to only draw labels of the respective coordinate in the CRS.
       - list: Can contain the side identifiers and/or coordinate types to select which ones to draw. For all labels one would use `["x", "y", "top", "bottom", "left", "right", "geo"]`.
       - dict: The keys are the side identifiers `("top", "bottom", "left", "right")` and the values are the coordinates `("x", "y")`; this way you can precisely decide what kind of label to draw and where. For x labels on the bottom and y labels on the right you could pass in `{"bottom": "x", "left": "y"}`.

       Note that, by default, x and y labels are not drawn on left/right and top/bottom edges respectively unless explicitly requested.
   gridlines_color: :py:class:`str <str>`, default: `grey`.
       The parameter `color` for `ax.gridlines`.
   gridlines_alpha: :py:class:`float <float>`, default: `0.5`.
       The parameter `alpha` for `ax.gridlines`.
   gridlines_linestyle: :py:class:`str <str>`, default: `"--"`.
       The parameter `linestyle` for `ax.gridlines`.
   coastlines_edgecolor: :py:class:`str <str>`, default: `"black"`.
       The parameter `color` for `ax.coastlines`.
   coastlines_kwargs: :py:class:`float <float>`, default: ``{"lw": 0.5}``.
       The kwargs for `ax.coastlines`.

   Returns
   -------
   - fig: :py:class:`Figure <matplotlib:matplotlib.figure.Figure>`
   - ax: :py:class:`Axes <matplotlib:matplotlib.axes.Axes>` or array of Axes: `ax` can be either a single :py:class:`Axes <matplotlib:matplotlib.axes.Axes>` object, or an array of Axes objects if more than one subplot was created. The dimensions of the resulting array can be controlled with the squeeze keyword.

   .. seealso::
       :py:func:`matplotlib.pyplot.subplots <matplotlib:matplotlib.pyplot.subplots>`
       :py:func:`cartopy.mpl.geoaxes.GeoAxes.gridlines <cartopy:cartopy.mpl.geoaxes.GeoAxes.gridlines>`
       :py:func:`cartopy.mpl.geoaxes.GeoAxes.coastlines <cartopy:cartopy.mpl.geoaxes.GeoAxes.coastlines>`


.. py:function:: quick_draw_rectangular_box(lon1: float, lon2: float, lat1: float, lat2: float, ax: matplotlib.axes.Axes = None, **patches_kwargs)

   Create geographical rectangular box.

   Parameters
   ----------
   lon1, lon2: :py:class:`float <float>`.
       Rectangular box longitude point. The applicable value should be between -180 :math:`^\circ` and 360 :math:`^\circ`.
       `lon1` and `lon2` must have a certain difference, should not be equal,
       do not strictly require the size relationship between `lon1` and `lon2`.
   lat1, lat2: :py:class:`float <float>`.
       Rectangular box latitude point. The applicable value should be between -90 :math:`^\circ` and 90 :math:`^\circ`.
       `lat1` and `lat2` must have a certain difference, should not be equal,
       do not strictly require the size relationship between `lat1` and `lat2`.
   ax : :py:class:`matplotlib.axes.Axes`, optional.
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   **patches_kwargs:
       Patch properties. see more in :py:class:`matplotlib.patches.Patch <matplotlib.patches.Patch>`

   .. seealso::
       :py:class:`matplotlib.patches.Rectangle <matplotlib.patches.Rectangle>`


.. py:function:: curved_quiver(ds: xarray.Dataset, x: collections.abc.Hashable, y: collections.abc.Hashable, u: collections.abc.Hashable, v: collections.abc.Hashable, ax: matplotlib.axes.Axes | None = None, density=1, linewidth=None, color=None, cmap=None, norm=None, arrowsize=1, arrowstyle='-|>', transform=None, zorder=None, start_points=None, integration_direction='both', grains=15, broken_streamlines=True) -> easyclimate.plot.modplot.CurvedQuiverplotSet

   Plot streamlines of a vector flow.

   .. warning::

       This function is experimental and the API is subject to change. Please use with caution.

   Parameters
   ----------
   ds : :py:class:`xarray.Dataset`.
       Wind dataset.
   x : Hashable or None, optional.
       Variable name for x-axis.
   y : Hashable or None, optional.
       Variable name for y-axis.
   u : Hashable or None, optional.
       Variable name for the u velocity (in `x` direction).
   v : Hashable or None, optional.
       Variable name for the v velocity (in `y` direction).
   ax : :py:class:`matplotlib.axes.Axes`, optional.
       Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
   density : float or (float, float)
       Controls the closeness of streamlines. When ``density = 1``, the domain
       is divided into a 30x30 grid. *density* linearly scales this grid.
       Each cell in the grid can have, at most, one traversing streamline.
       For different densities in each direction, use a tuple
       (density_x, density_y).
   linewidth : float or 2D array
       The width of the streamlines. With a 2D array the line width can be
       varied across the grid. The array must have the same shape as *u*
       and *v*.
   color : color or 2D array
       The streamline color. If given an array, its values are converted to
       colors using *cmap* and *norm*.  The array must have the same shape
       as *u* and *v*.
   cmap, norm
       Data normalization and colormapping parameters for *color*; only used
       if *color* is an array of floats. See `~.Axes.imshow` for a detailed
       description.
   arrowsize : float
       Scaling factor for the arrow size.
   arrowstyle : str
       Arrow style specification.
       See `~matplotlib.patches.FancyArrowPatch`.
   start_points : (N, 2) array
       Coordinates of starting points for the streamlines in data coordinates
       (the same coordinates as the *x* and *y* arrays).
   zorder : float
       The zorder of the streamlines and arrows.
       Artists with lower zorder values are drawn first.
   integration_direction : {'forward', 'backward', 'both'}, default: 'both'
       Integrate the streamline in forward, backward or both directions.
   broken_streamlines : boolean, default: True
       If False, forces streamlines to continue until they
       leave the plot domain.  If True, they may be terminated if they
       come too close to another streamline.

   Returns
   -------
   CurvedQuiverplotSet
       Container object with attributes

       - ``lines``: `.LineCollection` of streamlines

       - ``arrows``: `.PatchCollection` containing `.FancyArrowPatch`
         objects representing the arrows half-way along streamlines.

           This container will probably change in the future to allow changes
           to the colormap, alpha, etc. for both lines and arrows, but these
           changes should be backward compatible.

   .. seealso::
       - https://github.com/matplotlib/matplotlib/issues/20038
       - https://github.com/kieranmrhunt/curved-quivers
       - https://github.com/Deltares/dfm_tools/issues/483
       - https://github.com/NCAR/geocat-viz/issues/4
       - https://docs.xarray.dev/en/stable/generated/xarray.Dataset.plot.streamplot.html#xarray.Dataset.plot.streamplot


.. py:function:: add_curved_quiverkey(curved_quiver: easyclimate.plot.modplot.CurvedQuiverplotSet, X: float, Y: float, U: float, label: str, ax: matplotlib.axes.Axes = None, color: str = 'black', angle: float = 0.0, labelpos: Literal['N', 'S', 'E', 'W'] = 'N', labelsep: float = 0.02, labelcolor: str = None, fontproperties: matplotlib.font_manager.FontProperties = None, zorder: float = None)

   Add a key to a quiver plot.

   The positioning of the key depends on X, Y, coordinates, and labelpos.
   If labelpos is 'N' or 'S', X, Y give the position of the middle of the key arrow.
   If labelpos is 'E', X, Y positions the head, and if labelpos is 'W', X, Y positions the tail;
   in either of these two cases, X, Y is somewhere in the middle of the arrow+label key object.

   .. warning::

       This function is experimental and the API is subject to change. Please use with caution.

   Parameters
   ----------
   Q : :py:class:`easyclimate.modplot.CurvedQuiverplotSet`
       A `.CurvedQuiverplotSet` object as returned by a call to `curved_quiver()`.
   X, Y : float
       The location of the key.
   U : float
       The length of the key.
   label : str
       The key label (e.g., length and units of the key).
   ax : :py:class:`matplotlib.axes.Axes`, optional.
       Axes on which to plot.
   angle : float, default: `0.0`.
       The angle of the key arrow, in degrees anti-clockwise from the
       horizontal axis.
   labelpos : {'N', 'S', 'E', 'W'}, default: `N`.
       Position the label above, below, to the right, to the left of the
       arrow, respectively.
   labelsep : float, default: `0.02`.
       Distance in inches between the arrow and the label.
   labelcolor : str.
           Label color.
   fontproperties : dict, optional
       A dictionary with keyword arguments accepted by the
       `~matplotlib.font_manager.FontProperties` initializer:
       *family*, *style*, *variant*, *size*, *weight*.
   zorder : float
           The zorder of the key.


.. py:function:: bar_plot_with_threshold(da: xarray.DataArray, width=0.8, threshold: float = 0, pos_color: str = 'red', neg_color: str = 'blue', ax=None, **kwargs) -> matplotlib.container.BarContainer

   Plot a bar chart with time for a 1D :py:class:`xarray.DataArray <xarray.DataArray>` with bars colored based on a threshold value.

   Parameters:
   -----------
   da : :py:class:`xarray.DataArray <xarray.DataArray>`
       1-dimensional data array to plot
   width: :py:class:`float <float>` or array-like, default: 0.8
       The width(s) of the bars.

       .. note::

           If x has units (e.g., datetime), then the width is converted to a multiple of the width relative to the difference units of the x values (e.g., time difference).

   threshold : :py:class:`float <float>`, optional
       Threshold value for color separation (default: 0)
   pos_color : :py:class:`str <str>`, optional
       Color for bars ≥ threshold (default: 'red')
   neg_color : :py:class:`str <str>`, optional
       Color for bars < threshold (default: 'blue')
   ax : matplotlib axes, optional
       Axes object to plot on (uses current axes if None)
   **kwargs :
       Additional arguments passed to plt.bar

   Returns:
   --------
   matplotlib.container.BarContainer
       The bar plot object

   .. seealso::

       :py:func:`matplotlib.pyplot.bar <matplotlib.pyplot.bar>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_ao_index.py


.. py:function:: line_plot_with_threshold(da: xarray.DataArray, threshold: float = 0, pos_color: str = 'red', neg_color: str = 'blue', ax=None, line_plot: bool = True, fill_pos_plot: bool = True, fill_neg_plot: bool = True, line_kwargs=None, fill_kwargs=None) -> tuple

   Plot a line chart with proper shading at threshold crossings.

   Parameters:
   -----------
   da : :py:class:`xarray.DataArray <xarray.DataArray>`
       1-dimensional data array
   threshold : :py:class:`float <float>`, optional
       Color separation threshold (default: 0)
   pos_color : :py:class:`str <str>`, optional
       Color for values ≥ threshold (default: 'red')
   neg_color : :py:class:`str <str>`, optional
       Color for values < threshold (default: 'blue')
   ax : matplotlib axes, optional
       Axes to plot on (default: current axes)
   line_kwargs : :py:class:`dict <dict>`, optional
       Arguments for plt.plot
   fill_kwargs : :py:class:`dict <dict>`, optional
       Arguments for plt.fill_between

   Returns:
   --------
   tuple
       (line plot, fill objects)

   .. seealso::

       :py:func:`matplotlib.lines.Line2D <matplotlib.lines.Line2D>`

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_basic_statistical_analysis.py
       ./dynamic_docs/plot_corr_reg.py


