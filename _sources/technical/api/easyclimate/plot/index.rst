easyclimate.plot
================

.. py:module:: easyclimate.plot


Submodules
----------

.. toctree::
   :maxdepth: 1

   /technical/api/easyclimate/plot/axisticker/index
   /technical/api/easyclimate/plot/projection/index
   /technical/api/easyclimate/plot/quick_draw/index
   /technical/api/easyclimate/plot/significance_plot/index
   /technical/api/easyclimate/plot/taylor_diagram/index


Functions
---------

.. autoapisummary::

   easyclimate.plot.draw_Circlemap_PolarStereo
   easyclimate.plot.add_lon_cyclic
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
   Taylor, K. E. (2001), Summarizing multiple aspects of model performance in a single diagram, J. Geophys. Res., 106(D7), 7183-7192, doi:`10.1029/2000JD900719 <https://doi.org/10.1029/2000JD900719>`__.



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


.. py:function:: quick_draw_spatial_basemap(nrows: int = 1, ncols: int = 1, central_longitude: float = 0.0, draw_labels: str | bool | list | dict = ['bottom', 'left'], gridlines_color: str = 'grey', gridlines_alpha: float = 0.5, gridlines_linestyle: str = '--', coastlines_edgecolor: str = 'black', coastlines_linewidths: float = 0.5)

   Create geographical and spatial base map.

   Parameters
   ----------
   nrows, ncols :py:class:`int <int>`, default: 1
       Number of rows/columns of the subplot grid.
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
       The parameter `edgecolor` for `ax.coastlines`.
   coastlines_linewidths: :py:class:`float <float>`, default: `0.5`.
       The parameter `linewidths` for `ax.coastlines`.

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


