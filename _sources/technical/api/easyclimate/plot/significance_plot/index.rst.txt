easyclimate.plot.significance_plot
==================================

.. py:module:: easyclimate.plot.significance_plot

.. autoapi-nested-parse::

   Mapping areas of significance



Functions
---------

.. autoapisummary::

   easyclimate.plot.significance_plot.draw_significant_area_contourf
   easyclimate.plot.significance_plot.get_significance_point
   easyclimate.plot.significance_plot.draw_significant_area_scatter


Module Contents
---------------

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


