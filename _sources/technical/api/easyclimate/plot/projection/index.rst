easyclimate.plot.projection
===========================

.. py:module:: easyclimate.plot.projection

.. autoapi-nested-parse::

   Graph processing related functions



Attributes
----------

.. autoapisummary::

   easyclimate.plot.projection.check_return


Functions
---------

.. autoapisummary::

   easyclimate.plot.projection.transfer_xarray_lon_from180TO360
   easyclimate.plot.projection.assert_compared_version
   easyclimate.plot.projection.draw_Circlemap_PolarStereo
   easyclimate.plot.projection.add_lon_cyclic


Module Contents
---------------

.. py:function:: transfer_xarray_lon_from180TO360(data_input: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon') -> xarray.DataArray | xarray.Dataset

   Longitude conversion -180-180 to 0-360.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
        The spatio-temporal data to be calculated.
   lon_dim: :py:class:`str <str>`, default: `lon`.
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
   - ver1: :py:class:`float <float>`, version number 1
   - ver2: :py:class:`float <float>`, version number 2

   Returns
   -------
   :py:class:`int<int>`.

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
   lon_dim: :py:class:`str<str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.

   .. seealso
       :py:func:`xarray.DataArray.pad <xarray:xarray.DataArray.pad>`, :py:func:`cartopy.util.add_cyclic_point <cartopy:cartopy.util.add_cyclic_point>`


