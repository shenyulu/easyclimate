easyclimate.interp.mesh2mesh
============================

.. py:module:: easyclimate.interp.mesh2mesh

.. autoapi-nested-parse::

   Regridding



Functions
---------

.. autoapisummary::

   easyclimate.interp.mesh2mesh.interp_mesh2mesh


Module Contents
---------------

.. py:function:: interp_mesh2mesh(data_input: xarray.DataArray | xarray.Dataset, target_grid: xarray.DataArray | xarray.Dataset, lon_dim: str = 'lon', lat_dim: str = 'lat', method: str = 'linear')

   Regridding regular or lat-lon grid data.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       The spatio-temporal data to be calculated.
   target_grid: :py:class:`xarray.DataArray<xarray.DataArray>` or :py:class:`xarray.Dataset<xarray.Dataset>`
       Target grid to be regridding.

       :py:class:`xarray.DataArray<xarray.DataArray>` version sample

       .. code:: python

           target_grid = xr.DataArray(
               dims=('lat', 'lon'),
               coords={'lat': np.arange(-89, 89, 3) + 1 / 1.0, 'lon': np.arange(-180, 180, 3) + 1 / 1.0}
           )

       :py:class:`xarray.Dataset<xarray.Dataset>` version sample

       .. code:: python

           target_grid = xr.Dataset()
           target_grid['lat'] = np.arange(-89, 89, 3) + 1 / 1.0
           target_grid['lon'] = np.arange(-180, 180, 3) + 1 / 1.0

   lon_dim: :py:class:`str <str>`, default: `lon`.
       Longitude coordinate dimension name. By default extracting is applied over the `lon` dimension.
   lat_dim: :py:class:`str <str>`, default: `lat`.
       Latitude coordinate dimension name. By default extracting is applied over the `lat` dimension.
   method: :py:class:`str <str>`, default: `linear`.
       The methods of regridding.

       - `linear`: linear, bilinear, or higher dimensional linear interpolation.
       - `nearest`: nearest-neighbor regridding.
       - `cubic`: cubic spline regridding.
       - `conservative`: conservative regridding.

   Reference
   --------------
   https://github.com/EXCITED-CO2/xarray-regrid


