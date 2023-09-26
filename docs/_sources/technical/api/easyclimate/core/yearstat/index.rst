:py:mod:`easyclimate.core.yearstat`
===================================

.. py:module:: easyclimate.core.yearstat

.. autoapi-nested-parse::

   This module computes statistical values over timesteps of the same year



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   easyclimate.core.yearstat.calc_yearmean
   easyclimate.core.yearstat.calc_yearsum
   easyclimate.core.yearstat.calc_yearstd
   easyclimate.core.yearstat.calc_yearvar
   easyclimate.core.yearstat.calc_yearmax
   easyclimate.core.yearstat.calc_yearmin



.. py:function:: calc_yearmean(data_input, dim='time', **kwargs)

   Calculate yearly mean.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{mean} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating mean on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.mean <numpy:numpy.mean>`, :py:func:`dask.array.mean <dask:dask.array.mean>`, 
       :py:meth:`xarray.DataArray.mean <xarray:xarray.DataArray.mean>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.mean <xarray:xarray.core.groupby.DataArrayGroupBy.mean>`.


.. py:function:: calc_yearsum(data_input, dim='time', **kwargs)

   Calculate yearly sum.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{sum} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating sum on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.sum <numpy:numpy.sum>`, :py:func:`dask.array.sum <dask:dask.array.sum>`, 
       :py:meth:`xarray.DataArray.sum <xarray:xarray.DataArray.sum>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.sum <xarray:xarray.core.groupby.DataArrayGroupBy.sum>`.


.. py:function:: calc_yearstd(data_input, dim='time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{std} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating std on this object's data. 
       These could include dask-specific kwargs like split_every.

   .. note::
       The parameter `ddof` is `Delta Degrees of Freedom`: the divisor used in the calculation is `N - ddof`, 
       where `N` represents the number of elements. If the data needs to be Normalize by `(n-1)`, then `ddof=1`.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.std <numpy:numpy.std>`, :py:func:`dask.array.std <dask:dask.array.std>`, 
       :py:meth:`xarray.DataArray.std <xarray:xarray.DataArray.std>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.std <xarray:xarray.core.groupby.DataArrayGroupBy.std>`.


.. py:function:: calc_yearvar(data_input, dim='time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{var} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating var on this object's data. 
       These could include dask-specific kwargs like split_every.

   .. note::
       The parameter `ddof` is `Delta Degrees of Freedom`: the divisor used in the calculation is `N - ddof`, 
       where `N` represents the number of elements. If the data needs to be Normalize by `(n-1)`, then `ddof=1`.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.var <numpy:numpy.var>`, :py:func:`dask.array.var <dask:dask.array.var>`, 
       :py:meth:`xarray.DataArray.var <xarray:xarray.DataArray.var>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.var <xarray:xarray.core.groupby.DataArrayGroupBy.var>`.


.. py:function:: calc_yearmax(data_input, dim='time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{max} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating max on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.maximum <numpy:numpy.maximum>`, :py:func:`dask.array.max <dask:dask.array.max>`, 
       :py:meth:`xarray.DataArray.max <xarray:xarray.DataArray.max>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.max <xarray:xarray.core.groupby.DataArrayGroupBy.max>`.


.. py:function:: calc_yearmin(data_input, dim='time', **kwargs)

   Calculate yearly standard deviation.

   For every adjacent sequence :math:`t_1, ..., t_n` of timesteps of the same year it is:

   .. math::
       o(t, x) = \mathrm{min} \left \lbrace i(t', x), t_1 < t' \leqslant t_n \right\rbrace

   .. tip::
       This function uses :py:meth:`xarray.DataArray.groupby<xarray:xarray.DataArray.groupby>` to implement the calculation.
       To substantially improve the performance of GroupBy operations, particularly with dask `install the flox package <https://flox.readthedocs.io/>`_. 
       flox `extends Xarray's in-built GroupBy capabilities <https://flox.readthedocs.io/en/latest/xarray.html>`_ by allowing 
       grouping by multiple variables, and lazy grouping by dask arrays. If installed, Xarray will automatically use flox by default.

   Parameters
   ----------
   data_input : :py:class:`xarray.DataArray<xarray.DataArray>`
       :py:class:`xarray.DataArray<xarray.DataArray>` to be calculated.
   dim : str
       Dimension(s) over which to apply extracting. By default extracting is applied over the `time` dimension.
   **kwargs:
       Additional keyword arguments passed on to the appropriate array function for calculating min on this object's data. 
       These could include dask-specific kwargs like split_every.

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>` with `time` dimension type of :py:class:`numpy.datetime64<numpy.datetime64>`.

   .. seealso::
       :py:func:`numpy.minimum <numpy:numpy.minimum>`, :py:func:`dask.array.min <dask:dask.array.min>`, 
       :py:meth:`xarray.DataArray.min <xarray:xarray.DataArray.min>`, 
       :py:meth:`xarray.core.groupby.DataArrayGroupBy.min <xarray:xarray.core.groupby.DataArrayGroupBy.min>`.


