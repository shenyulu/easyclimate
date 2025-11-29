easyclimate.filter.gaussian_filter
==================================

.. py:module:: easyclimate.filter.gaussian_filter

.. autoapi-nested-parse::

   Gaussian filter



Functions
---------

.. autoapisummary::

   easyclimate.filter.gaussian_filter.calc_gaussian_filter


Module Contents
---------------

.. py:function:: calc_gaussian_filter(data: xarray.DataArray, window_length: float, sigma: float = None, dim: str = 'time', keep_attrs: bool = False) -> xarray.DataArray

   Apply a Gaussian filter to data along a specified dimension.

   Parameters
   ----------
   da : :py:class:`xarray.DataArray<xarray.DataArray>`
       Input data array
   window_length: :py:class:`int <int>`, optional
       The window width.
   sigma : :py:class:`float <float>`, optional
       Standard deviation for Gaussian kernel. If None, calculated as :math:`\mathrm{window_length} / \sqrt{8 log2}`.
   dim : :py:class:`str <str>`, optional
       Dimension along which to filter (default: 'time')
   keep_attrs : :py:class:`bool <bool>`, optional
       Whether to preserve attributes (default: False)

   Returns
   -------
   :py:class:`xarray.DataArray<xarray.DataArray>`
       Smoothed data with the same dimensions as input.

   .. minigallery::
       :add-heading: Example(s) related to the function

       ./dynamic_docs/plot_ao_index.py


