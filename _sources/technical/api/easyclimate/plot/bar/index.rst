easyclimate.plot.bar
====================

.. py:module:: easyclimate.plot.bar

.. autoapi-nested-parse::

   Bar plot for xarray dataset



Functions
---------

.. autoapisummary::

   easyclimate.plot.bar.bar_plot_with_threshold


Module Contents
---------------

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


