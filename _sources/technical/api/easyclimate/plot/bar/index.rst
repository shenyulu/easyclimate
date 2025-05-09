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

.. py:function:: bar_plot_with_threshold(da: xarray.DataArray, threshold: float = 0, pos_color: str = 'red', neg_color: str = 'blue', ax=None, **kwargs) -> matplotlib.container.BarContainer

   Plot a bar chart for a 1D xarray.DataArray with bars colored based on a threshold value.

   Parameters:
   -----------
   da : xarray.DataArray
       1-dimensional data array to plot
   threshold : float, optional
       Threshold value for color separation (default: 0)
   pos_color : str, optional
       Color for bars ≥ threshold (default: 'red')
   neg_color : str, optional
       Color for bars < threshold (default: 'blue')
   ax : matplotlib axes, optional
       Axes object to plot on (uses current axes if None)
   **kwargs :
       Additional arguments passed to plt.bar

   Returns:
   --------
   matplotlib.container.BarContainer
       The bar plot object

   Raises:
   -------
   ValueError
       If input DataArray is not 1-dimensional


